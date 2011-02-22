#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceRack.hpp"
#include "genfile/get_set.hpp"
#include "genfile/get_list_of_snps_in_source.hpp"

namespace genfile {
	
	std::auto_ptr< SNPDataSourceRack > SNPDataSourceRack::create( std::vector< wildcard::FilenameMatch > const& filenames ) {
		std::auto_ptr< SNPDataSourceRack > rack( new SNPDataSourceRack() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			rack->add_source( SNPDataSource::create( filenames[i].filename(), filenames[i].match() )) ;
		}
		return rack ;
	}

	std::auto_ptr< SNPDataSourceRack > SNPDataSourceRack::create(
		std::vector< wildcard::FilenameMatch > const& filenames,
		std::string const& snp_match_fields
	) {
		std::auto_ptr< SNPDataSourceRack > rack( new SNPDataSourceRack( snp_match_fields ) ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			rack->add_source( SNPDataSource::create( filenames[i].filename(), filenames[i].match() )) ;
		}
		return rack ;
	}
	
	SNPDataSourceRack::SNPDataSourceRack()
		: m_number_of_samples(0),
		  m_read_past_end( false ),
		  m_comparator( "position,rsid,SNPID,alleles" )
	{
	}

	SNPDataSourceRack::SNPDataSourceRack( std::string const& snp_match_fields )
		: m_number_of_samples(0),
		  m_read_past_end( false ),
		  m_comparator( snp_match_fields )
	{
		assert( snp_match_fields.find( "position" ) != std::string::npos ) ;
	}

	SNPDataSourceRack::~SNPDataSourceRack() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			delete m_sources[i] ;
		}
	}

	void SNPDataSourceRack::add_source( std::auto_ptr< SNPDataSource > source ) {
		std::vector< SNPIdentifyingData > const snps = get_list_of_snps_in_source( *source ) ;
		add_source( source, snps ) ;
	}

	void SNPDataSourceRack::add_source(
		std::auto_ptr< SNPDataSource > source,
		std::vector< SNPIdentifyingData > const& snps
	) {
		m_sources.push_back( source.release() ) ;
		m_number_of_samples += m_sources.back()->number_of_samples() ;
		m_sources.back()->reset_to_start() ;
		
		if( m_sources.size() == 1 ) {
			m_included_snps = snps ;
		}
		else {
			m_included_snps = get_intersected_snps( m_included_snps, snps ) ;
		}
	}

	bool SNPDataSourceRack::check_snps_are_sorted_by_position(
		std::vector< SNPIdentifyingData > const& snps
	) {
		for( std::size_t i = 1; i < snps.size(); ++i ) {
			if( snps[i].get_position() < snps[i-1].get_position() ) {
				return false ;
			}
		}
		return true ;
	}

	std::vector< SNPIdentifyingData > SNPDataSourceRack::get_intersected_snps(
		std::vector< SNPIdentifyingData > const& snps1,
		std::vector< SNPIdentifyingData > const& snps2
	) const {
		if( !check_snps_are_sorted_by_position( snps1 )) {
			throw BadArgumentError( "SNPDataSourceRack::get_intersected_snps()", "snps1 should be in nondecreasing order of position" ) ;
		}
		if( !check_snps_are_sorted_by_position( snps2 )) {
			throw BadArgumentError( "SNPDataSourceRack::get_intersected_snps()", "snps2 should be in nondecreasing order of position" ) ;
		}

		//
		// The algorithm here must match that for get_snp_identifying_data_impl.
		// This has the following feature: each list can only be traversed once,
		// forwards.
		// For each snp in the first list we find the first SNP not already considered in
		// the second list which (has the same position and) matches according to our comparator.
		// 
		std::vector< SNPIdentifyingData >::const_iterator
			i1 = snps1.begin(),
			i1_end = snps1.end(),
			i2 = snps2.begin(),
			i2_end = snps2.end() ;

		SNPIdentifyingData::CompareFields position_comparator( "position" ) ;

		std::vector< SNPIdentifyingData > intersected_snps ;

		for( ; i1 != i1_end && i2 != i2_end ; ++i1 ) {
			// Find next SNP with matching position.
			std::vector< SNPIdentifyingData >::const_iterator
				i2_pos = std::lower_bound( i2, i2_end, *i1, position_comparator ) ;
			// Find next SNP with all fields matching, including position.
			for( ; i2_pos != i2_end && !m_comparator.are_equal( *i2_pos, *i1 ) && i2_pos->get_position() == i1->get_position(); ++i2_pos ) ;

			if( i2_pos != i2_end && i2_pos->get_position() == i1->get_position() ) {
				intersected_snps.push_back( *i1 ) ;
				++i2_pos ;
			}
			i2 = i2_pos ;
		}

		return intersected_snps ;
	}

	// Implicit conversion to bool.  Return true if there have been no errors so far.
	SNPDataSourceRack::operator bool() const {
		if( m_read_past_end ) {
			return false ;
		}

		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			if( !m_sources[i]->operator bool() ) {
				return false ;
			}
		}
		return !m_sources.empty() ;
	}
	
	std::string SNPDataSourceRack::get_source_spec() const {
		std::string result = "rack:" ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			result += "," + m_sources[i]->get_source_spec() ;
		}
		return result ;
	}
	
	std::string SNPDataSourceRack::get_summary( std::string const& prefix, std::size_t width ) const {
		std::ostringstream ostr ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			ostr << prefix << std::setw( width ) << "cohort " << (i+1) << ":\n" ;
			ostr << m_sources[i]->get_summary( prefix, width ) ;
			ostr << "\n" ;
		}
		ostr << prefix << std::setw( width ) << "Total all cohorts:" << " " << number_of_samples() << " samples, " << total_number_of_snps() << " overlap SNPs.\n" ;
		return ostr.str() ;
	}

	SNPDataSource& SNPDataSourceRack::get_source( std::size_t index ) const {
		assert( index < m_sources.size() ) ;
		return *m_sources[ index ] ;
	}

	// Return the number of samples represented in the snps in this source.
	unsigned int SNPDataSourceRack::number_of_samples() const {
		return m_number_of_samples ;
	}

	// Return the total number of snps the source contains.
	unsigned int SNPDataSourceRack::total_number_of_snps() const {
		return m_included_snps.size() ;
	}

	void SNPDataSourceRack::reset_to_start_impl() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			m_sources[i]->reset_to_start() ;
		}
		m_read_past_end = false ;
	}

	void SNPDataSourceRack::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		if( m_sources.size() == 0 ) {
			return ;
		}
		std::string SNPID ;
		std::string RSID ;
		Chromosome chromosome ;
		uint32_t SNP_position ;
		char allele1, allele2 ;
		
		if( number_of_snps_read() == m_included_snps.size() ) {
			m_read_past_end = true ;
		}
		else if(
			!m_sources[0]->get_next_snp_with_specified_position(
				ignore(),
				set_value( SNPID ),
				set_value( RSID ),
				set_value( chromosome ),
				set_value( SNP_position ),
				set_value( allele1 ),
				set_value( allele2 ),
				m_included_snps[ number_of_snps_read() ].get_position().chromosome(),
				m_included_snps[ number_of_snps_read() ].get_position().position()
			)
		) {
			throw MissingSNPError( 0, m_included_snps[ number_of_snps_read() ] ) ;
		}
		else {
			if( *this ) {
				for( std::size_t i = 1; i < m_sources.size(); ++i ) {
					move_source_to_snp_matching(
						i,
						SNPID,
						RSID,
						chromosome,
						SNP_position,
						allele1,
						allele2
					) ;
				}
			}
		
			if( *this ) {
				set_number_of_samples( m_number_of_samples ) ;
				set_SNPID( SNPID ) ;
				set_RSID( RSID ) ;
				set_chromosome( chromosome ) ;
				set_SNP_position( SNP_position ) ;
				set_allele1( allele1 ) ;
				set_allele2( allele2 ) ;
			}
		}
	}

	// Read the identifying info of the next SNP in the given source
	// which has the given chromosome and SNP position.
	// If no such SNP is found, throw MissingSNPError.
	// If such a SNP is found but it has different SNPID, RSID, or alleles, throw SNPMismatchError.
	void SNPDataSourceRack::move_source_to_snp_matching(
		std::size_t source_i,
		std::string SNPID,
		std::string RSID,
		Chromosome chromosome,
		uint32_t SNP_position,
		char allele1,
		char allele2
	) {
		assert( source_i > 0 ) ;
		SNPIdentifyingData reference_snp(
			SNPID,
			RSID,
			GenomePosition( chromosome, SNP_position ),
			allele1,
			allele2
		) ;
		SNPIdentifyingData this_snp ;
		
		while( m_sources[source_i]->get_next_snp_with_specified_position(
			ignore(),
			set_value( this_snp.SNPID() ),
			set_value( this_snp.rsid() ),
			set_value( this_snp.position().chromosome() ),
			set_value( this_snp.position().position() ),
			set_value( this_snp.first_allele() ),
			set_value( this_snp.second_allele() ),
			chromosome,
			SNP_position
		) ) {
			if( m_comparator.are_equal( this_snp, reference_snp )) {
				return ;
			}

			m_sources[source_i]->ignore_snp_probability_data() ;
		}
		throw MissingSNPError( source_i, reference_snp ) ;
	}
		

	void SNPDataSourceRack::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		uint32_t number_of_samples_so_far = 0u ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			m_sources[i]->read_snp_probability_data( RackGenotypeProbabilitySetter( set_genotype_probabilities, number_of_samples_so_far )) ;
			number_of_samples_so_far += m_sources[i]->number_of_samples() ;
		}
	}

	void SNPDataSourceRack::ignore_snp_probability_data_impl() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			m_sources[i]->ignore_snp_probability_data() ;
		}
	}
	
	SNPDataSourceRack::RackGenotypeProbabilitySetter::RackGenotypeProbabilitySetter(
		GenotypeProbabilitySetter const& base_setter,
		uint32_t index_of_first_sample
	) 
		: m_base_setter( base_setter ),
	  	  m_index_of_first_sample( index_of_first_sample )
	{}
	
	void SNPDataSourceRack::RackGenotypeProbabilitySetter::operator()(
		std::size_t i, double AA, double AB, double BB
	) const {
		m_base_setter( i + m_index_of_first_sample, AA, AB, BB ) ;
	}
}
