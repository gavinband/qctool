#include <vector>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceRack.hpp"

namespace genfile {
	
	std::auto_ptr< SNPDataSourceRack > SNPDataSourceRack::create( std::vector< wildcard::FilenameMatch > const& filenames ) {
		std::auto_ptr< SNPDataSourceRack > rack( new SNPDataSourceRack() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			rack->add_source( SNPDataSource::create( filenames[i].filename(), filenames[i].match() )) ;
		}
		return rack ;
	}
	
	SNPDataSourceRack::SNPDataSourceRack()
		: m_number_of_samples(0),
		  m_read_past_end( false )
	{
	}

	SNPDataSourceRack::~SNPDataSourceRack() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			delete m_sources[i] ;
		}
	}

	void SNPDataSourceRack::add_source( std::auto_ptr< SNPDataSource > source ) {
		m_sources.push_back( source.release() ) ;
		m_number_of_samples += m_sources.back()->number_of_samples() ;
		std::vector< GenomePosition > intersected_positions = get_intersected_snp_positions( m_positions_of_included_snps, *m_sources.back() ) ;
		m_positions_of_included_snps.swap( intersected_positions ) ;
		m_sources.back()->reset_to_start() ;
	}

	std::vector< GenomePosition > SNPDataSourceRack::get_intersected_snp_positions(
		std::vector< GenomePosition > const& snp_positions,
		SNPDataSource& source
	) const {
		std::vector< GenomePosition > source_positions = get_source_snp_positions( source ) ;
		std::vector< GenomePosition > intersected_positions ;
		if( m_sources.size() == 1 ) {
			intersected_positions.swap( source_positions ) ;
		}
		else {
			intersected_positions.reserve( source_positions.size() ) ;
			std::set_intersection(
				snp_positions.begin(), snp_positions.end(),
				source_positions.begin(), source_positions.end(),
				std::back_inserter( intersected_positions )
			) ;
		}

		return intersected_positions ;
	}

	std::vector< GenomePosition > SNPDataSourceRack::get_source_snp_positions( SNPDataSource& source ) const {
		std::vector< GenomePosition > result ;
		result.reserve( source.total_number_of_snps() ) ;
		for(
			GenomePosition position ;
			source.get_snp_identifying_data(
				genfile::ignore(),
				genfile::ignore(),
				genfile::ignore(),
				genfile::set_value( position.chromosome() ),
				genfile::set_value( position.position() ),
				genfile::ignore(),
				genfile::ignore()
			) ;
			source.ignore_snp_probability_data()
		) {
			result.push_back( position ) ;
		}
		return result ;
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
		return m_positions_of_included_snps.size() ;
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
		
		if( number_of_snps_read() == m_positions_of_included_snps.size() ) {
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
				m_positions_of_included_snps[ number_of_snps_read() ].chromosome(),
				m_positions_of_included_snps[ number_of_snps_read() ].position()
			)
		) {
			throw MissingSNPError( 0, m_positions_of_included_snps[ number_of_snps_read() ] ) ;
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
		std::string this_source_SNPID, this_source_RSID ;
		Chromosome this_source_chromosome ;
		uint32_t this_source_SNP_position ;
		char this_source_allele1, this_source_allele2 ;
		
		if(
			!m_sources[source_i]->get_next_snp_with_specified_position(
				ignore(),
				set_value( this_source_SNPID ),
				set_value( this_source_RSID ),
				set_value( this_source_chromosome ),
				set_value( this_source_SNP_position ),
				set_value( this_source_allele1 ),
				set_value( this_source_allele2 ),
				chromosome,
				SNP_position
			)
		) {
			throw MissingSNPError( source_i, GenomePosition( chromosome, SNP_position )) ;
		}
		else if( this_source_SNPID != SNPID
				|| this_source_RSID != RSID
				|| this_source_allele1 != allele1
				|| this_source_allele2 != allele2
		) {
				throw SNPMismatchError( source_i, GenomePosition( chromosome, SNP_position )) ;
		}
		else {
			return ;
		}
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