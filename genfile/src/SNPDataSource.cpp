#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/GenFileSNPDataSource.hpp"
#include "genfile/BGenFileSNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/VCFFormatSNPDataSource.hpp"
#include "genfile/get_set.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	std::auto_ptr< SNPDataSource > SNPDataSource::create( std::string const& filename, Chromosome chromosome_hint )
	{
		return SNPDataSource::create( filename, chromosome_hint, get_compression_type_indicated_by_filename( filename )) ;
	}

	std::auto_ptr< SNPDataSource > SNPDataSource::create(
		std::string const& filename,
		Chromosome chromosome_hint,
		CompressionType compression_type
	) {
		std::pair< std::string, std::string > uf = uniformise( filename ) ;
		if( uf.first == "bgen" ) {
			return std::auto_ptr< SNPDataSource >( new BGenFileSNPDataSource( uf.second, compression_type )) ;
		}
		else if( uf.first == "vcf" ) {
			try {
				{
					// test if there's a readable index file.
					std::auto_ptr< std::istream > index_filename = open_text_file_for_input( uf.second + ".index" ) ;
				}
				return SNPDataSource::UniquePtr( new VCFFormatSNPDataSource( uf.second, uf.second + ".index", "GT" )) ;
			}
			catch( ResourceNotOpenedError const& e ) {
				return SNPDataSource::UniquePtr( new VCFFormatSNPDataSource( uf.second, "GT" )) ;
			}
		}
		else if( uf.first == "gen" ) {
			return std::auto_ptr< SNPDataSource >( new GenFileSNPDataSource( uf.second, chromosome_hint, compression_type )) ;
		}
		else {
			// assume GEN format.
			return std::auto_ptr< SNPDataSource >( new GenFileSNPDataSource( uf.second, chromosome_hint, compression_type )) ;
		}
	}

	std::auto_ptr< SNPDataSource > SNPDataSource::create_chain(
		std::vector< wildcard::FilenameMatch > const& matches,
		NotifyProgress notify_progress
	) {
		std::auto_ptr< SNPDataSourceChain > ptr = SNPDataSourceChain::create( matches, notify_progress ) ;
		return std::auto_ptr< SNPDataSource >( ptr.release() ) ;
	}

	SNPDataSource::SNPDataSource()
	: m_number_of_snps_read(0), m_state( e_HaveNotReadIdentifyingData ), m_source_reset_callback(0)
	{}

	SNPDataSource::~SNPDataSource() {}

	void SNPDataSource::reset_to_start() {
		reset_to_start_impl() ;
		m_number_of_snps_read = 0 ;
		m_state = e_HaveNotReadIdentifyingData ;
		if( m_source_reset_callback ) {
			m_source_reset_callback() ;
		}
	}

	void SNPDataSource::set_source_reset_callback( source_reset_callback_t callback ) {
		m_source_reset_callback = callback ;
	}

	SNPDataSource& SNPDataSource::read_snp(
		IntegerSetter set_number_of_samples,
		StringSetter set_SNPID,
		StringSetter set_RSID,
		ChromosomeSetter set_chromosome,
		SNPPositionSetter set_SNP_position,
		AlleleSetter set_allele1,
		AlleleSetter set_allele2,
		GenotypeProbabilitySetter set_genotype_probabilities
	)
	{
		uint32_t this_number_of_samples ;
		get_snp_identifying_data( set_value( this_number_of_samples ), set_SNPID, set_RSID, set_chromosome, set_SNP_position, set_allele1, set_allele2 ) ;
		if( *this ) {
			assert( this_number_of_samples == number_of_samples() ) ;
			read_snp_probability_data( set_genotype_probabilities ) ;
			set_number_of_samples( this_number_of_samples ) ;
		}
		return *this ;
	}
	
	bool SNPDataSource::get_next_snp_with_specified_position(
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2,
		Chromosome specified_chromosome,
		uint32_t specified_SNP_position
	) {
		return get_next_snp_with_position_in_range(
			set_number_of_samples,
			set_SNPID,
			set_RSID,
			set_chromosome,
			set_SNP_position,
			set_allele1,
			set_allele2,
			specified_chromosome,
			specified_chromosome,
			specified_SNP_position,
			specified_SNP_position
		) ;
	}

	bool SNPDataSource::get_next_snp_with_specified_position(
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2,
		GenomePosition specified_position
	) {
		return get_next_snp_with_position_in_range(
			set_number_of_samples,
			set_SNPID,
			set_RSID,
			set_chromosome,
			set_SNP_position,
			set_allele1,
			set_allele2,
			specified_position.chromosome(),
			specified_position.chromosome(),
			specified_position.position(),
			specified_position.position()
		) ;
	}

	bool SNPDataSource::get_next_snp_with_position_in_range(
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2,
		Chromosome chromosome_lower_bound,
		Chromosome chromosome_upper_bound,
		uint32_t position_lower_bound,
		uint32_t position_upper_bound
	) {
		assert( chromosome_lower_bound <= chromosome_upper_bound ) ;
		assert( position_lower_bound <= position_upper_bound ) ;

		uint32_t this_number_of_samples = number_of_samples() ;
		std::string SNPID, RSID ;
		unsigned char chromosome ;
		uint32_t SNP_position ;
		char first_allele, second_allele ;

		while( *this )  {
			get_snp_identifying_data( set_value( this_number_of_samples ), set_value( SNPID ), set_value( RSID ), set_value( chromosome ), set_value( SNP_position ), set_value( first_allele ), set_value( second_allele ) ) ;

			if( *this ) {
				assert( this_number_of_samples == number_of_samples() ) ;
				if( chromosome < chromosome_lower_bound || ((chromosome <= chromosome_upper_bound) && (SNP_position < position_lower_bound ))) {
					ignore_snp_probability_data() ;
				}
				else if( chromosome > chromosome_upper_bound || SNP_position > position_upper_bound ) {
					return false ;
				}
				else {
					set_SNPID( SNPID ) ;
					set_RSID( RSID ) ;
					set_chromosome( Chromosome( chromosome )) ;
					set_SNP_position( SNP_position ) ;
					set_allele1( first_allele ) ;
					set_allele2( second_allele ) ;
					set_number_of_samples( this_number_of_samples ) ;

					return *this ;
				}
			}
		}

		return false ;	
	}

	SNPDataSource& SNPDataSource::get_snp_identifying_data(
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		get_snp_identifying_data_impl(
			set_number_of_samples,
			set_SNPID,
			set_RSID,
			set_chromosome,
			set_SNP_position,
			set_allele1,
			set_allele2
		) ;
		if( *this ) {
			m_state = e_HaveReadIdentifyingData ;
		}
		return *this ;
	}

	SNPDataSource& SNPDataSource::get_snp_identifying_data(
		SNPIdentifyingData& snp
	) {
		return get_snp_identifying_data(
			ignore(),
			set_value( snp.SNPID() ),
			set_value( snp.rsid() ),
			set_value( snp.position().chromosome() ),
			set_value( snp.position().position() ),
			set_value( snp.first_allele() ),
			set_value( snp.second_allele() )
		) ;
	}

	SNPDataSource& SNPDataSource::read_snp_probability_data(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		assert( m_state == e_HaveReadIdentifyingData ) ;
		read_snp_probability_data_impl( set_genotype_probabilities ) ;
		if( *this ) {
			m_state = e_HaveNotReadIdentifyingData ;
			++m_number_of_snps_read ;
		} else {
			throw MalformedInputError( get_source_spec(), m_number_of_snps_read ) ;
		}
		return *this ;
	}
	
	SNPDataSource& SNPDataSource::ignore_snp_probability_data() {
		assert( m_state == e_HaveReadIdentifyingData ) ;
		ignore_snp_probability_data_impl() ;
		if( *this ) {
			m_state = e_HaveNotReadIdentifyingData ;
			++m_number_of_snps_read ;
		}
		return *this ;
	}

	std::string SNPDataSource::get_summary( std::string const& prefix, std::size_t width ) const {
		std::ostringstream ostr ;
		ostr << prefix << std::setw( width ) << "Spec: " << get_source_spec() << "\n" ;
		ostr << prefix << std::setw( width ) << "Number of samples: " <<  number_of_samples() << "\n" ;
		ostr << prefix << std::setw( width ) << "Number of SNPs: " <<  number_of_samples() << "\n" ;
		return ostr.str() ;
	}	

	void IdentifyingDataCachingSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		if( state() != e_HaveReadIdentifyingData ) {
			read_snp_identifying_data_impl(
				&m_cached_number_of_samples,
				&m_cached_SNPID,
				&m_cached_RSID,
				&m_cached_chromosome,
				&m_cached_SNP_position,
				&m_cached_allele1,
				&m_cached_allele2
			) ;
		}

		set_number_of_samples( m_cached_number_of_samples ) ;
		set_SNPID( m_cached_SNPID ) ;
		set_RSID( m_cached_RSID ) ;
		set_chromosome( m_cached_chromosome ) ;
		set_SNP_position( m_cached_SNP_position ) ;
		set_allele1( m_cached_allele1 ) ;
		set_allele2( m_cached_allele2 ) ;
	}
	
}

