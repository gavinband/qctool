#include <iostream>
#include <string>

#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"
#include "GenFileSNPDataSource.hpp"
#include "BGenFileSNPDataSource.hpp"
#include "SNPDataSourceChain.hpp"

namespace genfile {
	std::auto_ptr< SNPDataSource > SNPDataSource::create( std::string const& filename ) {
		return SNPDataSource::create( filename, get_compression_type_indicated_by_filename( filename )) ;
	}

	std::auto_ptr< SNPDataSource > SNPDataSource::create( std::string const& filename, CompressionType compression_type ) {
		if( filename_indicates_bgen_format( filename )) {
			return std::auto_ptr< SNPDataSource >( new BGenFileSNPDataSource( filename, compression_type )) ;
		}
		else {
			return std::auto_ptr< SNPDataSource >( new GenFileSNPDataSource( filename, compression_type )) ;
		}
	}

	std::auto_ptr< SNPDataSource > SNPDataSource::create( std::vector< std::string > const& filenames ) {
		std::auto_ptr< SNPDataSourceChain > chain( new SNPDataSourceChain() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			chain->add_source( SNPDataSource::create( filenames[i] )) ;
		}
		return std::auto_ptr< SNPDataSource >( chain.release() ) ;
	}

	SNPDataSource::SNPDataSource()
	: m_number_of_snps_read(0), m_state( e_HaveNotReadIdentifyingData )
	{}

	SNPDataSource::~SNPDataSource() {} ;

	SNPDataSource& SNPDataSource::read_snp(
		IntegerSetter set_number_of_samples,
		StringSetter set_SNPID,
		StringSetter set_RSID,
		SNPPositionSetter set_SNP_position,
		AlleleSetter set_allele1,
		AlleleSetter set_allele2,
		GenotypeProbabilitySetter set_genotype_probabilities
	)
	{
		uint32_t number_of_samples ;
		get_snp_identifying_data( set_value( number_of_samples ), set_SNPID, set_RSID, set_SNP_position, set_allele1, set_allele2 ) ;
		if( *this ) {
			read_snp_probability_data( &number_of_samples, set_genotype_probabilities ) ;
			set_number_of_samples( number_of_samples ) ;
		}
		return *this ;
	}
	
	bool SNPDataSource::read_next_snp_with_specified_position(
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2,
		GenotypeProbabilitySetter const& set_genotype_probabilities,
		uint32_t specified_SNP_position
	) {
		return read_next_snp_with_position_in_range(
			set_number_of_samples,
			set_SNPID,
			set_RSID,
			set_SNP_position,
			set_allele1,
			set_allele2,
			set_genotype_probabilities,
			specified_SNP_position,
			specified_SNP_position
		) ;
	}

	bool SNPDataSource::read_next_snp_with_position_in_range(
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2,
		GenotypeProbabilitySetter const& set_genotype_probabilities,
		uint32_t position_lower_bound,
		uint32_t position_upper_bound
	) {
		assert( position_lower_bound <= position_upper_bound ) ;

		uint32_t number_of_samples ;
		std::string SNPID, RSID ;
		uint32_t SNP_position ;
		char first_allele, second_allele ;

		while( *this )  {
			get_snp_identifying_data( set_value( number_of_samples ), set_value( SNPID ), set_value( RSID ), set_value( SNP_position ), set_value( first_allele ), set_value( second_allele ) ) ;

			if( *this ) {
				if( SNP_position < position_lower_bound ) {
					ignore_snp_probability_data() ;
				}
				else if( SNP_position > position_upper_bound ) {
					return false ;
				}
				else {
					set_SNPID( SNPID ) ;
					set_RSID( RSID ) ;
					set_SNP_position( SNP_position ) ;
					set_allele1( first_allele ) ;
					set_allele2( second_allele ) ;

					read_snp_probability_data( &number_of_samples, set_genotype_probabilities ) ;
					set_number_of_samples( number_of_samples ) ;
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
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		get_snp_identifying_data_impl(
			set_number_of_samples,
			set_SNPID,
			set_RSID,
			set_SNP_position,
			set_allele1,
			set_allele2
		) ;
		if( *this ) {
			m_state = e_HaveReadIdentifyingData ;
		}
		return *this ;
	}

	SNPDataSource& SNPDataSource::read_snp_probability_data(
		uint32_t* number_of_samples,
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		assert( m_state == e_HaveReadIdentifyingData ) ;
		read_snp_probability_data_impl( number_of_samples, set_genotype_probabilities ) ;
		if( *this ) {
			m_state = e_HaveNotReadIdentifyingData ;
			++m_number_of_snps_read ;
		}
		return *this ;
	} ;

	SNPDataSource& SNPDataSource::ignore_snp_probability_data() {
		assert( m_state == e_HaveReadIdentifyingData ) ;
		ignore_snp_probability_data_impl() ;
		if( *this ) {
			m_state = e_HaveNotReadIdentifyingData ;
			++m_number_of_snps_read ;
		}
		return *this ;
	}

	void IdentifyingDataCachingSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		if( state() != e_HaveReadIdentifyingData ) {
			read_snp_identifying_data_impl(
				&m_cached_number_of_samples,
				&m_cached_SNPID,
				&m_cached_RSID,
				&m_cached_SNP_position,
				&m_cached_allele1,
				&m_cached_allele2
			) ;
		}

		set_number_of_samples( m_cached_number_of_samples ) ;
		set_SNPID( m_cached_SNPID ) ;
		set_RSID( m_cached_RSID ) ;
		set_SNP_position( m_cached_SNP_position ) ;
		set_allele1( m_cached_allele1 ) ;
		set_allele2( m_cached_allele2 ) ;
	}
}

