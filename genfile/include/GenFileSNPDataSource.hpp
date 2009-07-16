#ifndef GENFILESNPDATAPROVIDER_HPP
#define GENFILESNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "gen.hpp"
#include "SNPDataSource.hpp"

namespace genfile {
	// This class represents a SNPDataSource which reads its data
	// from a plain GEN file.
	class GenFileSNPDataSource: public SNPDataSource
	{
	public:
		GenFileSNPDataSource( std::string const& filename )
			: m_filename( filename )
		{
			setup( filename, get_compression_type_indicated_by_filename( filename )) ;
		}

		GenFileSNPDataSource( std::string const& filename, CompressionType compression_type )
			: m_filename( filename )
		{
			setup( filename, compression_type ) ;
		}

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const { return m_total_number_of_snps ; }
		
		operator bool() const { return *m_stream_ptr ; }
		std::istream& stream() { return *m_stream_ptr ; }
		std::istream const& stream() const { return *m_stream_ptr ; }

	private:
		
		void read_snp_impl(
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2,
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) {
			gen::read_snp_block( stream(), set_number_of_samples, set_SNPID, set_RSID, set_SNP_position, set_allele1, set_allele2, set_genotype_probabilities ) ;
		} ;

		void read_next_matching_snp_impl(
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2,
			GenotypeProbabilitySetter const& set_genotype_probabilities,
			SNPMatcher const& snp_matcher
		) {
			std::string SNPID, RSID ;
			uint32_t SNP_position ;
			char first_allele, second_allele ;

			do {
				gen::impl::read_snp_identifying_data( stream(), &SNPID, &RSID, &SNP_position, &first_allele, &second_allele ) ;
				if( snp_matcher( SNPID, RSID, SNP_position, first_allele, second_allele )) {
					set_SNPID( SNPID ) ;
					set_RSID( RSID ) ;
					set_SNP_position( SNP_position ) ;
					set_allele1( first_allele ) ;
					set_allele2( second_allele ) ;
					gen::impl::read_snp_probability_data( stream(), set_number_of_samples, set_genotype_probabilities ) ;

					return ;
				}
				else {
					// discard the rest of the line.
					std::string line ;
					std::getline( stream(), line ) ;
				}
			}
			while( *this ) ;
		}

	private:
		std::string m_filename ;
		unsigned int m_number_of_samples, m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;

		void setup( std::string const& filename, CompressionType compression_type ) {
			m_stream_ptr = open_text_file_for_input( filename, compression_type ) ;
			if( !(*m_stream_ptr)) {
				throw FileNotOpenedError() ;
			}
			read_header_data() ;
			// That will have consumed the file, so re-open it.
			m_stream_ptr = open_text_file_for_input( filename, compression_type ) ;
			assert( *m_stream_ptr ) ;
		}

		void read_header_data() {
			gen::read_header_information(
				*m_stream_ptr,
				set_value( m_total_number_of_snps ),
				set_value( m_number_of_samples ),
				ignore()
			) ;
		}
	} ;
}

#endif