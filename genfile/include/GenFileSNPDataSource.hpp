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
	class GenFileSNPDataSource: public IdentifyingDataCachingSNPDataSource
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
		void read_snp_identifying_data_impl( 
			uint32_t*, // number_of_samples is unused.
			std::string* SNPID,
			std::string* RSID,
			uint32_t* SNP_position,
			char* allele1,
			char* allele2
		) {
			gen::impl::read_snp_identifying_data( stream(), SNPID, RSID, SNP_position, allele1, allele2 ) ;
		}

		void read_snp_probability_data_impl(
			uint32_t* number_of_samples,
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) {
			gen::impl::read_snp_probability_data( stream(), set_value( *number_of_samples ), set_genotype_probabilities ) ;
		}

		void ignore_snp_probability_data_impl( uint32_t number_of_samples ) {
			std::string line ;
			std::getline( stream(), line ) ;
		}

	private:
		std::string m_filename ;
		unsigned int m_number_of_samples, m_total_number_of_snps ;
		std::auto_ptr< std::istream > m_stream_ptr ;

		void setup( std::string const& filename, CompressionType compression_type ) {
			m_stream_ptr = open_text_file_for_input( filename, compression_type ) ;
			read_header_data() ;				
			// That will have consumed the file, so re-open it.
			m_stream_ptr = open_text_file_for_input( filename, compression_type ) ;
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