#ifndef GENFILESNPDATASINK_HPP
#define GENFILESNPDATASINK_HPP

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/SNPDataSource.hpp"

namespace genfile {
	
	// This class represents a SNPDataSink which writes its data
	// to a plain GEN file.
	class GenFileSNPDataSink: public SNPDataSink
	{
	public:
		GenFileSNPDataSink( std::string const& filename, Chromosome chromosome )
			: m_filename( filename )
		{
			setup( filename, get_compression_type_indicated_by_filename( filename )) ;
		}

		GenFileSNPDataSink( std::string const& filename, Chromosome chromosome, CompressionType compression_type )
			: m_filename( filename )
		{
			setup( filename, compression_type ) ;
		}

	private:
		
		// Methods required by SNPDataSink
		operator bool() const { return *m_stream_ptr ; }
		
		void write_snp_impl(
			uint32_t number_of_samples,
			std::string SNPID,
			std::string RSID,
			Chromosome chromosome,
			uint32_t SNP_position,
			char first_allele,
			char second_allele,
			GenotypeProbabilityGetter const& get_AA_probability,
			GenotypeProbabilityGetter const& get_AB_probability,
			GenotypeProbabilityGetter const& get_BB_probability
		) {
			if( number_of_snps_written() == 0 ) {
				m_chromosome = chromosome ;
			}
			else if( m_chromosome != chromosome ) {
				throw BadArgumentError( "GenFileSNPDataSink::write_snp_impl()", "chromosome=" + string_utils::to_string( chromosome ) ) ;
			}

			gen::write_snp_block( *m_stream_ptr, number_of_samples, SNPID, RSID, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability ) ;
		} ;

	private:
		// This class's implementation
		void setup( std::string const& filename, CompressionType compression_type ) {
			m_stream_ptr = open_text_file_for_output( filename, compression_type ) ;
			if( !(*m_stream_ptr)) {
				throw ResourceNotOpenedError( m_filename ) ;
			}
			assert( *m_stream_ptr ) ;
		}

		std::string m_filename ;
		std::auto_ptr< std::ostream > m_stream_ptr ;
		Chromosome m_chromosome ;
	} ;
}

#endif