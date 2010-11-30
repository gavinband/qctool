#include <iostream>
#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/GenFileSNPDataSink.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"

namespace genfile {

	std::auto_ptr< SNPDataSink > SNPDataSink::create(
		std::string const& filename,
		std::string const& free_data
	) {
		return SNPDataSink::create( filename, get_compression_type_indicated_by_filename( filename ), free_data ) ;
	}

	std::auto_ptr< SNPDataSink > SNPDataSink::create(
		std::string const& filename,
		CompressionType compression_type,
		std::string const& free_data
	) {
		if( filename_indicates_bgen_format( filename )) {
			if( compression_type == e_GzipCompression ) {
				return std::auto_ptr< SNPDataSink >( new ZippedBGenFileSNPDataSink( filename, free_data )) ;
			}
			else {
				return std::auto_ptr< SNPDataSink >( new BGenFileSNPDataSink( filename, free_data, bgen::e_CompressedSNPBlocks )) ;
			}
		}
		else {
			return std::auto_ptr< SNPDataSink >( new GenFileSNPDataSink( filename, get_chromosome_indicated_by_filename( filename ), compression_type )) ;
		}
	}
	
	SNPDataSink::SNPDataSink():
		m_number_of_samples(0u), m_number_of_snps_written(0u)
	{}

	SNPDataSink::~SNPDataSink()
	{}
	
	SNPDataSink& SNPDataSink::write_snp(
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
		if( m_number_of_samples == 0 ) {
			m_number_of_samples = number_of_samples ;
		}
		else {
			assert( number_of_samples == m_number_of_samples ) ;
		}
		write_snp_impl( number_of_samples, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability ) ;
		if( *this ) {
			++m_number_of_snps_written ;
		}
		return *this ;
	}
}
