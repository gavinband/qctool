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

	SNPDataSink::UniquePtr SNPDataSink::create(
		std::string const& filename,
		std::string const& free_data
	) {
		return SNPDataSink::create_impl( filename, get_compression_type_indicated_by_filename( filename ), free_data ) ;
	}

	SNPDataSink::UniquePtr SNPDataSink::create_impl(
		std::string const& filename,
		CompressionType compression_type,
		std::string const& free_data
	) {
		std::pair< std::string, std::string > d = uniformise( filename ) ;
		if( d.first == "bgen" ) {
			if( compression_type == "gzip_compression" ) {
				return SNPDataSink::UniquePtr( new ZippedBGenFileSNPDataSink( filename, free_data )) ;
			}
			else {
				return SNPDataSink::UniquePtr( new BGenFileSNPDataSink( filename, free_data, bgen::e_CompressedSNPBlocks )) ;
			}
		}
		else {
			return SNPDataSink::UniquePtr( new GenFileSNPDataSink( filename, get_chromosome_indicated_by_filename( filename ), compression_type )) ;
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

	SNPDataSink& SNPDataSink::write_variant_data(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader
	) {
		assert(0) ;
		/*
		if( m_number_of_samples == 0 ) {
			m_number_of_samples = id_data.get_number_of_samples() ;
		}
		else {
			assert( id_data.get_number_of_samples() == m_number_of_samples ) ;
		}
		write_variant_data( id_data, data_reader ) ;
		if( *this ) {
			++m_number_of_snps_written ;
		}
		*/
		return *this ;
	}
	
	void SNPDataSink::write_variant_data_impl(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader
	) {
		assert(0) ;
	}
}
