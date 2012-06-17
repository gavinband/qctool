
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/GenFileSNPDataSink.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"
#include "genfile/VCFFormatSNPDataSink.hpp"

namespace genfile {

	SNPDataSink::UniquePtr SNPDataSink::create(
		std::string const& filename,
		Metadata const& metadata
	) {
		return SNPDataSink::create_impl( filename, get_compression_type_indicated_by_filename( filename ), metadata ) ;
	}

	SNPDataSink::UniquePtr SNPDataSink::create_impl(
		std::string const& filename,
		CompressionType compression_type,
		Metadata const& metadata
	) {
		std::pair< std::string, std::string > d = uniformise( filename ) ;
		if( d.first == "bgen" ) {
			return SNPDataSink::UniquePtr( new BGenFileSNPDataSink( filename, metadata )) ;
		}
		else if( d.first == "vcf" ) {
			return SNPDataSink::UniquePtr( new VCFFormatSNPDataSink( filename )) ;
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
		SNPIdentifyingData const& snp,
		GenotypeProbabilityGetter const& get_AA_probability,
		GenotypeProbabilityGetter const& get_AB_probability,
		GenotypeProbabilityGetter const& get_BB_probability,
		Info const& info
	) {
		return write_snp(
			number_of_samples,
			snp.get_SNPID(),
			snp.get_rsid(),
			snp.get_position().chromosome(),
			snp.get_position().position(),
			snp.get_first_allele(),
			snp.get_second_allele(),
			get_AA_probability,
			get_AB_probability,
			get_BB_probability,
			info
		) ;
	}

	SNPDataSink& SNPDataSink::write_snp(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		Chromosome chromosome,
		uint32_t SNP_position,
		std::string first_allele,
		std::string second_allele,
		GenotypeProbabilityGetter const& get_AA_probability,
		GenotypeProbabilityGetter const& get_AB_probability,
		GenotypeProbabilityGetter const& get_BB_probability,
		Info const& info
	) {
		if( m_number_of_samples == 0 ) {
			m_number_of_samples = number_of_samples ;
		}
		else {
			assert( number_of_samples == m_number_of_samples ) ;
		}
		write_snp_impl( number_of_samples, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability, info ) ;
		if( *this ) {
			++m_number_of_snps_written ;
		}
		return *this ;
	}

	SNPDataSink& SNPDataSink::write_variant_data(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		if( m_number_of_samples == 0 ) {
			m_number_of_samples = data_reader.get_number_of_samples() ;
		}
		else {
			assert( data_reader.get_number_of_samples() == m_number_of_samples ) ;
		}
		write_variant_data_impl( id_data, data_reader, info ) ;
		if( *this ) {
			++m_number_of_snps_written ;
		}
		return *this ;
	}
	
	void SNPDataSink::write_variant_data_impl(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		throw OperationUnsupportedError( "genfile::SNPDataSink::write_variant_data_impl()", "call", get_spec() ) ;
	}
	
	void SNPDataSink::set_sample_names( SampleNameGetter ) {}
}
