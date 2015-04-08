
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <map>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <Eigen/Core>

#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/GenFileSNPDataSink.hpp"
#include "genfile/GenDosageFileSNPDataSink.hpp"
#include "genfile/GenIntensityFileSNPDataSink.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"
#include "genfile/PedigreeMappingBedFileSNPDataSink.hpp"
#include "genfile/FromSamplesPedigree.hpp"
#include "genfile/ShapeITHaplotypesSNPDataSink.hpp"
#include "genfile/VCFFormatSNPDataSink.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "genfile/get_set_eigen.hpp"
#include "genfile/get_set.hpp"

namespace genfile {

	std::vector< std::string > SNPDataSink::get_file_types() {
		std::vector< std::string > result ;
		result.push_back( "gen" ) ;
		result.push_back( "bgen" ) ;
		result.push_back( "bgen_v12" ) ;
		result.push_back( "vcf" ) ;
		result.push_back( "shapeit_haplotypes" ) ;
		result.push_back( "shapeit" ) ;
		result.push_back( "dosage" ) ;
		result.push_back( "intensity" ) ;
		return result ;
	}

	SNPDataSink::UniquePtr SNPDataSink::create(
		std::string const& filename,
		Metadata const& metadata,
		std::string const& filetype_hint
	) {
		return SNPDataSink::create_impl( filename, get_compression_type_indicated_by_filename( filename ), metadata, filetype_hint ) ;
	}

	SNPDataSink::UniquePtr SNPDataSink::create_impl(
		std::string const& filename,
		CompressionType compression_type,
		Metadata const& metadata,
		std::string const& filetype_hint
	) {
		std::pair< std::string, std::string > d = uniformise( filename ) ;
		if( filetype_hint != "guess" ) {
			d.first = filetype_hint ;
		}
		if( d.first == "bgen" ) {
			return SNPDataSink::UniquePtr( new BGenFileSNPDataSink( filename, metadata )) ;
		}
		else if( d.first == "bgen_v12" ) {
			return SNPDataSink::UniquePtr( new BGenFileSNPDataSink( filename, metadata, "v12", 16 )) ;
		}
		else if( d.first == "vcf" ) {
			return SNPDataSink::UniquePtr( new VCFFormatSNPDataSink( filename )) ;
		}
		else if( d.first == "shapeit_haplotypes" || d.first == "shapeit" ) {
			return SNPDataSink::UniquePtr( new ShapeITHaplotypesSNPDataSink( filename, get_chromosome_indicated_by_filename( filename ), compression_type )) ;
		}
		else if( d.first == "dosage" ) {
			return SNPDataSink::UniquePtr( new GenDosageFileSNPDataSink( filename, get_chromosome_indicated_by_filename( filename ), compression_type )) ;
		}
		else if( d.first == "intensity" ) {
			return SNPDataSink::UniquePtr( new GenIntensityFileSNPDataSink( filename, get_chromosome_indicated_by_filename( filename ), compression_type )) ;
		}
		else if( d.first == "bed" ) {
			return SNPDataSink::UniquePtr(
				PedigreeMappingBedFileSNPDataSink( filename, 0.9 )
			) ;
		}
		else {
			return SNPDataSink::UniquePtr( new GenFileSNPDataSink( filename, get_chromosome_indicated_by_filename( filename ), compression_type )) ;
		}
	}
	
	SNPDataSink::SNPDataSink():
		m_number_of_samples(0u), m_samples_have_been_set( false ), m_number_of_snps_written(0u)
	{}

	SNPDataSink::~SNPDataSink()
	{}
	
	SNPDataSink& SNPDataSink::set_sample_names( std::size_t number_of_samples, SampleNameGetter name_getter ) {
		set_sample_names_impl( number_of_samples, name_getter ) ;
		if( *this ) {
			m_samples_have_been_set = (*this) ;
			m_number_of_samples = number_of_samples ;
		}
		return (*this) ;
	}
	
	SNPDataSink& SNPDataSink::set_metadata( Metadata const& metadata ) {
		set_metadata_impl( metadata ) ;
		return (*this) ;
	}
	
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
		assert( m_samples_have_been_set ) ;
		assert( number_of_samples == m_number_of_samples ) ;
		write_snp_impl( number_of_samples, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability, info ) ;
		if( *this ) {
			++m_number_of_snps_written ;
		}
		return *this ;
	}

	void SNPDataSink::write_snp_impl(
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
		assert(0) ;
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
		m_genotypes.setZero() ;
		{
			vcf::GenotypeSetter< Eigen::MatrixBase< Eigen::MatrixXd > > setter( m_genotypes ) ;
            if( data_reader.supports( ":genotypes:" )) {
                data_reader.get( ":genotypes:", setter ) ;
            }
		}

		write_snp_impl(
			m_genotypes.rows(),
			id_data.get_SNPID(),
			id_data.get_rsid(),
			id_data.get_position().chromosome(),
			id_data.get_position().position(),
			id_data.get_first_allele(),
			id_data.get_second_allele(),
			genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes, 0ul ),
			genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes, 1ul ),
			genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes, 2ul ),
			info
		) ;
	}
	
	SNPDataSink& SNPDataSink::finalise() {
		finalise_impl() ;
		return *this ;
	}
}
