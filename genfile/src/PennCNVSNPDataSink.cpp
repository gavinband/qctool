
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/PennCNVSNPDataSink.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/vcf/get_set_eigen.hpp"

namespace genfile {
	
	PennCNVSNPDataSink::PennCNVSNPDataSink( std::string const& filename, double const threshhold ):
		m_filename( filename ),
		m_threshhold( threshhold ),
		m_genotype_alleles( 4 )
	{
		setup( filename ) ;
	}

	PennCNVSNPDataSink::SinkPos PennCNVSNPDataSink::get_stream_pos() const {
		return SinkPos( this, m_stream_ptr->tellp() ) ;
	}

	std::string PennCNVSNPDataSink::get_spec() const { return m_filename ; }

	PennCNVSNPDataSink::operator bool() const { return m_stream_ptr->good() ; }
	
	void PennCNVSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) {
		char const tab = '\t' ;
		stream() << "Name\tChr\tPosition" ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			stream()
				<< tab << getter(i) << ".GType"
				<< tab << getter(i) << ".Log R Ratio"
				<< tab << getter(i) << ".B Allele Freq" ;
		}
		stream() << "\n" ;
		m_genotypes.resize( number_of_samples ) ;
		m_intensities.resize( number_of_samples, 2 ) ;
		m_nonmissingness.resize( number_of_samples, 2 ) ;
		m_intensities.setZero() ;
		m_nonmissingness.setZero() ;
	}
	
	void PennCNVSNPDataSink::write_variant_data_impl(
		VariantIdentifyingData const& variant,
		VariantDataReader& data_reader,
		Info const& info
	) {
		assert( variant.number_of_alleles() == 2 ) ;
		assert( variant.get_allele(0).size() == 1 ) ;
		assert( variant.get_allele(1).size() == 1 ) ;

		{
			data_reader.get( ":genotypes:", genfile::vcf::get_threshholded_calls( m_genotypes, m_threshhold, 0, 1, 2, 3 ) ) ;
			
			m_intensities.setZero() ;
			m_nonmissingness.setZero() ;
			vcf::MatrixSetter< Eigen::MatrixXd > intensity_setter( m_intensities, m_nonmissingness ) ;
			data_reader.get( "LRRB", intensity_setter ) ;
			assert( m_intensities.rows() == int( m_genotypes.size() ) ) ;
			assert( m_intensities.cols() == 2 ) ;
			assert( m_nonmissingness.rows() == int( m_genotypes.size() ) ) ;
			assert( m_nonmissingness.cols() == 2 ) ;
		}

		char const tab = '\t' ;
		stream() << variant.get_primary_id() << tab
			<< variant.get_position().chromosome() << tab
			<< variant.get_position().position() ;

		m_genotype_alleles[0] = "NA" ;
		m_genotype_alleles[1] = variant.get_allele(0) + variant.get_allele(0) ;
		m_genotype_alleles[2] = variant.get_allele(0) + variant.get_allele(1) ;
		m_genotype_alleles[3] = variant.get_allele(1) + variant.get_allele(1) ;

		for( int i = 0; i < m_intensities.rows(); ++i ) {
			// PennCNV / QuantiSNP-style files are genotype, log R ratio, then B Allele Freq.
			stream() << tab << m_genotype_alleles[ m_genotypes[i] ] ;
			if( m_nonmissingness(i,0) != 0 ) {
				stream() << tab << m_intensities(i,0) ;
			} else {
				stream() << tab << "NA" ;
			}
			if( m_nonmissingness(i,1) != 0 ) {
				stream() << tab << m_intensities(i,1) ;
			} else {
				stream() << tab << "NA" ;
			}
		}
		stream() << "\n" ;
	}

	void PennCNVSNPDataSink::setup( std::string const& filename ) {
		m_stream_ptr = open_text_file_for_output( filename, get_compression_type_indicated_by_filename( filename )) ;
		if( !(*m_stream_ptr)) {
			throw ResourceNotOpenedError( m_filename ) ;
		}
		assert( *m_stream_ptr ) ;
	}
}

