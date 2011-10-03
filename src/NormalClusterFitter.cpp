#include <memory>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "genfile/VariantEntry.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/Error.hpp"
#include "genfile/vcf/get_set.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "ClusterFitter.hpp"
#include "NormalClusterFitter.hpp"

namespace {
	std::vector< std::pair< std::string, std::string > > parse_spec( std::string spec ) {
		std::vector< std::pair< std::string, std::string > > result ;
		std::vector< std::string > elts = genfile::string_utils::split_and_strip_discarding_empty_entries( spec, "," ) ;
		for( std::size_t i = 0; i < elts.size(); ++i ) {
			std::vector< std::string > these_elts = genfile::string_utils::split_and_strip( elts[i], "/" ) ;
			if( these_elts.size() != 2 ) {
				throw genfile::BadArgumentError( "(ClusterFitter.cpp) impl::parse_spec()", "spec=\"" + elts[i] + "\"" ) ;
			}
			result.push_back( std::make_pair( these_elts[0], these_elts[1] )) ;
		}
		return result ;
	}
}

NormalClusterFitter::NormalClusterFitter( appcontext::OptionProcessor const& options ):
	m_options( options ),
	m_spec( parse_spec( options.get_value< std::string >( "-fit-clusters" ))),
	m_call_threshhold( options.get_value< double >( "-call-threshhold" ))
{}

void NormalClusterFitter::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps = number_of_snps ;
}

void NormalClusterFitter::processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& data_reader ) {
	Base::processed_snp( id_data, data_reader ) ;

	std::vector< genfile::VariantEntry > genotypes ;
	Eigen::MatrixXd intensities ;

	Eigen::MatrixXd fit( 2, 9 ) ;
	std::vector< std::size_t > non_missing_counts( 3 ) ;

	for( std::size_t spec_i = 0; spec_i < m_spec.size(); ++spec_i ) {
		genfile::vcf::GenotypeSetter< std::vector< genfile::VariantEntry > > genotype_setter( genotypes, m_call_threshhold ) ;
		genfile::vcf::MatrixSetter< Eigen::MatrixXd > intensity_setter( intensities ) ;
		data_reader
			.get( m_spec[spec_i].first, genotype_setter )
			.get( m_spec[spec_i].second, intensity_setter )
		;
		get_cluster_fit( genotypes, intensities, fit, non_missing_counts ) ;
		send_results( "NormalClusterFit:" + m_spec[ spec_i ].first + "/" + m_spec[ spec_i ].second, non_missing_counts, fit ) ;
	}
}

void NormalClusterFitter::get_cluster_fit(
	std::vector< genfile::VariantEntry > const& genotypes,
	Eigen::MatrixXd const& intensities,
	Eigen::MatrixXd& fit,
	std::vector< std::size_t >& non_missing_counts
) {
	assert( genotypes.size() == m_number_of_samples ) ;
	assert( std::size_t( intensities.cols() ) == m_number_of_samples ) ;
	assert( intensities.rows() == 2 ) ;
	assert( fit.rows() == 2 ) ;
	assert( fit.cols() == 9 ) ;
	assert( std::size_t( intensities.cols() ) == m_number_of_samples ) ;
	assert( non_missing_counts.size() == 3 ) ;
	fit.setZero() ;
	std::vector< Eigen::MatrixXd > vars( 3 ) ;
	for( std::size_t i = 0; i < 3; ++i ) { vars[i].resize( 2, 2 ) ; vars[i].setZero() ; }
	for( std::size_t i = 0; i < m_number_of_samples; ++i ) {
		if( genotypes[i].is_int() && intensities( 0, i ) == intensities( 0, i ) && intensities( 1, i ) == intensities( 1, i ) ) {
			std::size_t g = genotypes[i].as< int >() ;
			fit.col( g*3 ) += intensities.col( i ) ;
			vars[g] += intensities.col( i ) * intensities.col( i ).transpose() ;
			fit.block( 0,g*3+1,2,2 ).noalias() += intensities.col( i ) * intensities.col( i ).transpose() ;
			++non_missing_counts[g] ;
		} else {
			// something missing; ignore this.
		}
	}
	for( std::size_t g = 0; g < 3; ++g ) {
		fit.col( g*3 ) /= non_missing_counts[g] ;
		fit.block( 0, g*3+1, 2, 2 ) -= ( fit.col( g*3 ) * fit.col( g*3 ).transpose() ) ;
		fit.block( 0, g*3+1, 2, 2 ) /= ( non_missing_counts[g] - 1.0 ) ;
	}
	
}

void NormalClusterFitter::end_processing_snps() {
	// nothing to do.
}
