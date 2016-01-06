
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <utility>
#include <string>
#include <boost/bind.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/get_set.hpp"
#include "components/SNPSummaryComponent/ClusterFitComputation.hpp"
#include "metro/DataRange.hpp"
#include "metro/DataSubset.hpp"
#include "metro/likelihood/MultivariateT.hpp"
#include "metro/likelihood/Mixture.hpp"
#include "metro/ValueStabilisesStoppingCondition.hpp"
#include "metro/log_sum_exp.hpp"

// #define DEBUG_CLUSTERFITCOMPUTATION 1

namespace snp_summary_component {
	ClusterFitComputation::ClusterFitComputation(
		double nu,
		double regularisationVariance,
		double regularisationWeight,
		double call_threshhold
	):
		m_call_threshhold( call_threshhold ),
		m_nu( nu ),
		m_scale( "X:Y" ),
		m_xAxisName( "X" ),
		m_yAxisName( "Y" ),
		m_regularisingSigma( regularisationVariance * Eigen::MatrixXd::Identity( 2, 2 ) ),
		m_regularisingWeight( regularisationWeight )
	{}

	namespace {
		metro::DataSubset compute_data_subset( std::size_t N, boost::function< bool( std::size_t ) > predicate ) {
			metro::DataSubset result ;
			bool in_range = false ;
			int start_of_range = 0 ;
			for( int i = 0; i < N; ++i ) {
				bool in = predicate( i ) ;
				if( in ) {
					if( !in_range ) {
						start_of_range = i ;
						in_range = true ;
					}
				} else {
					if( in_range ) {
						result.add( metro::DataRange( start_of_range, i )) ;
						in_range = false ;
					}
				}
			}
			if( in_range ) {
				result.add( metro::DataRange( start_of_range, N )) ;
			}
			return( result ) ;
		}
		
		bool nonmissing_genotypes_and_intensities(
			std::size_t i,
			int g,
			Eigen::MatrixXd const& genotypes,
			Eigen::MatrixXd const& nonmissingness,
			double call_threshhold
		) {
			// return true if the intensities are there
			// and the genotype meets a hard call threshhold.
			return ( nonmissingness.row(i).sum() == 2 )
				&& (
					g < 3
					? ( genotypes( i, g ) > call_threshhold )
					: ( genotypes.row(i).maxCoeff() < call_threshhold )
				)
			;
		}

		bool nonmissing_intensities(
			std::size_t i,
			Eigen::MatrixXd const& nonmissingness
		) {
			return ( nonmissingness.row(i).sum() == 2 ) ;
		}
	}
	
	void ClusterFitComputation::set_scale( std::string const& scale ) {
		if( scale == "X:Y" ) {
			m_xAxisName = "X" ;
			m_yAxisName = "Y" ;
		} else if( scale == "contrast:logR" ) {
			m_xAxisName = "contrast" ;
			m_yAxisName = "logR" ;
		} else {
			throw genfile::BadArgumentError(
				"ClusterFitComputation::set_transform()",
				"transform=\"" + scale + "\"",
				"Unrecognised transform.  Use \"X:Y\" or \"contrast:logR\""
			) ;
		}
		m_scale = scale ;
	}
	
	void ClusterFitComputation::operator()(
		SNPIdentifyingData const& snp,
		Genotypes const& genotypes,
		SampleSexes const&,
		genfile::VariantDataReader& data_reader,
		ResultCallback callback
	) {
		if( !data_reader.supports( "XY" )) {
			return ;
		}
		
		int const N = genotypes.rows() ;
		
		{
			m_intensities.setZero() ;
			m_nonmissingness.setZero() ;
			genfile::vcf::MatrixSetter< IntensityMatrix > intensity_setter( m_intensities, m_nonmissingness ) ;
			data_reader.get( "XY", intensity_setter ) ;
			assert( m_intensities.rows() == N ) ;
			assert( m_intensities.cols() == 2 ) ;
			assert( m_nonmissingness.rows() == N ) ;
			assert( m_nonmissingness.cols() == 2 ) ;
		}
	
#if DEBUG_CLUSTERFITCOMPUTATION
		std::cerr << "At " << snp << ":\n"
			<< "genotypes =\n"
			<< genotypes.block( 0, 0, 10, 3 )
			<< "\n, intensities =\n"
			<< m_intensities.block( 0, 0, 10, 2 )
			<< ".\n" ;
#endif

		if( m_scale == "X:Y" ) {
			// do nothing.
		} else if( m_scale == "contrast:logR" ) {
			Eigen::MatrixXd transform( 2, 2 ) ;
			transform <<
				1, 1,
				-1, 1
			;
			m_intensities *= transform ;
			m_intensities.col(0).array() /= m_intensities.col(1).array() ;
			m_intensities.col(1).array() = m_intensities.col(1).array().log() ;
		} else {
			assert(0) ;
		}
	
		// We compute log-likelihood under a model which is a equally-weighted
		// mixture of the three clusters.
		// For this purpose we record the parameters etc.
		Eigen::MatrixXd parameters = Eigen::MatrixXd::Constant( 5, 3, std::numeric_limits< double >::quiet_NaN() ) ;
		Eigen::VectorXd counts = Eigen::VectorXd::Zero( 3 ) ;
		Eigen::VectorXd genotypeLLs = Eigen::VectorXd::Zero( 3 ) ;
		metro::DataSubset const nonMissingIntensitiesSubset = compute_data_subset(
			N,
			boost::bind( &nonmissing_intensities, _1, m_nonmissingness )
		) ;
		
		// We also compute the log-likelihood P( intensities | clusters ) on the set of samples that have genotypes
		// and compare it with the loglikelihood
		// P( genotypes, intensities | clusters ) = P( intensities | genotypes, clusters ) P( genotypes | clusters )
		// If the latter is much bigger than the former, this indicates the genotypes add a lot of information, which
		// might indicate a badly clustered SNP.
		metro::DataSubset nonMissingGenotypesAndIntensitySubset ;
		double intensity_genotype_ll = 0 ;
		double intensity_genotype_ll_weight = 0 ;

		metro::likelihood::Mixture< double, Eigen::VectorXd, Eigen::MatrixXd > mixture( m_intensities ) ;

		int numberOfClusters = 0 ;

		callback( "clustering-scale", m_scale ) ;
		
#if DEBUG_CLUSTERFITCOMPUTATION
				std::cerr << "snp: " << snp << ".\n" ;
#endif

		for( int g = 0; g < 3; ++g ) {
			std::string const stub = "g=" + ( g == 3 ? std::string( "NA" ) : genfile::string_utils::to_string( g ) ) ;
			metro::DataSubset subset = compute_data_subset(
				N,
				boost::bind(
					&nonmissing_genotypes_and_intensities,
					_1, g, genotypes, m_nonmissingness, m_call_threshhold
				)
			) ;
			
			counts(g) = subset.size() ;
			callback( stub + ":count", genfile::VariantEntry::Integer( subset.size() ) ) ;
			
			typedef metro::likelihood::MultivariateT< double, Eigen::VectorXd, Eigen::MatrixXd > Cluster ;
			Cluster::UniquePtr cluster( new Cluster( m_intensities, m_nu ) ) ;
			metro::ValueStabilisesStoppingCondition stoppingCondition( 0.01, 100 ) ;
			
			if( cluster->estimate_by_em( subset, stoppingCondition, m_regularisingSigma, m_regularisingWeight ) ) {
				callback( stub + ":nu", m_nu ) ;
				callback( stub + ":mu_" + m_xAxisName, cluster->mean()(0) ) ;
				callback( stub + ":mu_" + m_yAxisName, cluster->mean()(1) ) ;
				callback( stub + ":sigma_" + m_xAxisName + m_xAxisName, cluster->sigma()(0,0) ) ;
				callback( stub + ":sigma_" + m_xAxisName + m_yAxisName, cluster->sigma()(1,0) ) ;
				callback( stub + ":sigma_" + m_yAxisName + m_yAxisName, cluster->sigma()(1,1) ) ;
				
				nonMissingGenotypesAndIntensitySubset.add( subset ) ;

				++numberOfClusters ;

				genotypeLLs(g) = cluster->get_value_of_function() ;

#if DEBUG_CLUSTERFITCOMPUTATION
				std::cerr << "cluster " << g << ", "
					<< "subset: " << subset << "\n"
					<< "count = " << subset.size()
					<< ", parameters = " << std::setprecision( 4 ) << cluster->parameters().transpose()
					<< ", ll = " << cluster->get_value_of_function() << ".\n" ;
#endif

				mixture.add_component( stub, 1, cluster ) ;
			} else {
#if DEBUG_CLUSTERFITCOMPUTATION
				std::cerr << "!! No convergence.\n" ;
#endif
				if( subset.size() > 0 ) {
					std::cerr << "!! For SNP " << snp << ", cluster " << g << " has size " << subset.size() << " but distribution did not converge.\n" ;
				}
				genfile::MissingValue const NA = genfile::MissingValue();
				callback( stub + ":nu", NA ) ;
				callback( stub + ":mu_" + m_xAxisName, NA ) ;
				callback( stub + ":mu_" + m_yAxisName, NA ) ;
				callback( stub + ":sigma_" + m_xAxisName + m_xAxisName, NA ) ;
				callback( stub + ":sigma_" + m_xAxisName + m_yAxisName, NA ) ;
				callback( stub + ":sigma_" + m_yAxisName + m_yAxisName, NA ) ;
			}

			callback( stub + ":iterations", genfile::VariantEntry::Integer( stoppingCondition.iterations() ) ) ;
		}

		// Now output the loglikelihoods under a model conditional on genotype...
		callback( "number-of-clusters", numberOfClusters ) ;
		callback( "informative-sample-count", genfile::VariantEntry::Integer( nonMissingGenotypesAndIntensitySubset.size() ) ) ;
		callback( "total-sample-count", genotypes.rows() ) ;
		callback( "ll-given-genotype",  genotypeLLs.sum() ) ;
		// And under equal-weighted mixtures...
		mixture.set_data( m_intensities ) ;
		mixture.evaluate_at( mixture.parameters(), nonMissingGenotypesAndIntensitySubset ) ;
		double const mixtureLL = mixture.get_value_of_function() ;
#if DEBUG_CLUSTERFITCOMPUTATION
		std::cerr << "numberOfClusters = " << numberOfClusters << "\n"
			<< "informative sample count = " << nonMissingGenotypesAndIntensitySubset.size() << "\n" ;
		std::cerr << "Mixture parameters = " << mixture.parameters().transpose() << ".\n" ;
		std::cerr << "mixtureLL = " << mixtureLL << ".\n" ;
		Eigen::VectorXd terms( m_intensities.rows() ) ;
		mixture.get_terms_of_function( terms ) ;
		std::cerr << "terms = " << terms.transpose() << ".\n" ;
#endif
		callback( "equal-weighted-mixture:genotyped_samples:ll", mixtureLL ) ;
		mixture.evaluate_at( mixture.parameters(), nonMissingIntensitiesSubset ) ;
		callback( "equal-weighted-mixture:all_samples:ll",  mixture.get_value_of_function() ) ;
		callback(
			"ll-comparison",
			( genotypeLLs.sum() - ( nonMissingGenotypesAndIntensitySubset.size() * std::log( numberOfClusters )) - mixtureLL )
		) ;
	}

	std::string ClusterFitComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "ClusterFitComputation" ;
	}
}
