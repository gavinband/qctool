
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

//#define DEBUG_CLUSTERFITCOMPUTATION 1

namespace snp_summary_component {
	ClusterFitComputation::ClusterFitComputation(
		double nu,
		double regularisationVariance,
		double regularisationWeight,
		double call_threshhold
	):
		m_call_threshhold( call_threshhold ),
		m_nu( nu ),
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
			genfile::vcf::MatrixSetter< IntensityMatrix > intensity_setter( m_intensities, m_nonmissingness ) ;
			data_reader.get( "XY", intensity_setter ) ;
			assert( m_intensities.rows() == N ) ;
			assert( m_intensities.cols() == 2 ) ;
			assert( m_nonmissingness.rows() == N ) ;
			assert( m_nonmissingness.cols() == 2 ) ;
		}
	
		// We compute log-likelihood under a model which is a equally-weighted
		// mixture of the three clusters.
		// For this purpose we record the parameters etc.
		Eigen::MatrixXd parameters = Eigen::MatrixXd::Constant( 5, 3, std::numeric_limits< double >::quiet_NaN() ) ;
		Eigen::VectorXd linearParameters = Eigen::VectorXd::Constant( 15, std::numeric_limits< double >::quiet_NaN() ) ;
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

		for( int g = 0; g < 3; ++g ) {
			metro::DataSubset subset = compute_data_subset(
				N,
				boost::bind(
					&nonmissing_genotypes_and_intensities,
					_1, g, genotypes, m_nonmissingness, m_call_threshhold
				)
			) ;
			
			counts(g) = subset.size() ;
			
			typedef metro::likelihood::MultivariateT< double, Eigen::VectorXd, Eigen::MatrixXd > Cluster ;
			Cluster::UniquePtr cluster( new Cluster( m_intensities, m_nu ) ) ;
			metro::ValueStabilisesStoppingCondition stoppingCondition( 0.01, 100 ) ;
			
			if( cluster->estimate_by_em( subset, stoppingCondition, m_regularisingSigma, m_regularisingWeight ) ) {
				std::string const stub = "g=" + ( g == 3 ? std::string( "NA" ) : genfile::string_utils::to_string( g ) ) ;
				callback( stub + ":count", genfile::VariantEntry::Integer( subset.size() ) ) ;
				callback( stub + ":nu", m_nu ) ;
				callback( stub + ":mu_X", cluster->get_mean()(0) ) ;
				callback( stub + ":mu_Y", cluster->get_mean()(1) ) ;
				callback( stub + ":sigma_XX", cluster->get_sigma()(0,0) ) ;
				callback( stub + ":sigma_XY", cluster->get_sigma()(1,0) ) ;
				callback( stub + ":sigma_YY", cluster->get_sigma()(1,1) ) ;
				callback( stub + ":iterations", genfile::VariantEntry::Integer( stoppingCondition.iterations() ) ) ;
				
				nonMissingGenotypesAndIntensitySubset.add( subset ) ;

				linearParameters.segment( numberOfClusters * 5, 5 ) = cluster->get_parameters().transpose() ;
				++numberOfClusters ;

				genotypeLLs(g) = cluster->get_value_of_function() ;

#if DEBUG_CLUSTERFITCOMPUTATION
				std::cerr << "cluster " << g << ": count = " << subset.size()
					<< ", parameters = " << std::setprecision( 4 ) << cluster->get_parameters().transpose()
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
			}
			
			// Now output the loglikelihoods under a model conditional on genotype...
			callback( "number-of-clusters", numberOfClusters ) ;
			callback( "cluster-informative-sample-count", genfile::VariantEntry::Integer( nonMissingGenotypesAndIntensitySubset.size() ) ) ;
			callback( "per-cluster:average_ll",  genotypeLLs.sum() ) ;
			// And under equal-weighted mixtures...
			mixture.evaluate_at( linearParameters.segment( 0, numberOfClusters * 5 ), nonMissingGenotypesAndIntensitySubset ) ;
			callback( "equal-weighted-mixture:ll",  mixture.get_value_of_function() ) ;
			mixture.evaluate_at( linearParameters.segment( 0, numberOfClusters * 5 ), nonMissingIntensitiesSubset ) ;
			callback( "equal-weighted-mixture:intensities-ll",  mixture.get_value_of_function() ) ;
		}
	}

	std::string ClusterFitComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "ClusterFitComputation" ;
	}
}
