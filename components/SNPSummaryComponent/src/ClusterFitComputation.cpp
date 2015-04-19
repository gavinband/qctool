
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <utility>
#include <string>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/get_set.hpp"
#include "components/SNPSummaryComponent/ClusterFitComputation.hpp"
#include "metro/DataRange.hpp"
#include "metro/DataSubset.hpp"
#include "metro/likelihood/MultivariateT.hpp"
#include "metro/ValueStabilisesStoppingCondition.hpp"

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
		m_regularisingSigma( regularisationVariance * Eigen::MatrixXd::Identity( 2, 2 ) ),
		m_regularisingWeight( regularisationWeight )
	{}

	void ClusterFitComputation::operator()(
		SNPIdentifyingData const&,
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
	
		Eigen::RowVectorXd mean ;
		Eigen::MatrixXd covariance ;


		// We compute log-likelihood under a model which is a equally-weighted
		// mixture of the three clusters.
		// For this purpose we record the parameters etc.
		Eigen::MatrixXd parameters = Eigen::MatrixXd::Constant( 3, 5, std::numeric_limits< double >::quiet_NaN() ) ;
		Eigen::VectorXd counts = Eigen::VectorXd::Zero( 3 ) ;
		metro::DataSubset nonMissingSubset ;
		
		for( int g = 0; g < 3; ++g ) {
			metro::DataSubset subset ;
				
			// We 
			{
				bool in_range = false ;
				int start_of_range = 0 ;
				for( int i = 0; i < N; ++i ) {
					bool inCluster
						= ( m_nonmissingness.row(i).sum() == 2 )
						&& (
							g < 3
							? ( genotypes( i, g ) < m_call_threshhold )
							: ( genotypes.row(i).maxCoeff() < m_call_threshhold )
						)
					;
					if( inCluster ) {
						if( !in_range ) {
							start_of_range = i ;
							in_range = true ;
						}
					} else {
						if( in_range ) {
							subset.add( metro::DataRange( start_of_range, i )) ;
							in_range = false ;
						}
					}
				}
				if( in_range ) {
					subset.add( metro::DataRange( start_of_range, N )) ;
				}
			}
			
			counts(g) = subset.size() ;
			
			metro::likelihood::MultivariateT< double, Eigen::VectorXd, Eigen::MatrixXd > cluster( m_intensities, subset, m_nu ) ;
			metro::ValueStabilisesStoppingCondition stoppingCondition( 0.01, 100 ) ;
			
			if( cluster.estimate_by_em( stoppingCondition, m_regularisingSigma, m_regularisingWeight ) ) {
				std::string const stub = "g=" + ( g == 3 ? std::string( "NA" ) : genfile::string_utils::to_string( g ) ) ;
				callback( stub + ":count", genfile::VariantEntry::Integer( subset.size() ) ) ;
				callback( stub + ":nu", m_nu ) ;
				callback( stub + ":mu_X", cluster.get_mean()(0) ) ;
				callback( stub + ":mu_Y", cluster.get_mean()(1) ) ;
				callback( stub + ":sigma_XX", cluster.get_sigma()(0,0) ) ;
				callback( stub + ":sigma_XY", cluster.get_sigma()(1,0) ) ;
				callback( stub + ":sigma_YY", cluster.get_sigma()(1,1) ) ;
				callback( stub + ":iterations", genfile::VariantEntry::Integer( stoppingCondition.iterations() ) ) ;
				
#if DEBUG_CLUSTERFITCOMPUTATION
				std::cerr << "cluster " << g << ": count = " << subset.size() << ", parameters = " << cluster.get_parameters().transpose() << ", ll = " << cluster.get_value_of_function() << ".\n" ;
#endif
				nonMissingSubset.add( subset ) ;
				parameters.row(g) = cluster.get_parameters().transpose() ;
			} else {
				std::cerr << "!! No convergence.\n" ;
			}
		}

		// Compute the log-likelihood of all the intensity data under a model
		// which is an equally-weighted mixture of the three clusters.
		{
			double ll = 0 ;
			double ll_weight = 0 ;
			for( int g = 0; g < 3; ++g ) {
				if( counts(g) > 0 ) {
					metro::likelihood::MultivariateT< double, Eigen::VectorXd, Eigen::MatrixXd > cluster( m_intensities, nonMissingSubset, m_nu ) ;
					cluster.evaluate_at( parameters.row(g).transpose() ) ;
					ll += cluster.get_value_of_function() ;
					ll_weight += 1 ;
				}
			}

			callback( "cluster-model:number-of-clusters", ll_weight ) ;
			callback( "cluster-model:loglikelihood", ll ) ;
		}
		
	}

	std::string ClusterFitComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "ClusterFitComputation" ;
	}
}
