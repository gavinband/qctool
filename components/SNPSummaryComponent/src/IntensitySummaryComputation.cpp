
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
#include "components/SNPSummaryComponent/IntensitySummaryComputation.hpp"
#include "metro/mean_and_covariance.hpp"

IntensitySummaryComputation::IntensitySummaryComputation( double call_threshhold ):
	m_call_threshhold( call_threshhold )
{}

namespace {
	template< typename Data >
	std::pair< double, double > compute_mean_and_variance_compensated( Data const& data ) {
		double const mean = data.sum() / data.size() ;
		double compensation = ( data.array() - mean ).sum() ;
		compensation = ( compensation * compensation ) / data.size() ;
		double const variance = ( ( ( data.array() - mean ).square() ).sum() * compensation ) / ( data.size() - 1 ) ;
		return std::make_pair( mean, variance ) ;
	}
}

void IntensitySummaryComputation::operator()(
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
		genfile::vcf::MatrixSetter< IntensityMatrix > intensity_setter( m_intensities ) ;
		data_reader.get( "XY", intensity_setter ) ;
		assert( m_intensities.rows() == N ) ;
		assert( m_intensities.cols() == 2 ) ;
		m_nonmissingness.resize( m_intensities.rows(), m_intensities.cols() ) ;
		m_nonmissingness = ( m_intensities.array() == m_intensities.array() ).cast< double >() ;
	}
	
	Eigen::RowVectorXd mean ;
	Eigen::MatrixXd covariance ;

	metro::compute_mean_and_covariance(
		m_intensities,
		m_nonmissingness,
		mean,
		covariance
	) ;
	
	callback( "mean_X", mean(0) ) ;
	callback( "mean_Y", mean(1) ) ;
	callback( "var_X", covariance(0,0) ) ;
	callback( "var_Y", covariance(1,1) ) ;
	callback( "cov_XY",covariance(0,1) ) ;
	
	for( int g = 0; g < 3; ++g ) {
		m_intensities_by_genotype = m_intensities.array() ;
		for( int i = 0; i < N; ++i ) {
			
		}
		m_intensities_by_genotype.setConstant( std::numeric_limits< double >::quiet_NaN() ) ;
		m_nonmissingness_by_genotype.resize( N, 6 ) ;
		m_nonmissingness_by_genotype.setConstant( 0 ) ;
		
		metro::compute_mean_and_covariance(
			m_intensities_by_genotype.block( 2*g, 0, 2, N ).transpose(),
			m_nonmissingness.block( 2*g, 0, 2, N  ).transpose(),
			mean,
			covariance
		) ;
		std::string const stub = "g=" + genfile::string_utils::to_string( g ) ;
		callback( stub + ":mean_X", mean(0) ) ;
		callback( stub + ":mean_Y", mean(1) ) ;
		callback( stub + ":var_X", covariance(0,0) ) ;
		callback( stub + ":var_Y", covariance(1,1) ) ;
		callback( stub + ":cov_XY",covariance(0,1) ) ;
	}
}

std::string IntensitySummaryComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
	return prefix + "IntensitySummaryComputation" ;
}
