
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
#include "components/SNPSummaryComponent/IntensityReporter.hpp"
#include "metro/mean_and_covariance.hpp"

namespace snp_summary_component {
	IntensityReporter::IntensityReporter(
		genfile::CohortIndividualSource const& samples,
		double call_threshhold
	):
		m_samples( samples ),
		m_call_threshhold( call_threshhold )
	{
		m_samples.get_column_values(
			std::string( "ID_1" ),
			boost::bind(
				&std::vector< genfile::VariantEntry >::push_back,
				&m_sample_ids,
				_2
			)
		) ;
	}

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

	void IntensityReporter::operator()(
		SNPIdentifyingData const&,
		Genotypes const& genotypes,
		SampleSexes const&,
		genfile::VariantDataReader& data_reader,
		ResultCallback callback
	) {
		if( !data_reader.supports( "XY" )) {
			return ;
		}
		
		int const N = m_samples.get_number_of_individuals() ;
		
		{
			genfile::vcf::MatrixSetter< IntensityMatrix > intensity_setter( m_intensities, m_nonmissingness ) ;
			data_reader.get( "XY", intensity_setter ) ;
			assert( m_intensities.rows() == N ) ;
			assert( m_intensities.cols() == 2 ) ;
			assert( m_nonmissingness.rows() == N ) ;
			assert( m_nonmissingness.cols() == 2 ) ;
		}
		
		for( int i = 0; i < N; ++i ) {
			for( int h = 0; h < 2; ++h ) {
				std::string const variable_stub =  m_sample_ids[i].as< std::string >() + ":" ;
				Genotypes::Index g ;
				if( genotypes.row( i ).maxCoeff( &g ) >= m_call_threshhold ) {
					callback( variable_stub + "threshholded_genotype", int( g ) ) ;
				} else {
					callback( variable_stub + "threshholded_genotype", genfile::MissingValue() ) ;
				}
				if( m_nonmissingness( i, h ) ) {
					callback( variable_stub + (( h == 0 ) ? "x" : "y" ) + "_intensity", m_intensities( i, h ) ) ;
				} else {
					callback( variable_stub + (( h == 0 ) ? "x" : "y" ) + "_intensity", genfile::MissingValue() ) ;
				}
			}
		}
	}

	std::string IntensityReporter::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "IntensityReporter" ;
	}
}
