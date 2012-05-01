
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/math/distributions/chi_squared.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "snptest/case_control/LogLikelihood.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"
#include "integration/maximisation.hpp"
#include "components/SNPSummaryComponent/AssociationTest.hpp"
#include "components/SNPSummaryComponent/XChromosomeFrequentistCaseControlAssociationTest.hpp"

XChromosomeFrequentistCaseControlAssociationTest::XChromosomeFrequentistCaseControlAssociationTest(
	genfile::CohortIndividualSource const& samples,
	Vector const& phenotypes,
	Matrix const& covariates,
	bool const with_X_inactivation
):
	FrequentistCaseControlAssociationTest( phenotypes, covariates ),
	m_samples( samples ),
	m_sexes( get_sexes( m_samples ) ),
	m_samples_by_sex( get_samples_by_sex( m_sexes )),
	m_with_X_inactivation( with_X_inactivation ),
	m_single_organ_model( false )
{}

std::vector< char > XChromosomeFrequentistCaseControlAssociationTest::get_sexes( genfile::CohortIndividualSource const& samples ) const {
	std::vector< char > result( samples.get_number_of_individuals(), '.' ) ;
	genfile::CohortIndividualSource::ColumnSpec column_spec = samples.get_column_spec() ;
	if( column_spec.check_for_column( "sex" ) ) {
		if( column_spec[ "sex" ].is_discrete() ) {
			for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
				genfile::CohortIndividualSource::Entry const entry = samples.get_entry( i, "sex" ) ;
				if( !entry.is_missing() ) {
					std::string const sex = genfile::string_utils::to_lower( entry.as< std::string >() ) ;
					if( sex == "m" || sex == "male" ) {
						result[i] = 'm' ;
					} else if( sex == "f" || sex == "female" ) {
						result[i] = 'f' ;
					}
					else {
						throw genfile::MalformedInputError( samples.get_source_spec(), i+2, column_spec.find_column( "sex" )) ;
					}
				}
			}
		}
		else {
			std::cerr << "!! (SNPSummaryComputationManager::get_sexes): sex column found but it has the wrong type!\n" ;
			throw genfile::MalformedInputError( samples.get_source_spec(), 1, column_spec.find_column( "sex" )) ;
		}
	}
	return result ;
}

std::map< char, std::vector< int > > XChromosomeFrequentistCaseControlAssociationTest::get_samples_by_sex( std::vector< char > const& sexes ) const {
	std::map< char, std::vector< int > > result ;
	result[ 'm' ] ;
	result[ 'f' ] ;
	result[ '.' ] ;

	for( std::size_t i = 0; i < sexes.size(); ++i ) {
		result[ sexes[i] ].push_back( i ) ;
	}

	return result ;
}

void XChromosomeFrequentistCaseControlAssociationTest::operator()(
	SNPIdentifyingData const& snp,
	Matrix const& genotypes,
	SampleSexes const& sexes, 
	genfile::VariantDataReader&,
	ResultCallback callback
) {
	assert( sexes.size() == genotypes.rows() ) ;
	if( snp.get_position().chromosome() != genfile::Chromosome( "0X" ) ) {
		return ;
	}
	
	Matrix predictor_probs = genotypes ;
	Vector predictor_levels ;

	// 1. We assume males are coded as 0/1 genotypes and samples with missing sex have missing genotypes.
	// (This is done in SNPSummaryComputationManager for this.)

	// 2. If we are doing X chromosome inactivation, replace female call appropriately.
	if( m_with_X_inactivation ) {
		if( m_single_organ_model ) {
			// We assume all the cells affecting the trait have the same active allele (but we don't know which one it is.)
			// Thus, we must integrate out the uncertainty in the active allele.
			// This boils down to recoding females as 0 or 1 with probability
			// c_0 + 1/2 c_1  and  1/2 c_1 + c_2
			// where c_0, c_1, c_2, are the three call probabilities.
			for( int i = 0; i < genotypes.rows(); ++i ) {
				if( sexes[i] == 'f' ) {
					predictor_probs( i, 1 ) *= 0.5 ;
					predictor_probs( i, 0 ) += predictor_probs( i, 1 ) ;
					predictor_probs( i, 1 ) += predictor_probs( i, 2 ) ;
					predictor_probs( i, 2 ) = 0 ;
				} else if( sexes[i] == '.' ) {
					predictor_probs.row( i ).setZero() ;
				}
			}
			predictor_probs.resize( predictor_probs.rows(), 2 ) ;
			predictor_levels = Vector::LinSpaced( 2, 0, 1 ) ;
		} else {
			// We model hemizygote males like homozygote females, and heterozygote females half-way in-between.
			// Replace male genotypes with homozygous ones:
			for( int i = 0; i < genotypes.rows(); ++i ) {
				if( sexes[i] == 'm' ) {
					predictor_probs( i, 2 ) = predictor_probs( i, 1 ) ;
					predictor_probs( i, 1 ) = 0 ;
				} else if( sexes[i] == '.' ) {
					predictor_probs.row( i ).setZero() ;
				}
			}
			predictor_levels = Vector::LinSpaced( 3, 0, 1 ) ;
		}
	} else {
		predictor_levels = Vector::LinSpaced( 3, 0, 2 ) ;
	}

	test( snp, predictor_probs, predictor_levels, callback ) ;
}

void XChromosomeFrequentistCaseControlAssociationTest::set_model( std::string const& model ) {
	if( model == "single-organ" ) {
		m_single_organ_model = true ;
	} else if( model == "normal" ) {
		m_single_organ_model = false ;
	} else {
		throw genfile::BadArgumentError( "XChromosomeFrequentistCaseControlAssociationTest::set_model()", "model=\"" + model + "\"" ) ;
	}
}


