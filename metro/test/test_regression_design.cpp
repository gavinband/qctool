#include <vector>
#include <iostream>
#include <boost/format.hpp>
#include "test_case.hpp"
#include "metro/SampleRange.hpp"
#include "metro/union_ranges.hpp"
#include "metro/RegressionDesign.hpp"

#define DEBUG 1

BOOST_AUTO_TEST_SUITE( test_regressiondesign )
AUTO_TEST_CASE( test_construction ) {
	using metro::SampleRange ;
	using metro::RegressionDesign ;
	typedef RegressionDesign::Matrix Matrix ;
	typedef RegressionDesign::Vector Vector ;

	for( int nSamples = 0; nSamples < 100; ++nSamples ) {
		Vector const outcome_nonmissingness = Vector::Constant( nSamples, 1.0 ) ;
		Vector outcome = Vector::Zero( nSamples ) ;
		for( int i = 1; i < nSamples; i += 2 ) {
			outcome(i) = 1 ;
		}
		
		std::vector< metro::SampleRange > included_samples ;
		included_samples.push_back( metro::SampleRange( 0, nSamples )) ;
	
		for( int nCovariates = 0; nCovariates < 10; ++nCovariates ) {
			Matrix const covariates = Matrix::Constant( nSamples, nCovariates, 1001.0 ) ;
			Matrix const covariate_nonmissingness = Matrix::Constant( nSamples, nCovariates, 1.0 ) ;
			for( int nPredictors = 0; nPredictors < 10; ++nPredictors ) {
				Matrix predictor_levels = Matrix::Constant( 1, nPredictors, 5.0 ) ;
				Matrix predictor_probabilities = Matrix::Constant( nSamples, 1, 1.0 ) ;
			
				RegressionDesign design(
					outcome, outcome_nonmissingness, "outcome",
					covariates, covariate_nonmissingness, std::vector< std::string >( nCovariates, "cov" ),
					std::vector< std::string >( nPredictors, "predictor" )
				) ;
			
				design.set_predictor_levels( predictor_levels, predictor_probabilities, included_samples ) ;

				// check the outcome
				BOOST_CHECK_EQUAL( design.outcome().size(), nSamples ) ;
				for( int i = 0; i < nSamples; ++i ) {
					BOOST_CHECK_EQUAL( design.outcome()[i], ( i % 2 ) ) ;
				}

				Matrix const& matrix = design.set_predictor_level( 0 ).matrix() ;
				BOOST_CHECK_EQUAL( matrix.rows(), nSamples ) ;
				BOOST_CHECK_EQUAL( matrix.cols(), nCovariates + nPredictors + 1 ) ;

				if( nSamples > 0 ) {
					// check the baseline column
					BOOST_CHECK_EQUAL( matrix.col(0).maxCoeff(), 1.0 ) ;
					BOOST_CHECK_EQUAL( matrix.col(0).minCoeff(), 1.0 ) ;

					// check the predictor columns
					for( int i = 0; i < nPredictors; ++i ) {
						BOOST_CHECK_EQUAL( matrix.col(i+1).maxCoeff(), 5.0 ) ;
						BOOST_CHECK_EQUAL( matrix.col(i+1).minCoeff(), 5.0 ) ;
					}

					// check the covariate columns
					for( int i = 0; i < nCovariates; ++i ) {
						BOOST_CHECK_EQUAL( matrix.col(i+1+nPredictors).maxCoeff(), 1001.0 ) ;
						BOOST_CHECK_EQUAL( matrix.col(i+1+nPredictors).minCoeff(), 1001.0 ) ;
					}
				}
			}
		}
	}
}

AUTO_TEST_CASE( test_predictor_levels ) {
	using metro::SampleRange ;
	using metro::RegressionDesign ;
	typedef RegressionDesign::Matrix Matrix ;
	typedef RegressionDesign::Vector Vector ;

	int nSamples = 10;
	int nCovariates = 0 ;
	
	Vector const outcome_nonmissingness = Vector::Constant( nSamples, 1.0 ) ;
	Vector outcome = Vector::Zero( nSamples ) ;
	for( int i = 0; i < nSamples; i += 2 ) {
		outcome(i) = 1 ;
	}
	
	std::vector< metro::SampleRange > included_samples ;
	included_samples.push_back( metro::SampleRange( 0, nSamples )) ;
	
	Matrix const covariates = Matrix::Constant( nSamples, nCovariates, 1001.0 ) ;
	Matrix const covariate_nonmissingness = Matrix::Constant( nSamples, nCovariates, 1.0 ) ;

	for( int nPredictors = 0; nPredictors < 10; ++nPredictors ) {
		RegressionDesign design(
			outcome, outcome_nonmissingness, "outcome",
			covariates, covariate_nonmissingness, std::vector< std::string >( nCovariates, "cov" ),
			std::vector< std::string >( nPredictors, "predictor" )
		) ;

		for( int nPredictorLevels = 0; nPredictorLevels < 10; ++nPredictorLevels ) {
			Matrix predictor_levels = Matrix::Zero( nPredictorLevels, nPredictors ) ;
			for( int i = 0; i < nPredictorLevels; ++i ) {
				for( int j = 0; j < nPredictors; ++j ) {
					predictor_levels( i, j ) = double( (i+1)*(j+1) ) ;
				}
			}

			Matrix predictor_probabilities = Matrix::Constant( nSamples, nPredictorLevels, 1 / double( nPredictorLevels ) ) ;
			design.set_predictor_levels( predictor_levels, predictor_probabilities, included_samples ) ;
			
			for( int iLevel = 0; iLevel < nPredictorLevels; ++iLevel ) {
				Matrix const& matrix = design.set_predictor_level( iLevel ).matrix() ;
			
				if( nSamples > 0 ) {
					for( int j = 0; j < nPredictors; ++j ) {
							BOOST_CHECK_EQUAL( matrix.col(j+1).maxCoeff(), double((iLevel+1)*(j+1)) ) ;
							BOOST_CHECK_EQUAL( matrix.col(j+1).minCoeff(), double((iLevel+1)*(j+1)) ) ;
					}
				}
			}
		}
	}
}

AUTO_TEST_CASE( test_names ) {
	using metro::SampleRange ;
	using metro::RegressionDesign ;

	for( int numberOfPredictors = 0; numberOfPredictors < 10; ++numberOfPredictors ) {
		std::vector< std::string > predictor_names( numberOfPredictors ) ;
		for( int i = 0; i < numberOfPredictors; ++i ) {
			predictor_names[i] = ( boost::format( "predictor%d" ) % (i+1) ).str() ;
		}
		
		for( int numberOfCovariates = 0; numberOfCovariates < 10; ++numberOfCovariates ) {
			Eigen::VectorXd outcome = Eigen::VectorXd::Zero( 5 ), outcome_nm = Eigen::VectorXd::Constant(5, 1) ;
			Eigen::MatrixXd cov = Eigen::MatrixXd::Constant( 5, numberOfCovariates, 2 ), cov_nm = Eigen::MatrixXd::Constant(5, numberOfCovariates, 1 ) ;

			std::vector< std::string > cov_names( numberOfCovariates ) ;
			for( int i = 0; i < numberOfCovariates; ++i ) {
				cov_names[i] = ( boost::format( "cov%d" ) % (i+1)).str() ;
				cov.col(i).setConstant( i ) ;
			}
		
			RegressionDesign design( outcome, outcome_nm, "outcome", cov, cov_nm, cov_names, predictor_names ) ;
			std::vector< std::string > names = design.design_matrix_column_names() ;
			TEST_ASSERT( names.size() == 1 + numberOfPredictors + numberOfCovariates ) ;
			TEST_ASSERT( names[0] == "baseline" ) ;
			for( int i = 0; i < numberOfPredictors; ++i ) {
				TEST_ASSERT( names[1+i] == predictor_names[i] ) ;
			}
			for( int i = 0; i < numberOfCovariates; ++i ) {
				TEST_ASSERT( names[1+numberOfPredictors+i] == cov_names[i] ) ;
			}
			
#if DEBUG
			std::cerr << design.get_summary() << "\n" ;
#endif
			
		}
	}
}

AUTO_TEST_CASE( test_names_interaction ) {
	using metro::SampleRange ;
	using metro::RegressionDesign ;

	for( int numberOfPredictors = 0; numberOfPredictors < 10; ++numberOfPredictors ) {
		std::vector< std::string > predictor_names( numberOfPredictors ) ;
		for( int i = 0; i < numberOfPredictors; ++i ) {
			predictor_names[i] = ( boost::format( "predictor%d" ) % (i+1) ).str() ;
		}
		
		for( int numberOfCovariates = 0; numberOfCovariates < 10; ++numberOfCovariates ) {
			Eigen::VectorXd outcome = Eigen::VectorXd::Zero( 5 ), outcome_nm = Eigen::VectorXd::Constant(5, 1) ;
			Eigen::MatrixXd cov = Eigen::MatrixXd::Constant( 5, numberOfCovariates, 2 ), cov_nm = Eigen::MatrixXd::Constant(5, numberOfCovariates, 1 ) ;

			std::vector< std::string > cov_names( numberOfCovariates ) ;
			for( int i = 0; i < numberOfCovariates; ++i ) {
				cov_names[i] = ( boost::format( "cov%d" ) % (i+1)).str() ;
				cov.col(i).setConstant( i ) ;
			}

			std::vector< int > interactions ;
			for( int interaction = 0; interaction < numberOfCovariates; interaction +=2 ) {
				interactions.push_back( interaction ) ;
			}
		
			RegressionDesign design( outcome, outcome_nm, "outcome", cov, cov_nm, cov_names, predictor_names, RegressionDesign::eIdentity, interactions ) ;

#if DEBUG
			std::cerr << design.get_summary() << "\n" ;
#endif

			std::vector< std::string > names = design.design_matrix_column_names() ;
			BOOST_CHECK_EQUAL( names.size(), 1 + numberOfPredictors * ( 1 + interactions.size() ) + numberOfCovariates ) ;
			BOOST_CHECK_EQUAL( names[0], "baseline" ) ;
			for( int i = 0; i < numberOfPredictors; ++i ) {
				BOOST_CHECK_EQUAL( names[1+i], predictor_names[i] ) ;
			}
			for( int j = 0; j < interactions.size(); ++j ) {
				for( int i = 0; i < numberOfPredictors; ++i ) {
					BOOST_CHECK_EQUAL( names[1+((j+1)*numberOfPredictors)+i], predictor_names[i] + "x" + cov_names[interactions[j]] ) ;
				}
			}
			for( int i = 0; i < numberOfCovariates; ++i ) {
				BOOST_CHECK_EQUAL( names[1+(numberOfPredictors*(1+interactions.size()))+i], cov_names[i] ) ;
			}
			
			
		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
