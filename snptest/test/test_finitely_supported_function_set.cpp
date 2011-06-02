#include "Eigen/Core"
#include "test_case.hpp"
#include "snptest/FinitelySupportedFunctionSet.hpp"

using namespace snptest ;

typedef Eigen::VectorXd Vector ;
typedef Eigen::MatrixXd Matrix ;

AUTO_TEST_CASE( test_sizes ) {
	std::cerr << "test_sizes..." ;
	for( std::size_t support = 0; support < 10; ++support ) {
		Vector v( support ) ;
		for( std::size_t functions = 0; functions < 10; ++functions ) {
			Matrix m( functions, support ) ;
			FinitelySupportedFunctionSet fsfs( v, m ) ;
			TEST_ASSERT( fsfs.get_size_of_support() == support ) ;
			TEST_ASSERT( fsfs.get_number_of_functions() == functions ) ;
			TEST_ASSERT( fsfs.get_support() == v ) ;
			for( std::size_t row = 0; row < functions; ++row ) {
				TEST_ASSERT( fsfs.get_values( row ) == m.row( row )) ;
			}
		}
	}
	std::cerr << "done.\n" ;
}

AUTO_TEST_CASE( test_support ) {
	std::cerr << "test_support..." ;
	for( std::size_t support = 0; support < 10; ++support ) {
		Vector v( support ) ;
		for( std::size_t i = 0; i < support; ++i ) {
			v[i] = i*2 ;
		}
		for( std::size_t functions = 0; functions < 10; ++functions ) {
			Matrix m( functions, support ) ;
			FinitelySupportedFunctionSet fsfs( v, m ) ;
			TEST_ASSERT( fsfs.get_size_of_support() == support ) ;
			TEST_ASSERT( fsfs.get_number_of_functions() == functions ) ;
			TEST_ASSERT( fsfs.get_support() == v ) ;
			for( std::size_t row = 0; row < functions; ++row ) {
				TEST_ASSERT( fsfs.get_values( row ) == m.row( row )) ;
			}
		}
	}
	std::cerr << "done.\n" ;
}

AUTO_TEST_CASE( test_values ) {
	std::cerr << "test_values..." ;
	for( std::size_t support = 0; support < 10; ++support ) {
		Vector v( support ) ;
		for( std::size_t functions = 0; functions < 10; ++functions ) {
			Matrix m( functions, support ) ;
			for( std::size_t i = 0; i < support; ++i ) {
				for( std::size_t j = 0; j < functions; ++j ) {
					m( j, i ) = 50.0 + (j*i)*2.0 ;
				}
			}
			FinitelySupportedFunctionSet fsfs( v, m ) ;
			TEST_ASSERT( fsfs.get_size_of_support() == support ) ;
			TEST_ASSERT( fsfs.get_number_of_functions() == functions ) ;
			TEST_ASSERT( fsfs.get_support() == v ) ;
			for( std::size_t row = 0; row < functions; ++row ) {
				TEST_ASSERT( fsfs.get_values( row ) == m.row( row )) ;
			}
		}
	}
	std::cerr << "done.\n" ;
}
