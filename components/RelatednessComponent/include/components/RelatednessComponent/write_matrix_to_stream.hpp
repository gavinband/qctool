
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef RELATEDNESS_COMPONENT_WRITE_MATRIX_HPP
#define RELATEDNESS_COMPONENT_WRITE_MATRIX_HPP

#include <iostream>
#include <iomanip>

namespace pca {
	template< typename Stream, typename Matrix >
	void write_matrix_to_stream( Stream& stream, Matrix const& matrix ) {
		stream << std::resetiosflags( std::ios::floatfield ) ;
		stream << std::setprecision( 3 ) ;
		for( int i = 0; i < matrix.rows(); ++i ) {
			for( int j = 0; j < matrix.rows(); ++j ) {
				if( matrix( i,j ) != matrix( i,j )) {
					stream << std::setw( 10 ) << "NA" ;
				} else if( std::abs( matrix( i, j ) ) < 0.00001 ) {
					stream << std::setw( 10 ) << "0" ;
				}
				else {
					stream << std::setw( 10 ) << matrix(i,j) ;
				}
			}
			stream << "\n" ;
		}
		stream << "\n" ;
	}
}

#endif
