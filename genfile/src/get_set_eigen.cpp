
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "../config.hpp"
#if HAVE_EIGEN
#include "genfile/get_set_eigen.hpp"

namespace genfile {
	GenotypeSetter< Eigen::MatrixXd >::GenotypeSetter( Eigen::MatrixXd& matrix ):
		m_matrix( matrix )
	{}
	
	void GenotypeSetter< Eigen::MatrixXd >::operator()( std::size_t i, double AA, double AB, double BB ) const {
		assert( i < static_cast< std::size_t >( m_matrix.rows() ) ) ;
		m_matrix( i, 0 ) = AA ;
		m_matrix( i, 1 ) = AB ;
		m_matrix( i, 2 ) = BB ;
	}
	
	GenotypeGetter< Eigen::MatrixXd >::GenotypeGetter( Eigen::MatrixXd const& matrix, std::size_t g ):
		m_matrix( matrix ),
		m_g( g )
	{
		assert( m_matrix.cols() == 3 ) ;
		assert( g < 3 ) ;
	}
	
	double GenotypeGetter< Eigen::MatrixXd >::operator()( std::size_t i ) const {
		assert( i < std::size_t( m_matrix.rows() ) ) ;
		return m_matrix( i, m_g ) ;
	}
}

#endif

