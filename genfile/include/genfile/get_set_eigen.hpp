
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_GET_SET_EIGEN_HPP
#define GENFILE_GET_SET_EIGEN_HPP

#include "../config.hpp"
#if HAVE_EIGEN
#include "Eigen/Core"
#include "genfile/get_set.hpp"

namespace genfile {
	template<>
	struct GenotypeSetter< Eigen::MatrixXd >
	{
		GenotypeSetter( Eigen::MatrixXd& matrix ) ;
		void operator()( std::size_t i, double, double, double ) const ;	
		private:
			Eigen::MatrixXd& m_matrix ;
	} ;

	template<>
	struct GenotypeGetter< Eigen::MatrixXd >
	{
		GenotypeGetter( Eigen::MatrixXd& matrix, std::size_t g ) ;
		double operator()( std::size_t i ) const ;
		private:
			Eigen::MatrixXd& m_matrix ;
			std::size_t const m_g ;
	} ;
}

#endif
#endif
