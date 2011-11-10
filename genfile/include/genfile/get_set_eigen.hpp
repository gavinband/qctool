#ifndef GENFILE_GET_SET_EIGEN_HPP
#define GENFILE_GET_SET_EIGEN_HPP

#include "../config.hpp"
#if HAVE_EIGEN
#include "Eigen/Core"
#include "genfile/get_set.hpp"

namespace genfile {
	template<>
	struct GenotypeSetter< ::Eigen::MatrixXd >
	{
		GenotypeSetter( ::Eigen::MatrixXd& matrix ) ;
		void operator()( std::size_t i, double, double, double ) const ;	
		private:
			::Eigen::MatrixXd& m_matrix ;
	} ;
}

#endif
#endif
