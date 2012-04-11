
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_NORMAL_CLUSTER_FIT_COMPARER_HPP
#define QCTOOL_NORMAL_CLUSTER_FIT_COMPARER_HPP

#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "ClusterFitter.hpp"

struct NormalClusterFitComparer: public boost::noncopyable
{
	typedef std::auto_ptr< NormalClusterFitComparer > UniquePtr ;
	virtual ~NormalClusterFitComparer() {}
	virtual Eigen::MatrixXd compare( Eigen::MatrixXd const& intensities, Eigen::MatrixXd const& fit1, Eigen::MatrixXd const& fit2 ) const = 0 ;
} ;


#endif
