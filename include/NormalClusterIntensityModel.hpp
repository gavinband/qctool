
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_NORMAL_CLUSTER_INTENSITY_MODEL_HPP
#define QCTOOL_NORMAL_CLUSTER_INTENSITY_MODEL_HPP

#include <memory>
#include <vector>
#include <Eigen/Core>
#include "IntensityModel.hpp"
#include "MixtureModel.hpp"

struct NormalClusterIntensityModel: public IntensityModel {
public:
	typedef std::auto_ptr< NormalClusterIntensityModel > UniquePtr ;
public:
	typedef Eigen::Vector2d Vector ;
	typedef Eigen::Matrix2d Matrix ;

	NormalClusterIntensityModel() ;

	void add_cluster( Vector const& mean, Matrix const& variance, double weight ) ;
	
	std::size_t get_number_of_genotype_classes() const ;
	double get_log_density_at_intensity( Vector const& ) const ;
	double get_log_density_at_intensity_given_genotype( std::size_t genotype, Vector const& v ) const ;

private:
	typedef MixtureModel< Vector, Matrix > Model ;
	Model m_model ;
} ;
#endif
