
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_INTENSITY_MODEL_HPP
#define QCTOOL_INTENSITY_MODEL_HPP

#include <memory>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include "Distribution.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"

// An IntensityModel represents a model for a set of X/Y intensities at a variant.
// It provides a way to obtain the 
struct IntensityModel {
public:
	typedef std::auto_ptr< IntensityModel > UniquePtr ;
	typedef boost::shared_ptr< IntensityModel > SharedPtr ;
	static UniquePtr estimate( Eigen::MatrixXd const& intensities, genfile::SingleSNPGenotypeProbabilities const& calls, double threshhold = 0.9 ) ;

public:
	typedef Eigen::Vector2d Vector ;
	typedef Eigen::Matrix2d Matrix ;
	virtual ~IntensityModel() {}
	virtual std::size_t get_number_of_genotype_classes() const = 0 ;
	virtual double get_log_density_at_intensity( Vector const& ) const = 0 ;
	virtual double get_log_density_at_intensity_given_genotype( std::size_t genotype, Vector const& ) const = 0 ;
} ;

#endif