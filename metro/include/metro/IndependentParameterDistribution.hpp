
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNPTEST_REGRESSION_INDEPENDENT_INDEPENDENTPARAMETERDISTRIBUTION_HPP
#define SNPTEST_REGRESSION_INDEPENDENT_INDEPENDENTPARAMETERDISTRIBUTION_HPP

#include <memory>
#include <unordered_map>
#include <Eigen/Core>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include "metro/UnivariateLogDensity.hpp"
#include "metro/SmoothFunction.hpp"

namespace metro {
	struct IndependentParameterDistribution: public SmoothFunction {
	public:
		typedef std::auto_ptr< IndependentParameterDistribution > UniquePtr ;
		typedef Eigen::VectorXd Point ;
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
		typedef Eigen::Matrix< int, Eigen::Dynamic, Eigen::Dynamic > IntegerMatrix ;

		static UniquePtr create( std::vector< std::string > const& parameter_names ) ;

	public:
		IndependentParameterDistribution( std::vector< std::string > const& parameter_names ) ;

		void set_prior( std::string const& name, UnivariateLogDensity::UniquePtr pdf ) ;
		int number_of_parameters() const ;
		void evaluate_at(
			Point const& parameters,
			int const numberOfDerivatives = 2
		) ;
		double get_value_of_function() const ;
		Vector get_value_of_first_derivative() const ;
		Matrix get_value_of_second_derivative() const ;
		
		std::string get_summary() const ;
		
	private:
		std::size_t m_number_of_parameters ;
		std::vector< std::string > const m_parameter_names ;
		typedef std::map< std::string, std::size_t > ParameterNameIndices ;
		ParameterNameIndices m_parameter_name_indices ;
		Vector m_parameters ;
		boost::ptr_vector< UnivariateLogDensity > m_densities ;
	} ;
}

#endif

