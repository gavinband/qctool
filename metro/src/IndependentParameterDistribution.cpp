
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <exception>
#include "Eigen/Core"
#include "metro/distributions/Flat.hpp"
#include "metro/IndependentParameterDistribution.hpp"

namespace metro {
	IndependentParameterDistribution::UniquePtr IndependentParameterDistribution::create(
		std::vector< std::string > const& parameter_names
	) {
		return IndependentParameterDistribution::UniquePtr(
			new IndependentParameterDistribution( parameter_names )
		) ;
	}
	
	IndependentParameterDistribution::IndependentParameterDistribution(
		std::vector< std::string > const& parameter_names
	):
		m_parameter_names( parameter_names )
	{
		for( std::size_t i = 0; i < m_parameter_names.size(); ++i ) {
			m_parameter_name_indices[ m_parameter_names[i] ] = i ;
			m_densities.push_back( new metro::distributions::Flat() ) ;
		}
	}
	
	void IndependentParameterDistribution::set_prior( std::string const& name, UnivariateLogDensity::UniquePtr pdf ) {
		ParameterNameIndices::const_iterator where = m_parameter_name_indices.find( name ) ;
		if( where == m_parameter_name_indices.end() ) {
			throw std::invalid_argument( "Unexpected parameter \"" + name + "\"." ) ;
		}
		std::size_t const i = m_parameter_name_indices[name] ;
		m_densities.replace( i, pdf.release() ) ;
	}
	
	int IndependentParameterDistribution::number_of_parameters() const {
		return m_parameter_names.size() ;
	};
	
	void IndependentParameterDistribution::evaluate_at(
		Point const& parameters,
		int const numberOfDerivatives
	) {
		m_parameters = parameters ;
		for( std::size_t i = 0; i < m_densities.size(); ++i ) {
			m_densities[i].evaluate_at( parameters(i) ) ;
		}
	}
	
	double IndependentParameterDistribution::get_value_of_function() const {
		double result = 0.0 ;
		for( std::size_t i = 0; i < m_densities.size(); ++i ) {
			result += m_densities[i].get_value_of_function() ;
		}
		return result ;
	}
	
	IndependentParameterDistribution::Vector IndependentParameterDistribution::get_value_of_first_derivative() const {
		Vector result = Vector::Zero( m_parameters.size() ) ;
		for( std::size_t i = 0; i < m_densities.size(); ++i ) {
			result(i) = m_densities[i].get_value_of_first_derivative() ;
		}
		return result ;
	}
	
	IndependentParameterDistribution::Matrix IndependentParameterDistribution::get_value_of_second_derivative() const {
		Matrix result = Matrix::Zero( m_parameters.size(), m_parameters.size() ) ;
		for( std::size_t i = 0; i < m_densities.size(); ++i ) {
			result(i,i) = m_densities[i].get_value_of_second_derivative() ;
		}
		return result ;
	}
	
	std::string IndependentParameterDistribution::get_summary() const {
		std::string result = "per-parameter priors:\n" ;
		for( std::size_t i = 0; i < m_parameter_names.size(); ++i ) {
			std::string const& name = m_parameter_names[i] ;
			ParameterNameIndices::const_iterator where = m_parameter_name_indices.find( name ) ;
			result += ((i>0) ? "\n        " : "        " ) + name + ": " + m_densities[where->second].get_summary() ;
		}
		return result ;
	}
}
