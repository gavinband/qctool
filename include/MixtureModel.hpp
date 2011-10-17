#ifndef QCTOOL_MIXTURE_MODEL_HPP
#define QCTOOL_MIXTURE_MODEL_HPP

#include <vector>
#include <numeric>
#include <boost/ptr_container/ptr_vector.hpp>
#include "Distribution.hpp"

template< typename Vector, typename Matrix >
struct MixtureModel {
public:
	typedef Distribution< Vector, Matrix > ComponentDistribution ;

public:
	std::size_t get_number_of_components() const { return m_components.size() ; }
	
	void add_component( typename ComponentDistribution::UniquePtr distribution, double weight ) {
		assert( weight == weight ) ;
		m_components.push_back( distribution ) ;
		m_weights.push_back( weight ) ;
		double total_weight = std::accumulate( m_weights.begin(), m_weights.end(), 0.0 ) ;
		for( std::size_t i = 0; i < m_weights.size(); ++i ) {
			m_weights[i] /= total_weight ;
		}
	}

	double get_log_density( Vector const& v ) const {
		double result = 0.0 ;
		for( std::size_t i = 0; i < m_components.size(); ++i ) {
			result += std::log( m_weights[i] ) + m_components[i].get_log_density_at( v ) ;
		}
		return result ;
	}

	double get_log_density_given_component( std::size_t component, Vector const& v ) const {
		return m_components[ component ].get_log_density_at( v ) ;
	}

private:
	boost::ptr_vector< ComponentDistribution > m_components ;
	std::vector< double > m_weights ;
} ;

#endif
