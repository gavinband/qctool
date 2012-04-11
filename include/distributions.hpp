
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef __GTOOL__DISTRIBUTIONS__HPP__
#define __GTOOL__DISTRIBUTIONS__HPP__

/*
 * Note: classes in this file are designed to match boost math toolkit naming scheme.
 * There are two reasons for this choice:
 * 1. At some point, these may be replaced with the ones from boost.
 * 2. Even if not, we may as well follow an already-defined, well-thought-out api such as the one in boost.
 */

class chi_squared_distribution 
{
	public:
		chi_squared_distribution( double degrees_of_freedom ): m_degrees_of_freedom( degrees_of_freedom ) {} ;
		double degrees_of_freedom() const { return m_degrees_of_freedom ; }
	private:
		double m_degrees_of_freedom ;
} ;

template< typename Distribution >
struct complement_type
{
	complement_type( Distribution const& distribution )
		: m_distribution( distribution )
	{}
	
	double degrees_of_freedom() const { return m_distribution.degrees_of_freedom() ; }
	private:
		Distribution const& m_distribution;
} ;

template< typename Distribution >
complement_type< Distribution > complement( Distribution const& distribution ) {
	return complement_type< Distribution >( distribution ) ;
}

template< typename Distribution >
double quantile( Distribution const& distribution, double x ) ;

#endif
