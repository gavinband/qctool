
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_DISTRIBUTION_HPP
#define QCTOOL_DISTRIBUTION_HPP

#include <string>
#include <memory>
#include <boost/function.hpp>

template< typename Vector, typename Matrix >
struct Distribution {
public:
	typedef std::auto_ptr< Distribution > UniquePtr ;
public:
	virtual ~Distribution() {}
	virtual double get_log_density_at( Vector const& ) const = 0 ;
} ;

template< typename Vector, typename Matrix >
struct ParameterisedDistribution: public Distribution< Vector, Matrix > {
public:
	typedef boost::function< void( std::string const& ) > NameSetter ;
	virtual void get_parameter_names( NameSetter ) const = 0 ;
	virtual void get_parameter( std::string const& name, boost::function< void( Matrix const& ) > ) const = 0 ;
	virtual void set_parameter( std::string const& name, Matrix const& ) = 0 ;
	
} ;


#endif
