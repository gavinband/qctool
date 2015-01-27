
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef EFFECTPARAMETERNAMEPACK_HPP
#define EFFECTPARAMETERNAMEPACK_HPP

#include <string>
#include <vector>

struct EffectParameterNamePack {
public:
	EffectParameterNamePack() ;
	EffectParameterNamePack(
		std::vector< std::string > betas,
		std::vector< std::string > ses,
		std::vector< std::string > cov
	) ;

	EffectParameterNamePack( EffectParameterNamePack const& other ) ;
	EffectParameterNamePack& operator=( EffectParameterNamePack const& other ) ;

	std::size_t size() const ;
	std::string const& parameter_name( std::size_t i ) const ;
	std::string const& se_name( std::size_t i ) const ;
	std::string const& covariance_name( std::size_t i, std::size_t j ) const ;
	
	bool operator==( EffectParameterNamePack const& other ) const ;
	bool operator!=( EffectParameterNamePack const& other ) const ;
	
	std::string get_summary() const ;

private:
	std::vector< std::string > m_betas ;
	std::vector< std::string > m_ses ;
	std::vector< std::string > m_cov ;
} ;

#endif
