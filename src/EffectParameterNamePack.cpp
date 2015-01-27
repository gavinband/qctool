
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <cassert>
#include "genfile/string_utils/string_utils.hpp"
#include "EffectParameterNamePack.hpp"

EffectParameterNamePack::EffectParameterNamePack() {}
EffectParameterNamePack::EffectParameterNamePack(
	std::vector< std::string > betas,
	std::vector< std::string > ses,
	std::vector< std::string > cov
):
	m_betas( betas ),
	m_ses( ses ),
	m_cov( cov )
{
	assert( m_ses.size() == m_betas.size() ) ;
	assert( m_cov.size() == ( m_betas.size() * ( m_betas.size() - 1 )) / 2 ) ;
}

EffectParameterNamePack::EffectParameterNamePack( EffectParameterNamePack const& other ):
	m_betas( other.m_betas ),
	m_ses( other.m_ses ),
	m_cov( other.m_cov )
{
}

EffectParameterNamePack& EffectParameterNamePack::operator=( EffectParameterNamePack const& other ) {
	m_betas = other.m_betas ;
	m_ses = other.m_ses ;
	m_cov = other.m_cov ;
	return *this ;
}

std::size_t EffectParameterNamePack::size() const {
	return m_betas.size() ;
}

std::string const& EffectParameterNamePack::parameter_name( std::size_t i ) const {
	return m_betas[i] ;
}

std::string const& EffectParameterNamePack::se_name( std::size_t i ) const {
	return m_ses[i] ;
}

std::string const& EffectParameterNamePack::covariance_name( std::size_t i, std::size_t j ) const {
	assert( j > i ) ;
	// index in cov skips the lower diagonal.
	// Lower diagonal has ((i+1)*(i+2)/2) entries since i is a 0-based index.
	// index in full array would be i*N+j
	std::size_t const N = m_betas.size() ;
	std::size_t index = i*N + j - ( (i+1)*(i+2)/2 ) ;
	return m_cov[index] ;
}

bool EffectParameterNamePack::operator==( EffectParameterNamePack const& other ) const {
	return m_betas == other.m_betas && m_ses == other.m_ses && m_cov == other.m_cov ;
}

bool EffectParameterNamePack::operator!=( EffectParameterNamePack const& other ) const {
	return !operator==( other ) ;
}

std::string EffectParameterNamePack::get_summary() const {
	using genfile::string_utils::join ;
	std::ostringstream ostr ;
	ostr << "(" << join( m_betas, "," ) << "; " << join( m_ses, "," ) << "; " << join( m_cov, ", " ) << ")" ;
	return ostr.str() ;
}
