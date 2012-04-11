
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "genfile/get_set.hpp"

namespace genfile {
	Ignorer ignore() { return Ignorer() ; }	

	void ValueSetter< std::string >::operator()( char const c ) {
		m_t.assign( std::size_t(1), c ) ;
	}

	GenotypeSetter< std::vector< double > >::GenotypeSetter( std::vector< double >& genotypes ):
		m_genotypes( genotypes )
	{}

	void GenotypeSetter< std::vector< double > >::operator()( std::size_t i, double AA, double AB, double BB ) const {
		if( m_genotypes.size() < (3 * (i+1))) {
			m_genotypes.resize( 3 * (i+1), 0.0 ) ;
		}
		m_genotypes[3*i + 0] = AA ;
		m_genotypes[3*i + 1] = AB ;
		m_genotypes[3*i + 2] = BB ;
	}
}
