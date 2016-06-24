
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include "genfile/Chromosome.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils/string_utils.hpp"

namespace genfile {
	std::ostream& operator<<( std::ostream& oStream, Chromosome const& chromosome ) {
		return oStream << static_cast< std::string >( chromosome ) ;
	}

	std::istream& operator>>( std::istream& inStream, Chromosome& chromosome ) {
		std::string s;
		inStream >> s ;
		chromosome = Chromosome( s ) ;
		return inStream ;
	}

	Chromosome::operator std::string () const {
		if( m_repr ) {
			return *m_repr ;
		} else {
			return "NA" ;
		}
	}

	bool Chromosome::is_sex_determining() const {
        return m_repr && ((*m_repr) == "0X" || (*m_repr) == "0Y" || (*m_repr) == "X" || (*m_repr) == "Y" ) ;
    }

	bool Chromosome::is_autosome() const {
		bool result = false ;
		if( m_repr ) {
			try {
				int c = string_utils::to_repr< int > ( *m_repr ) ;
				result = (c > 0 && c <= 22) ;
			}
			catch( string_utils::StringConversionError const& ) {
			}
		}
		return result ;
	}

	bool Chromosome::is_missing() const {
		return (!m_repr) ;
	}
}
