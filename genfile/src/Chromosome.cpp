
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
		if( m_chromosome_e == XYPseudoAutosomalDNA ) {
			return "XY" ;
		}
		else if( m_chromosome_e == MitochondrialDNA ) {
			return "MT" ;
		}
		else if( m_chromosome_e == XChromosome ) {
			return "0X" ;
		}
		else if( m_chromosome_e == YChromosome ) {
			return "0Y" ;
		}
		else if( m_chromosome_e == UnidentifiedChromosome ) {
			return "NA" ;
		}
		else {
			std::ostringstream oStream ;
			oStream << std::right << std::setfill('0') << std::setw(2) << int( m_chromosome_e ) ;
			return oStream.str() ;
		}
	}

	Chromosome::Chromosome( std::string const& s ) {
		if( s == "??" || s == "NA" ) {
			m_chromosome_e = UnidentifiedChromosome ;
		}
		else if( s == "MT" ) {
			m_chromosome_e = MitochondrialDNA ;
		}
		else if( s == "XY" ) {
			m_chromosome_e = XYPseudoAutosomalDNA ;
		}
		else if( s == "0X" || s == "X" ) {
			m_chromosome_e = XChromosome ;
		}
		else if( s == "0Y" || s == "Y" ) {
			m_chromosome_e = YChromosome ;
		}
		else {
			int i = 0 ;
			std::istringstream istr( s ) ;
			istr >> i ;
			istr.peek() ;
			if( !istr.eof() || i < 1 || i > 22 ) {
				m_chromosome_e = UnidentifiedChromosome ;
				// throw BadArgumentError( "Chromosome::Chromosome()", s ) ;
			}
			else {
				m_chromosome_e = Chromosome( i ) ;
			}
		}
	}
	
	bool Chromosome::operator==( std::string const& other ) const {
		return *this == Chromosome( other ) ;
	}
	
	Chromosome& Chromosome::operator++() {
		if( m_chromosome_e == UnidentifiedChromosome ) {
			m_chromosome_e = Chromosome1 ;
		}
		else if( m_chromosome_e == YChromosome ) {
			m_chromosome_e = XYPseudoAutosomalDNA ;
		}
		else {
			m_chromosome_e = ChromosomeEnum( m_chromosome_e + 1 ) ;
		}
			
		return *this ;
	}
	
	bool Chromosome::is_sex_determining() const {
        return m_chromosome_e == XChromosome || m_chromosome_e == YChromosome ;
    }

	bool Chromosome::is_autosome() const {
		return m_chromosome_e <= Chromosome22 ;
	}

	bool Chromosome::is_missing() const {
		return m_chromosome_e == UnidentifiedChromosome ;
	}
}
