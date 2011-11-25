#include <iostream>
#include <string>
#include <sstream>
#include "genfile/Chromosome.hpp"

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
		else if( m_chromosome_e == UnidentifiedChromosome ) {
			return "Unknown" ;
		}
		else {
			std::ostringstream oStream ;
			oStream << int( m_chromosome_e ) ;
			return oStream.str() ;
		}
	}

	Chromosome::Chromosome( std::string const& s ) {
		std::cerr << "Constructing chromosome from \"" << s << "\".\n" ;
		if( s == "Unknown" || s == "" ) {
			m_chromosome_e = UnidentifiedChromosome ;
		}
		else if( s == "MT" ) {
			m_chromosome_e = MitochondrialDNA ;
		}
		else if( s == "XY" ) {
			m_chromosome_e = XYPseudoAutosomalDNA ;
		}
		else if( s == "0X" ) {
			m_chromosome_e = XChromosome ;
		}
		else if( s == "0Y" ) {
			m_chromosome_e = YChromosome ;
		}
		else {
			int i ;
			std::istringstream istr( s ) ;
			istr >> i ;
			if( i > 0 && i < 23 ) {
				m_chromosome_e = Chromosome( i ) ;
			}
			else {
				throw ChromosomeNotRecognisedError( s ) ;
			}
		}
	}
}
