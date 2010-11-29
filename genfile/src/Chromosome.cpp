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
			return "Unknown" ;
		}
		else {
			std::ostringstream oStream ;
			oStream << std::right << std::setfill('0') << std::setw(2) << int( m_chromosome_e ) ;
			return oStream.str() ;
		}
	}

	Chromosome::Chromosome( std::string const& s ) {
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
				throw BadArgumentError( Chromosome::Chromosome(), s ) ;
			}
		}
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
}
