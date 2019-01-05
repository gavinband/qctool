
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
	namespace impl {
		namespace {
			char const* humanChromosomes [] = {
				"1", "01", "chr1",
				"2", "02", "chr2",
				"3", "03", "chr3",
				"4", "04", "chr4",
				"5", "05", "chr5",
				"6", "06", "chr6",
				"7", "07", "chr7",
				"8", "08", "chr8",
				"9", "09", "chr9",
				"10", "chr10",
				"11", "chr11",
				"12", "chr12",
				"13", "chr13",
				"14", "chr14",
				"15", "chr15",
				"16", "chr16",
				"17", "chr17",
				"18", "chr18",
				"19", "chr19",
				"20", "chr20",
				"21", "chr21",
				"22", "chr22",
				"X", "0X", "chrX", "23",
				"Y", "0Y", "chrY",
				"MT", "chrMT"
			} ;
		}
		MapGenome human(
			humanChromosomes,
			3*9 + 13*2 + 4 + 3 + 2
		) ;
	}

	namespace impl {
		MapGenome::MapGenome( char const** strings, std::size_t N ) {
			for( std::size_t i = 0; i < N; ++i ) {
				m_map[ strings[i] ] = i ;
			}
		}

		MapGenome::~MapGenome() {
		}
		
		bool MapGenome::compare( Chromosome const& left, Chromosome const& right ) const {
			bool const leftIsMissing = !left.m_repr ;
			bool rightIsMissing = !right.m_repr ;
			if( leftIsMissing || rightIsMissing ) {
				return (!leftIsMissing) && rightIsMissing ;
			}
			MapType::const_iterator li = m_map.find( *left.m_repr ) ;
			MapType::const_iterator ri = m_map.find( *right.m_repr ) ;
			return (
				(( li != m_map.end() ) && (( ri == m_map.end() ) || ( li->second < ri->second )))
				||
				(( li == m_map.end() ) && (( ri == m_map.end() ) && ( left.m_repr < right.m_repr )))
			) ;		
		}
	}
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
