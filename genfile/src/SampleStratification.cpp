
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include "genfile/Error.hpp"
#include "genfile/SampleRange.hpp"
#include "genfile/SampleStratification.hpp"

namespace genfile {
	SampleStratification::SampleStratification() {}
	SampleStratification::SampleStratification( SampleStratification const& other ):
		m_strata_names( other.m_strata_names ),
		m_index_by_strata( other.m_index_by_strata ),
		m_sample_ranges( other.m_sample_ranges )
	{}
	
	SampleStratification& SampleStratification::operator=( SampleStratification const& other ) {
		m_strata_names = other.m_strata_names ;
		m_index_by_strata = other.m_index_by_strata ;
		m_sample_ranges = other.m_sample_ranges ;
		return *this ;
	}

	void SampleStratification::add_stratum( std::string const& name ) {
		if( m_index_by_strata.find( name ) != m_index_by_strata.end() ) {
			throw genfile::BadArgumentError( "core::SampleStratification::add_stratum()", "name=\"" + name + "\" (already in index)" ) ;
		}
		m_strata_names.push_back( name ) ;
		m_sample_ranges.push_back( std::vector< SampleRange >() ) ;
		m_index_by_strata[ name ] = m_strata_names.size() - 1 ;
	}

	void SampleStratification::add_sample( std::string const& strata_name, int sample ) {
		return add_sample_range( strata_name, SampleRange( sample, sample + 1 ) ) ;
	}

	namespace {
		namespace impl {
			std::vector< SampleRange > union_ranges(
					std::vector< SampleRange > const& left,
					std::vector< SampleRange > const& right
			) {
				std::vector< SampleRange > tmp ;
				tmp.reserve( left.size() + right.size() ) ;
				tmp.insert( tmp.end(), left.begin(), left.end() ) ;
				tmp.insert( tmp.end(), right.begin(), right.end() ) ;
				std::sort( tmp.begin(), tmp.end() ) ;

				std::vector< SampleRange > result ;
				std::vector< SampleRange >::iterator i = tmp.begin() ;
				while( i != tmp.end() ) {
						std::vector< SampleRange >::iterator rightmost_range = i ;
						std::vector< SampleRange >::iterator j = i ;
						// walk right through the regions that overlap or abut this one.
						for( ++j; j != tmp.end() && j->begin() <= rightmost_range->end(); rightmost_range = j++ ) {
								// nothing to do.
						}
						result.push_back( SampleRange( i->begin(), rightmost_range->end() ) ) ;
						i = j ;
				}
				return result ;
			}

			std::vector< SampleRange > union_ranges(
					std::vector< SampleRange > const& left,
					SampleRange const& right
			) {
					return union_ranges( left, std::vector< SampleRange >( 1, right ) ) ;
			}
		}
	}

	void SampleStratification::add_sample_range( std::string const& strata_name, SampleRange const& range ) {
		std::map< std::string, std::size_t >::const_iterator where = m_index_by_strata.find( strata_name ) ;
		if( where == m_index_by_strata.end() ) {
			add_stratum( strata_name ) ;
			where = m_index_by_strata.find( strata_name ) ;
			assert( where != m_index_by_strata.end() ) ;
		}
		std::size_t const stratum_i = where->second ;
		m_sample_ranges[ stratum_i ] = impl::union_ranges( m_sample_ranges[ stratum_i ], range ) ;
	}
	
	std::string const& SampleStratification::stratum_name( std::size_t i ) const {
		assert( i < m_strata_names.size() ) ;
		return m_strata_names[i] ;
	}
	
	std::vector< SampleRange > const& SampleStratification::stratum( std::size_t i ) const {
		assert( i < m_strata_names.size() ) ;
		return m_sample_ranges[ i ] ;
	}

	std::vector< SampleRange > const& SampleStratification::stratum( std::string const& strata_name ) const {
		std::map< std::string, std::size_t >::const_iterator where = m_index_by_strata.find( strata_name ) ;
		if( where == m_index_by_strata.end() ) {
			throw genfile::BadArgumentError( "core::SampleStratification::stratum()", "strata_name=\"" + strata_name + "\"", "no such strata" ) ;
		}
		return m_sample_ranges[ where->second ] ;
	}
	
	std::ostream& operator<<( std::ostream& out, std::vector< SampleRange > const& range ) {
		for( std::size_t i = 0; i < range.size(); ++i ) {
			out << ((i>0) ? " " : "" ) << range[i] ;
		}
		return out ;
	}
	
}
