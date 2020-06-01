
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include "boost/optional.hpp"
#include "boost/noncopyable.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/CartesianProductVisitor.hpp"

namespace genfile {
	CartesianProductVisitor::CartesianProductVisitor( bool lower_triangle ):
		m_lower_triangle( lower_triangle ),
		m_min_distance(0)
	{}

	void CartesianProductVisitor::add_source(
		std::string const& name,
		genfile::SNPDataSource* source
	) {
		m_names.push_back( name ) ;
		m_sources.push_back( source ) ;
	}

	void CartesianProductVisitor::set_min_distance( int64_t distance ) {
		assert( distance >= 0 ) ;
		m_min_distance = distance ;
	}

	bool CartesianProductVisitor::step(  Callback callback ) {
		int64_t const distant = std::numeric_limits< int64_t >::max() ;
		std::vector< int > changed( m_sources.size(), 0 ) ;
		while( step_impl() ) {
			// keep track of cumulative changes
			for( std::size_t i = 0; i < changed.size(); ++i ) {
				changed[i] = std::max( changed[i], m_changed[i] ) ;
			}
			int64_t const min_distance = (m_min_distance == 0) ? distant : compute_min_distance( m_variants ) ;
			if( min_distance >= m_min_distance ) {
				callback(
					changed,
					m_variants,
					m_readers
				) ;
				return true ;
			}
		}
		return false ;
	}

	boost::optional< std::size_t > CartesianProductVisitor::count() const {
		boost::optional< std::size_t > result ;
		if( m_sources.size() > 0 && (result = m_sources[0]->size() )) {
			if( m_lower_triangle ) {
				assert( m_sources.size() == 2 ) ;
				result = ((*result) * (*result + 1) / 2) ;
			} else {
				for( std::size_t i = 1; i < m_sources.size(); ++i ) {
					boost::optional< std::size_t > source_size = m_sources[i]->size() ;
					if( !source_size ) {
						return boost::optional< std::size_t >() ;
					} else {
						(*result) *= (*source_size) ;
					}
				}
			}
		}
		return result ;
	}

	bool CartesianProductVisitor::step_impl() {
		if( m_sources.empty() ) {
			return false ;
		}
		if( m_changed.empty() ) {
			// First call.  Populate our data structures.
			m_changed.assign( m_sources.size(), 1 ) ;
			m_variants.resize( m_sources.size() ) ;
			assert( m_readers.size() == 0 ) ;
			for( std::size_t i = 0; i < m_sources.size(); ++i ) {
				if( !m_sources[i]->get_snp_identifying_data( &m_variants[i] )) {
					return false ;
				} ;
				m_readers.push_back( m_sources[i]->read_variant_data() ) ;
			}
		} else {
			// Move to the next variant.
			// Use a linearised recursion to step through all combinations.
			// We start at the last source and step to the next variant.
			// If a source is exhausted we reset it to the start and recurse to
			// the preceding source.
			std::size_t i = m_sources.size() - 1 ;
			for( ; true; --i ) {
				if( !m_sources[i]->get_snp_identifying_data( &m_variants[i] )) {
					if( i == 0 ) {
						return false ;
					} else {
						reset_source(i) ;
						m_sources[i]->get_snp_identifying_data( &m_variants[i] ) ;
					}
				} else {
					break ;
				}
			}
			m_changed.assign( m_sources.size(), 0 ) ;
			for( std::size_t j = i; j < m_sources.size(); ++j ) {
				m_changed[j] = 1 ;
				m_readers[j] = m_sources[j]->read_variant_data() ;
			}
		}
		return true ;
	}

	void CartesianProductVisitor::reset_source( std::size_t i ) {
		m_sources[i]->reset_to_start() ;
		if( m_lower_triangle ) {
			assert( m_sources.size() == 2 ) ;
			if( i == 1 ) {
				genfile::VariantIdentifyingData variant ;
				while( m_sources[1]->number_of_snps_read() < m_sources[0]->number_of_snps_read() ) {
					m_sources[1]->get_snp_identifying_data( &variant ) ;
					m_sources[1]->ignore_snp_probability_data() ;
				}
			}
		}
	}

	// return minimum physical distance between any pair of variants on same chromosome,
	// or maximum +pve value of int64 if all on different chromosomes
	int64_t CartesianProductVisitor::compute_min_distance( std::vector< genfile::VariantIdentifyingData > const& variants ) const {
		int64_t const distant = std::numeric_limits< int64_t >::max() ;
		int64_t result = distant ;
		for( std::size_t i = 0; i < ( variants.size()-1 ); ++i ) {
			for( std::size_t j = i+1; j < variants.size(); ++j ) {
				genfile::GenomePosition const& pos1 = variants[i].get_position() ;
				genfile::GenomePosition const& pos2 = variants[j].get_position() ;
				int64_t distance = (
					(pos1.chromosome() != pos2.chromosome())
					?
					distant
					:
					std::abs( int64_t(pos1.position()) - int64_t(pos2.position()))
				) ;
				result = std::min( result, distance ) ;
			}
		}
		return result ;
	}
}

