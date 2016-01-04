
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include <boost/bind.hpp>
#include "genfile/GenomePosition.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantIdentifyingData.hpp"

namespace genfile {
	namespace {
		std::size_t const MAX_SIZE = std::numeric_limits< uint32_t >::max() ;
	}

	VariantIdentifyingData::VariantIdentifyingData( std::string const& rsid ):
		m_data( rsid ),
		m_rsid_start( 0 ),
		m_allele_starts( 1, rsid.size() ),
		m_position( Chromosome(), 0 )
	{}

	namespace {
		std::vector< uint32_t > compute_allele_starts(
			uint32_t const rsid_start,
			std::string const& rsid,
			std::string const& first_allele,
			std::string const& second_allele
		) {
			std::vector< uint32_t > result( 3 ) ;
			result[0] = rsid_start + rsid.size() ;
			result[1] = result[0] + first_allele.size() ;
			result[2] = result[1] + second_allele.size() ;
			return result ;
		}
	}

	VariantIdentifyingData::VariantIdentifyingData(
		std::string const& rsid,
		GenomePosition const& position,
		std::string const& first_allele,
		std::string const& second_allele
	):
		m_data( rsid + first_allele + second_allele ),
		m_rsid_start( 0 ),
		m_allele_starts( compute_allele_starts( m_rsid_start, rsid, first_allele, second_allele ) ),
		m_position( position )
	{
		assert( rsid.size() < MAX_SIZE ) ; 
		assert( first_allele.size() < MAX_SIZE ); 
	}

 	VariantIdentifyingData::VariantIdentifyingData(
		std::string const& SNPID,
		std::string const& rsid,
		GenomePosition const& position,
		std::string const& first_allele,
		std::string const& second_allele
	):
		m_data( rsid + first_allele + second_allele + SNPID ),
		m_rsid_start( 0 ),
		m_allele_starts( compute_allele_starts( m_rsid_start, rsid, first_allele, second_allele ) ),
		m_position( position )
	{
		assert( rsid.size() < MAX_SIZE ) ;
		assert( SNPID.size() < MAX_SIZE ) ;
		assert( first_allele.size() < MAX_SIZE );
	}

	VariantIdentifyingData::VariantIdentifyingData(
		SNPIdentifyingData const& snp
	):
		m_data( snp.get_rsid() + snp.get_first_allele() + snp.get_second_allele() + snp.get_SNPID() ),
		m_rsid_start( 0 ),
		m_allele_starts( compute_allele_starts( m_rsid_start, snp.get_rsid(), snp.get_first_allele(), snp.get_second_allele() ) ),
		m_position( snp.get_position() )
	{}
	
	VariantIdentifyingData::VariantIdentifyingData( VariantIdentifyingData const& other ):
		m_data( other.m_data ),
		m_rsid_start( other.m_rsid_start ),
		m_allele_starts( other.m_allele_starts ),
		m_position( other.m_position )
	{}
	
	VariantIdentifyingData& VariantIdentifyingData::operator=( SNPIdentifyingData const& snp ) {
		assert( snp.get_SNPID().size() < MAX_SIZE ) ; 
		assert( snp.get_rsid().size() < MAX_SIZE ) ; 
		assert( snp.get_first_allele().size() < MAX_SIZE ) ; 
		assert( snp.get_second_allele().size() < MAX_SIZE ) ; 

		bool const use_SNPID = ( snp.get_SNPID() != snp.get_rsid() ) ;
		if( use_SNPID ) {
			std::string new_data = snp.get_rsid() + snp.get_first_allele() + snp.get_second_allele() + snp.get_SNPID() ;
			m_data.swap( new_data ) ;
		} else {
			std::string new_data = snp.get_rsid() + snp.get_first_allele() + snp.get_second_allele() ;
			m_data.swap( new_data ) ;
		}
		
		m_rsid_start = 0 ;
		m_allele_starts.resize(3) ;
		m_allele_starts[0] = m_rsid_start + snp.get_rsid().size() ;
		m_allele_starts[1] = m_allele_starts[0] + snp.get_first_allele().size() ;
		m_allele_starts[2] = m_allele_starts[1] + snp.get_second_allele().size() ;
		m_position = snp.get_position() ;
		return *this ;
	}

	VariantIdentifyingData& VariantIdentifyingData::operator=( VariantIdentifyingData const& other ) {
		m_data = other.m_data ;
		m_rsid_start = other.m_rsid_start ;
		m_allele_starts = other.m_allele_starts ;
		m_position = other.m_position ;
		return *this ;
	}
	
	VariantIdentifyingData::operator SNPIdentifyingData() const {
		assert( m_allele_starts.size() == 2 ) ;
		return SNPIdentifyingData(
			m_data.substr( m_allele_starts.back(), m_data.size() - m_allele_starts.back() ),
			get_rsid(),
			get_position(),
			get_allele(0),
			get_allele(1)
		) ;
	}

	void VariantIdentifyingData::set_primary_id( string_utils::slice const& rsid ) {
		std::string new_data ;
		std::size_t current_size = ( m_allele_starts.size() > 0 ) ? m_allele_starts[0] : 0 ;
		std::size_t const difference = rsid.size() - current_size ;
		new_data.reserve( m_data.size() + difference ) ;
		new_data.append( rsid ) ;
		new_data.append( m_data.begin() + current_size, m_data.end() ) ;
		m_data.swap( new_data ) ;
		for( std::size_t i = 0; i < m_allele_starts.size(); ++i ) {
			m_allele_starts[i] += difference ;
		}
	}
	
	void VariantIdentifyingData::set_allele( std::size_t i, slice const& allele ) {
		assert( i < m_allele_starts.size() ) ;
		std::size_t const old_end = m_allele_starts[i+1] ;
		std::size_t old_allele_size = old_end - m_allele_starts[i] ;
		if( old_allele_size == allele.size() ) {
			std::copy( &allele[0], &allele[0] + allele.size(), m_data.begin() + m_allele_starts[i] ) ;
		} else {
			std::string new_data ;
			new_data.reserve( (m_data.size() + allele.size()) - old_allele_size ) ;
			new_data.append( m_data.begin(), m_data.begin() + m_allele_starts[i] ) ;
			new_data.append( allele.begin(), allele.end() ) ;
			new_data.append( m_data.begin() + old_end, m_data.end() ) ;
			if( allele.size() > old_allele_size ) { // handle unsigned arithmetic.  Is this needed?
				for( std::size_t j = i+1; j < m_allele_starts.size(); ++j ) {
					m_allele_starts[j] += allele.size() - old_allele_size ;
				}
			} else {
				for( std::size_t j = i+1; j < m_allele_starts.size(); ++j ) {
					m_allele_starts[j] -= old_allele_size - allele.size() ;
				}
			}
			m_data.swap( new_data ) ;
		}
	}

	void VariantIdentifyingData::set_allele( std::size_t i, std::string const& allele ) {
		set_allele( i, slice( allele )) ;
	}

	void VariantIdentifyingData::set_allele( std::size_t i, char const* allele ) {
		set_allele( i, slice( allele )) ;
	}

	void VariantIdentifyingData::add_allele( slice const& allele ) {
		std::string new_data ;
		new_data.reserve( m_data.size() + allele.size() ) ;
		new_data.append( m_data.begin(), m_data.begin() + m_allele_starts.back() ) ;
		new_data.append( allele.begin(), allele.end() ) ;
		new_data.append( m_data.begin() + m_allele_starts.back(), m_data.end() ) ;
		m_data.swap( new_data ) ;
		m_allele_starts.push_back( m_allele_starts.back() + allele.size() ) ;
	}

	string_utils::slice VariantIdentifyingData::get_allele( std::size_t i ) const {
		assert( i+1 < m_allele_starts.size() ) ;
		return slice( m_data, m_allele_starts[i], m_allele_starts[i+1] ) ;
	}

	std::vector< string_utils::slice > VariantIdentifyingData::get_alleles(
		std::size_t start,
		std::size_t end
	) const {
		std::vector< slice > result ;
		end = std::min( end, (m_allele_starts.size()-1) ) ;
		for( std::size_t i = start; i < end; ++i ) {
			result.push_back( slice( m_data, m_allele_starts[i], m_allele_starts[i+1] )) ;
		}
		return result ;
	}

	void VariantIdentifyingData::get_alleles(
		boost::function< void( slice ) > callback,
		std::size_t start,
		std::size_t end
	) const {
		end = std::min( end, (m_allele_starts.size()-1) ) ;
		for( std::size_t i = start; i < end; ++i ) {
			callback( slice( m_data, m_allele_starts[i], m_allele_starts[i+1] )) ;
		}
	}

	void VariantIdentifyingData::swap_alleles() {
		assert( m_allele_starts.size() == 3 ) ;
		std::string const first_allele = get_allele(0) ;
		genfile::string_utils::slice const second_allele = get_allele(1) ;
		std::copy( second_allele.begin(), second_allele.end(), m_data.begin() + m_allele_starts[0] ) ;
		m_allele_starts[1] = m_allele_starts[0] + second_allele.size() ;
		std::copy( first_allele.begin(), first_allele.end(), m_data.begin() + m_allele_starts[1] ) ;
	}

	void VariantIdentifyingData::clear_identifiers() {
		m_data.resize( m_allele_starts.back() ) ;
	}

	void VariantIdentifyingData::add_identifier( slice const& id ) {
		assert( id.size() > 0 ) ;
		if( id != get_rsid() ) {
			std::vector< slice > const& ids = get_identifiers(1) ;
			if( ids.size() > 0 ) {
				if( std::find( ids.begin(), ids.end(), id ) == ids.end() ) {
					m_data += "\t" + std::string( id ) ;
				}
			} else {
				m_data += std::string( id ) ;
			}
		}
	}

	namespace {
		void push_back( std::vector< genfile::string_utils::slice >* target, genfile::string_utils::slice const& value ) {
			target->push_back( value ) ;
		}
	}

	std::size_t VariantIdentifyingData::number_of_identifiers() const {
		return get_identifiers().size() ;
	}
	std::vector< genfile::string_utils::slice > VariantIdentifyingData::get_identifiers(
		std::size_t start,
		std::size_t end
	) const {
		std::vector< genfile::string_utils::slice > result ;
		get_identifiers( boost::bind( &push_back, &result, _1 ), start, end ) ;
		return result ;
	}

	void VariantIdentifyingData::get_identifiers(
		boost::function< void( slice ) > callback,
		std::size_t start,
		std::size_t end
	) const {
		assert( start <= end ) ;
		std::vector< genfile::string_utils::slice > elts( 1, get_rsid() ) ;
		if( m_allele_starts.back() < m_data.size() ) {
			slice( m_data, m_allele_starts.back(), m_data.size() ).split( "\t", boost::bind( &push_back, &elts, _1 )) ;
		}
		assert( end == std::string::npos || end <= elts.size() ) ;
		end = std::min( end, elts.size() ) ;
		for( std::size_t i = start; i < end; ++i ) {
			callback( elts[i] ) ;
		}
	}

	std::string VariantIdentifyingData::get_identifiers_as_string(
		std::string const& separator,
		std::size_t start,
		std::size_t end
	) const {
		return join( get_identifiers( start, end ), separator ) ;
	}

	std::size_t VariantIdentifyingData::estimate_bytes_used() const {
		return sizeof( VariantIdentifyingData )
			+ m_data.size()
			+ ( 3 * sizeof( std::size_t ) ) // there is about this much overhead on the stack per string.
		;
	}

	std::ostream& operator<<( std::ostream& out, VariantIdentifyingData const& data ) {
		out << data.get_rsid() ;
		std::vector< genfile::string_utils::slice > const ids = data.get_identifiers( 1 ) ;
		if( ids.size() > 0 ) {
			out << " [" ;
			for( std::size_t i = 0; i < ids.size(); ++i ) {
				out << ( i > 0 ? "," : "" ) << ids[i] ;
			}
			out << "]" ;
		}
		out << " " << data.get_position().chromosome()
			<< " " << data.get_position().position()
			<< " " << data.get_allele(0)
			<< " " << data.get_allele(1) ;
		return out ;
	}

	std::ostream& operator<<( std::ostream& out, std::vector< VariantIdentifyingData > const& data ) {
		for( std::size_t i = 0; i < data.size(); ++i ) {
			out << data[i] << "\n" ;
		}
		return out ;
	}
	
	bool operator==( VariantIdentifyingData const& left, VariantIdentifyingData const& right ) {
		using genfile::string_utils::slice ;
		return left.m_data == right.m_data
			&& left.m_rsid_start == right.m_rsid_start
			&& left.m_allele_starts == right.m_allele_starts
			&& slice( left.m_data, left.m_allele_starts.back(), left.m_data.size() ) == slice( right.m_data, right.m_allele_starts.back(), right.m_data.size() )
			&& left.m_position == right.m_position ;
	}

	bool operator!=( VariantIdentifyingData const& left, VariantIdentifyingData const& right ) {
		return !( left == right ) ;
	}
	
	
    bool operator<( VariantIdentifyingData const& left, VariantIdentifyingData const& right ) {
		using genfile::string_utils::slice ;
		return (
			( left.get_position() < right.get_position() )
			||
			(
				(left.get_position() == right.get_position())
				&&
				(
					( left.get_rsid() < right.get_rsid() )
					||
					(
						(left.get_rsid() == right.get_rsid())
						&&
						(
							(left.get_allele(0) < right.get_allele(0))
							||
							(
								(left.get_allele(0) == right.get_allele(0))
								&&
								(
									(left.get_allele(1) < right.get_allele(1))
									||
									(
										(left.get_allele(1) == right.get_allele(1))
										&&
										(
											slice( left.m_data, left.m_allele_starts.back(), left.m_data.size() )
												< slice( right.m_data, right.m_allele_starts.back(), right.m_data.size() )
										)
									)
								)
							)
						)
					)
				)
			)
		) ;
	}

	VariantIdentifyingData::CompareFields::CompareFields(
		std::string const& fields_to_compare,
		bool flip_alleleles_if_necessary
	):
		m_fields_to_compare( parse_fields_to_compare( fields_to_compare ) ),
		m_flip_alleles_if_necessary( flip_alleleles_if_necessary )
	{
		assert( !m_fields_to_compare.empty() ) ;
	}
	
	VariantIdentifyingData::CompareFields::CompareFields( CompareFields const& other ):
		m_fields_to_compare( other.m_fields_to_compare ),
		m_flip_alleles_if_necessary( other.m_flip_alleles_if_necessary )
	{}

	VariantIdentifyingData::CompareFields& VariantIdentifyingData::CompareFields::operator=(
		VariantIdentifyingData::CompareFields const& other
	) {
		m_fields_to_compare = other.m_fields_to_compare ;
		m_flip_alleles_if_necessary = other.m_flip_alleles_if_necessary ;
		return *this ;
	}

	std::vector< int > VariantIdentifyingData::CompareFields::parse_fields_to_compare( std::string const& field_spec ) {
		std::vector< std::string > elts = string_utils::split_and_strip( field_spec, ",", " \t\n\r" ) ;
		assert( elts.size() > 0 ) ;
		std::vector< int > result ;
		for( std::size_t i = 0; i < elts.size(); ++i ) {
			if( elts[i] == "rsid" ) {
				result.push_back( eRSID ) ;
			}
			else if( elts[i] == "position" ) {
				result.push_back( ePosition ) ;
			}
			else if( elts[i] == "alleles" ) {
				result.push_back( eAlleles ) ;
			}
			else if( elts[i] == "IDs" || elts[i] == "SNPID" ) {
				result.push_back( eIDs ) ;
			}
			else {
				assert(0) ;
			}
		}
		return result ;
	}

	bool VariantIdentifyingData::CompareFields::operator()( VariantIdentifyingData const& left, VariantIdentifyingData const& right ) const {
		// Lexicographical compare in the same order as operator< above, but using only the specified fields.

		for( std::size_t i = 0; i < m_fields_to_compare.size(); ++i ) {
			switch( m_fields_to_compare[i] ) {
				case ePosition:
					if( left.get_position() > right.get_position() ) {
						return false ;
					}
					else if( left.get_position() < right.get_position() ) {
						return true ;
					}
					break ;
				case eRSID:
					if( left.get_rsid() > right.get_rsid() ) {
						return false ;
					}
					else if( left.get_rsid() < right.get_rsid() ) {
						return true ;
					}
					break ;
				case eIDs:
					if( left.get_identifiers() > right.get_identifiers() ) {
						return false ;
					}
					else if( left.get_identifiers() < right.get_identifiers() ) {
						return true ;
					}
					break ;
				case eAlleles:
					if( left.get_allele(0) > right.get_allele(0) ) {
						return false ;
					}
					else if( left.get_allele(0) < right.get_allele(0) ) {
						return true ;
					}
					else if( left.get_allele(1) > right.get_allele(1) ) {
						return false ;
					}
					else if( left.get_allele(1) < right.get_allele(1) ) {
						return true ;
					}
					break ;
				default:
					assert(0) ;
					break ;
			}
		}
		return false ;
	}
	
	bool VariantIdentifyingData::CompareFields::are_equal( VariantIdentifyingData const& left, VariantIdentifyingData const& right ) const {
		for( std::size_t i = 0; i < m_fields_to_compare.size(); ++i ) {
			switch( m_fields_to_compare[i] ) {
				case ePosition:
					if( left.get_position() != right.get_position() ) {
						return false ;
					}
					break ;
				case eRSID:
					if( left.get_rsid() != right.get_rsid() ) {
						return false ;
					}
					break ;
				case eIDs:
					if( left.get_identifiers() != right.get_identifiers() ) {
						return false ;
					}
					break ;
				case eAlleles:
					if( left.get_allele(0) != right.get_allele(0) ) {
						return false ;
					}
					if( left.get_allele(1) != right.get_allele(1) ) {
						return false ;
					}
					break ;
				default:
					assert(0) ;
					break ;
			}
		}
		return true ;
	}
	
	std::string VariantIdentifyingData::CompareFields::get_summary() const {
		std::string result = "comparing " ;
		for( std::size_t i = 0; i < m_fields_to_compare.size(); ++i ) {
			if( i > 0 ) {
				result += "," ;
			}
			switch( m_fields_to_compare[i] ) {
				case eIDs:
				result += "ids" ;
				break ;
				case eRSID:
				result += "rsid" ;
				break ;
				case ePosition:
				result += "position" ;
				break ;
				case eAlleles:
				result += "alleles" ;
				break ;
			}
		}
		return result ;
	}
	
	bool VariantIdentifyingData::CompareFields::check_if_comparable_fields_are_known( VariantIdentifyingData const& value ) const {
		if( std::find( m_fields_to_compare.begin(), m_fields_to_compare.end(), int( eAlleles ) ) != m_fields_to_compare.end() ) {
			return value.get_allele(0) != "?" && value.get_allele(1) != "?" ;
		}
		else {
			return true ;
		}
	}
}
