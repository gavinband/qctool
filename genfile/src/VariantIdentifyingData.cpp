
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <cassert>
#include "genfile/GenomePosition.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantIdentifyingData.hpp"

namespace genfile {
	namespace {
		std::size_t const MAX_SIZE = std::numeric_limits< uint32_t >::max() ;
	}

	VariantIdentifyingData::VariantIdentifyingData():
		m_rsid_start( 0 ),
		m_allele_starts(),
		m_identifiers_start( 0 ),
		m_position( Chromosome(), 0 )
	{}

	namespace {
		std::vector< uint32_t > compute_allele_starts( uint32_t const rsid_start, std::string const& rsid, std::string const& first_allele, std::string const& second_allele ) {
			std::vector< uint32_t > result( 2 ) ;
			result[0] = rsid_start + rsid.size() ;
			result[1] = result[0] + first_allele.size() ;
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
		m_identifiers_start( m_allele_starts[1] + second_allele.size() ),
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
		m_identifiers_start( m_allele_starts[1] + second_allele.size() ),
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
		m_identifiers_start( m_allele_starts[1] + snp.get_second_allele().size() ),
		m_position( snp.get_position() )
	{}
	
	VariantIdentifyingData::VariantIdentifyingData( VariantIdentifyingData const& other ):
		m_data( other.m_data ),
		m_rsid_start( other.m_rsid_start ),
		m_allele_starts( other.m_allele_starts ),
		m_identifiers_start( other.m_identifiers_start ),
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
		m_allele_starts.resize(2) ;
		m_allele_starts[0] = m_rsid_start + snp.get_rsid().size() ;
		m_allele_starts[1] = m_allele_starts[0] + snp.get_first_allele().size() ;
		m_identifiers_start = m_allele_starts[1] + snp.get_second_allele().size() ;
		m_position = snp.get_position() ;
		return *this ;
	}

	VariantIdentifyingData& VariantIdentifyingData::operator=( VariantIdentifyingData const& other ) {
		m_data = other.m_data ;
		m_rsid_start = other.m_rsid_start ;
		m_allele_starts = other.m_allele_starts ;
		m_identifiers_start = other.m_identifiers_start ;
		m_position = other.m_position ;
		return *this ;
	}
	
	VariantIdentifyingData::operator SNPIdentifyingData() const {
		assert( m_allele_starts.size() == 2 ) ;
		return SNPIdentifyingData(
			m_data.substr( m_identifiers_start, m_data.size() - m_identifiers_start ),
			get_rsid(),
			get_position(),
			get_allele(0),
			get_allele(1)
		) ;
	}

	void VariantIdentifyingData::set_rsid( string_utils::slice const& rsid ) {
		std::string new_data ;
		std::size_t const difference = rsid.size() - m_allele_starts[0] ;
		new_data.reserve( m_data.size() + difference ) ;
		new_data.append( rsid ) ;
		new_data.append( m_data.begin() + m_allele_starts[0], m_data.end() ) ;
		m_data.swap( new_data ) ;
		for( std::size_t i = 0; i < m_allele_starts.size(); ++i ) {
			m_allele_starts[i] += difference ;
		}
		m_identifiers_start += difference ;
	}

	void VariantIdentifyingData::set_first_allele( string_utils::slice const& allele ) {
		std::string new_data ;
		std::size_t const difference = allele.size() - ( m_allele_starts[1] - m_allele_starts[0] ) ;
		new_data.reserve( m_data.size() + difference ) ;
		new_data.append( m_data.begin(), m_data.begin() + m_allele_starts[0] ) ;
		new_data.append( allele.begin(), allele.end() ) ;
		new_data.append( m_data.begin() + m_allele_starts[1], m_data.end() ) ;
		m_data.swap( new_data ) ;
		for( std::size_t i = 1; i < m_allele_starts.size(); ++i ) {
			m_allele_starts[i] += difference ;
		}
		m_identifiers_start += difference ;
	}

	void VariantIdentifyingData::set_second_allele( string_utils::slice const& allele ) {
		std::string new_data ;
		std::size_t const difference = allele.size() - ( m_identifiers_start - m_allele_starts[1] ) ;
		new_data.reserve( m_data.size() + difference ) ;
		new_data.append( m_data.begin(), m_data.begin() + m_allele_starts[1] ) ;
		new_data.append( allele.begin(), allele.end() ) ;
		new_data.append( m_data.begin() + m_identifiers_start, m_data.end() ) ;
		m_data.swap( new_data ) ;
		m_identifiers_start += difference ;
	}

	void VariantIdentifyingData::add_allele( slice const& allele ) {
		std::string new_data ;
		new_data.reserve( m_data.size() + allele.size() ) ;
		new_data.append( m_data.begin(), m_data.begin() + m_identifiers_start ) ;
		new_data.append( allele.begin(), allele.end() ) ;
		new_data.append( m_data.begin() + m_identifiers_start, m_data.end() ) ;
		m_identifiers_start += allele.size() ;
	}

	void VariantIdentifyingData::get_alleles( boost::function< void( slice ) > callback ) const {
		if( m_allele_starts.size() > 0 ) {
			for( std::size_t i = 0; (i+1) < m_allele_starts.size(); ++i ) {
				callback( slice( m_data, m_allele_starts[i], m_allele_starts[i+1] )) ;
			}
			callback( slice( m_data, m_allele_starts.back(), m_identifiers_start )) ;
		}
	}

	void VariantIdentifyingData::swap_alleles() {
		assert( m_allele_starts.size() == 2 ) ;
		std::string const first_allele = get_allele(0) ;
		genfile::string_utils::slice const second_allele = get_allele(1) ;
		std::copy( second_allele.begin(), second_allele.end(), m_data.begin() + m_allele_starts[0] ) ;
		m_allele_starts[1] = m_allele_starts[0] + second_allele.size() ;
		std::copy( first_allele.begin(), first_allele.end(), m_data.begin() + m_allele_starts[1] ) ;
	}

	void VariantIdentifyingData::clear_identifiers() {
		m_data.resize( m_identifiers_start ) ;
	}

	void VariantIdentifyingData::add_identifier( slice const& id ) {
		// Deal with strange non-IDs
		assert( id.size() > 0 ) ;
		if( id != get_rsid() ) {
			std::vector< slice > ids = get_alternative_identifiers() ;
			if( ids.size() > 0 ) {
				if( std::find( ids.begin(), ids.end(), id ) == ids.end() ) {
					m_data += "\t" + std::string( id ) ;
				}
			} else {
				m_data += std::string( id ) ;
			}
		}
	}

	std::vector< genfile::string_utils::slice > VariantIdentifyingData::get_alternative_identifiers() const {
		if( m_identifiers_start == m_data.size() ) {
			return std::vector< slice >() ;
		}
		else {
			return slice( m_data, m_identifiers_start, m_data.size() ).split( "\t" ) ;
		}
	}

	void VariantIdentifyingData::get_alternative_identifiers( boost::function< void( slice ) > callback ) const {
		if( m_identifiers_start == m_data.size() ) {
			return ;
		}
		slice( m_data, m_identifiers_start, m_data.size() ).split( "\t", callback ) ;
	}

	std::string VariantIdentifyingData::get_alternate_identifiers_as_string() const {
		return join( get_alternative_identifiers(), "," ) ;
	}

	std::size_t VariantIdentifyingData::get_estimated_bytes_used() const {
		return sizeof( VariantIdentifyingData )
			+ m_data.size()
			+ ( 3 * sizeof( std::size_t ) ) // there is about this much overhead on the stack per string.
		;
	}

	std::ostream& operator<<( std::ostream& out, VariantIdentifyingData const& data ) {
		out << data.get_rsid() ;
		std::vector< genfile::string_utils::slice > const ids = data.get_alternative_identifiers() ;
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
			&& left.m_identifiers_start == right.m_identifiers_start
			&& slice( left.m_data, left.m_identifiers_start, left.m_data.size() ) == slice( right.m_data, right.m_identifiers_start, right.m_data.size() )
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
										( slice( left.m_data, left.m_identifiers_start, left.m_data.size() ) < slice( right.m_data, right.m_identifiers_start, right.m_data.size() ) )
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
					if( left.get_alternative_identifiers() > right.get_alternative_identifiers() ) {
						return false ;
					}
					else if( left.get_alternative_identifiers() < right.get_alternative_identifiers() ) {
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
					if( left.get_alternative_identifiers() != right.get_alternative_identifiers() ) {
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
