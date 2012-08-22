
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
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPIdentifyingData2.hpp"

namespace genfile {
	namespace {
		std::size_t const MAX_SIZE = std::numeric_limits< uint32_t >::max() ;
	}

	SNPIdentifyingData2::SNPIdentifyingData2():
		m_rsid_start( 0 ),
		m_first_allele_start( 0 ),
		m_second_allele_start( 0 ),
		m_identifiers_start( 0 ),
		m_position( Chromosome(), 0 )
	{}

	SNPIdentifyingData2::SNPIdentifyingData2(
		std::string const& rsid,
		GenomePosition const& position,
		std::string const& first_allele,
		std::string const& second_allele
	):
		m_data( rsid + first_allele + second_allele ),
		m_rsid_start( 0 ),
		m_first_allele_start( m_rsid_start + rsid.size() ),
		m_second_allele_start( m_first_allele_start + first_allele.size() ),
		m_identifiers_start( m_second_allele_start + second_allele.size() ),
		m_position( position )
	{
		assert( rsid.size() < MAX_SIZE ) ; 
		assert( first_allele.size() < MAX_SIZE ); 
	}

	SNPIdentifyingData2::SNPIdentifyingData2(
		SNPIdentifyingData const& snp
	):
		m_data( snp.get_rsid() + snp.get_first_allele() + snp.get_second_allele() + snp.get_SNPID() ),
		m_rsid_start( 0 ),
		m_first_allele_start( m_rsid_start + snp.get_rsid().size() ),
		m_second_allele_start( m_first_allele_start + snp.get_first_allele().size() ),
		m_identifiers_start( m_second_allele_start + snp.get_second_allele().size() ),
		m_position( snp.get_position() )
	{}
	
	SNPIdentifyingData2::SNPIdentifyingData2( SNPIdentifyingData2 const& other ):
		m_data( other.m_data ),
		m_rsid_start( other.m_rsid_start ),
		m_first_allele_start( other.m_first_allele_start ),
		m_second_allele_start( other.m_second_allele_start ),
		m_identifiers_start( other.m_identifiers_start ),
		m_position( other.m_position )
	{}
	
	SNPIdentifyingData2& SNPIdentifyingData2::operator=( SNPIdentifyingData const& snp ) {
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
		m_first_allele_start = m_rsid_start + snp.get_rsid().size() ;
		m_second_allele_start = m_first_allele_start + snp.get_first_allele().size() ;
		m_identifiers_start = m_second_allele_start + snp.get_second_allele().size() ;
		m_position = snp.get_position() ;
		return *this ;
	}

	SNPIdentifyingData2& SNPIdentifyingData2::operator=( SNPIdentifyingData2 const& other ) {
		m_data = other.m_data ;
		m_rsid_start = other.m_rsid_start ;
		m_first_allele_start = other.m_first_allele_start ;
		m_second_allele_start = other.m_second_allele_start ;
		m_identifiers_start = other.m_identifiers_start ;
		m_position = other.m_position ;
		return *this ;
	}
	
	SNPIdentifyingData2::operator SNPIdentifyingData() const {
		return SNPIdentifyingData( m_data.substr( m_identifiers_start, m_data.size() - m_identifiers_start ), get_rsid(), get_position(), get_first_allele(), get_second_allele() ) ;
	}

	void SNPIdentifyingData2::set_rsid( string_utils::slice const& rsid ) {
		assert( rsid.size() + std::size_t( m_rsid_start ) < MAX_SIZE ) ;
		std::string new_data ;
		new_data.reserve( m_rsid_start + m_data.size() - m_first_allele_start ) ;
		new_data.append( m_data.begin(), m_data.begin() + m_rsid_start ) ;
		new_data.append( rsid.begin(), rsid.end() ) ;
		new_data.append( m_data.begin() + m_first_allele_start, m_data.end() ) ;
		m_data.swap( new_data ) ;

		std::size_t difference = rsid.size() - m_first_allele_start + m_rsid_start ;
		m_first_allele_start += difference ;
		m_second_allele_start += difference ;
	}

	void SNPIdentifyingData2::set_first_allele( string_utils::slice const& allele ) {
		assert( allele.size() + std::size_t( m_rsid_start ) + std::size_t( m_first_allele_start ) <  MAX_SIZE ) ;
		std::string new_data ;
		new_data.reserve( m_first_allele_start + m_data.size() - m_second_allele_start ) ;
		new_data.append( m_data.begin(), m_data.begin() + m_first_allele_start ) ;
		new_data.append( allele.begin(), allele.end() ) ;
		new_data.append( m_data.begin() + m_second_allele_start, m_data.end() ) ;
		m_data.swap( new_data ) ;
		std::size_t difference = allele.size() - m_second_allele_start + m_first_allele_start ;
		m_second_allele_start += difference ;
	}

	void SNPIdentifyingData2::set_second_allele( string_utils::slice const& allele ) {
		std::string new_data ;
		new_data.reserve( m_second_allele_start + allele.size() ) ;
		new_data.append( m_data.begin(), m_data.begin() + m_second_allele_start ) ;
		new_data.append( allele.begin(), allele.end() ) ;
		m_data.swap( new_data ) ;
	}

	void SNPIdentifyingData2::add_identifier( slice const& id ) {
		assert( id.size() > 0 ) ;
		if( id != get_rsid() ) {
			std::vector< slice > ids = get_identifiers() ;
			if( ids.size() > 0 ) {
				if( std::find( ids.begin(), ids.end(), id ) == ids.end() ) {
					m_data += "\t" + std::string( id ) ;
				}
			} else {
				m_data += std::string( id ) ;
			}
		}
	}

	std::vector< genfile::string_utils::slice > SNPIdentifyingData2::get_identifiers() const {
		if( m_identifiers_start == m_data.size() ) {
			return std::vector< slice >() ;
		}
		else {
			return slice( m_data, m_identifiers_start, m_data.size() ).split( "\t" ) ;
		}
	}

	void SNPIdentifyingData2::get_identifiers( boost::function< void( slice ) > callback ) const {
		slice( m_data, m_identifiers_start, m_data.size() ).split( "\t", callback ) ;
	}

	std::size_t SNPIdentifyingData2::get_estimated_bytes_used() const {
		return sizeof( SNPIdentifyingData2 )
			+ m_data.size()
			+ ( 3 * sizeof( std::size_t ) ) // there is about this much overhead on the stack per string.
		;
	}

	std::ostream& operator<<( std::ostream& out, SNPIdentifyingData2 const& data ) {
		out << data.get_rsid() ;
		std::vector< genfile::string_utils::slice > const ids = data.get_identifiers() ;
		if( ids.size() > 0 ) {
			out << " [" ;
			for( std::size_t i = 0; i < ids.size(); ++i ) {
				out << ( i > 0 ? "," : "" ) << ids[i] ;
			}
			out << "]" ;
		}
		out << " " << data.get_position().chromosome()
			<< " " << data.get_position().position()
			<< " " << data.get_first_allele()
			<< " " << data.get_second_allele() ;
		return out ;
	}

	std::ostream& operator<<( std::ostream& out, std::vector< SNPIdentifyingData2 > const& data ) {
		for( std::size_t i = 0; i < data.size(); ++i ) {
			out << data[i] << "\n" ;
		}
		return out ;
	}
	
	bool operator==( SNPIdentifyingData2 const& left, SNPIdentifyingData2 const& right ) {
		using genfile::string_utils::slice ;
		return left.m_data == right.m_data
			&& left.m_rsid_start == right.m_rsid_start
			&& left.m_first_allele_start == right.m_first_allele_start
			&& left.m_second_allele_start == right.m_second_allele_start
			&& left.m_identifiers_start == right.m_identifiers_start
			&& slice( left.m_data, left.m_identifiers_start, left.m_data.size() ) == slice( right.m_data, right.m_identifiers_start, right.m_data.size() )
			&& left.m_position == right.m_position ;
	}

	bool operator!=( SNPIdentifyingData2 const& left, SNPIdentifyingData2 const& right ) {
		return !( left == right ) ;
	}
	
	
    bool operator<( SNPIdentifyingData2 const& left, SNPIdentifyingData2 const& right ) {
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
							(left.get_first_allele() < right.get_first_allele())
							||
							(
								(left.get_first_allele() == right.get_first_allele())
								&&
								(
									(left.get_second_allele() < right.get_second_allele())
									||
									(
										(left.get_second_allele() == right.get_second_allele())
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

	SNPIdentifyingData2::CompareFields::CompareFields( std::string const& fields_to_compare ):
		m_fields_to_compare( parse_fields_to_compare( fields_to_compare ) )
	{
		assert( !m_fields_to_compare.empty() ) ;
	}
	
	SNPIdentifyingData2::CompareFields::CompareFields( CompareFields const& other ):
		m_fields_to_compare( other.m_fields_to_compare )
	{}

	std::vector< int > SNPIdentifyingData2::CompareFields::parse_fields_to_compare( std::string const& field_spec ) {
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
			else if( elts[i] == "IDs" ) {
				result.push_back( eIDs ) ;
			}
			else {
				assert(0) ;
			}
		}
		return result ;
	}

	bool SNPIdentifyingData2::CompareFields::operator()( SNPIdentifyingData2 const& left, SNPIdentifyingData2 const& right ) const {
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
					if( left.get_first_allele() > right.get_first_allele() ) {
						return false ;
					}
					else if( left.get_first_allele() < right.get_first_allele() ) {
						return true ;
					}
					else if( left.get_second_allele() > right.get_second_allele() ) {
						return false ;
					}
					else if( left.get_second_allele() < right.get_second_allele() ) {
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
	
	bool SNPIdentifyingData2::CompareFields::are_equal( SNPIdentifyingData2 const& left, SNPIdentifyingData2 const& right ) const {
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
					if( left.get_first_allele() != right.get_first_allele() ) {
						return false ;
					}
					if( left.get_second_allele() != right.get_second_allele() ) {
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
	
	bool SNPIdentifyingData2::CompareFields::check_if_comparable_fields_are_known( SNPIdentifyingData2 const& value ) const {
		if( std::find( m_fields_to_compare.begin(), m_fields_to_compare.end(), int( eAlleles ) ) != m_fields_to_compare.end() ) {
			return value.get_first_allele() != "?" && value.get_second_allele() != "?" ;
		}
		else {
			return true ;
		}
	}
}
