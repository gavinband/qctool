#ifndef GENFILE_GET_SET_HPP
#define GENFILE_GET_SET_HPP

#include <string>
#include <vector>
#include "genfile/string_utils/string_utils.hpp"

namespace genfile {
	struct Ignorer
	{
		template< typename T > void operator()( T const& ) {}
		template< typename T > void operator()( T const&, T const& ) {}
		template< typename T > void operator()( T const&, T const&, T const& ) {}
		template< typename T1, typename T2 > void operator()( T1 const&, T2 const&, T2 const&, T2 const& ) {}
	}  ;

	Ignorer ignore() ;

	template< typename T >
	struct ValueSetter
	{
		ValueSetter( T& t ): m_t(t) {}
		template< typename T2 >
		void operator()( T2 const& t ) {
			m_t = T(t) ;
		}
	private:
		T& m_t ;
	} ;
	
	template<>
	struct ValueSetter< std::string >
	{
		ValueSetter< std::string >( std::string& t ): m_t(t) {}
		void operator()( char const c ) ;
		template< typename T2 > void operator()( T2 const& t ) {
			m_t = string_utils::to_string( t ) ;
		}
	private:
		std::string& m_t ;
	} ;

	template< typename T > ValueSetter< T > set_value( T& t ) {
		return ValueSetter< T >( t ) ;
	}

	template< typename T >
	struct GenotypeSetter
	{
		GenotypeSetter( T& t ): m_t( t ) {}
		void operator()( std::size_t i, double AA, double AB, double BB ) const {
			if( m_t.size() < i + 1 ) {
				m_t.resize( i + 1 ) ;
			}
			m_t[ i ][ 0 ] = AA ;
			m_t[ i ][ 1 ] = AB ;
			m_t[ i ][ 2 ] = BB ;
		}
	private:
		T& m_t ;
	} ;

	template<>
	struct GenotypeSetter< std::vector< double > >
	{
		GenotypeSetter( std::vector< double >& genotypes ) ;
		void operator()( std::size_t i, double, double, double ) const ;	

		private:
			std::vector< double >& m_genotypes ;
	} ;

	template< typename T > GenotypeSetter< T > set_genotypes( T& t ) {
		return GenotypeSetter< T >( t ) ;
	}
}

#endif
