
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_CHROMOSOME_HPP
#define GENFILE_CHROMOSOME_HPP

#include <string>
#include <iostream>
#include <boost/noncopyable.hpp>
#include <boost/optional.hpp>
#include <boost/unordered_map.hpp>

namespace genfile {
	struct Chromosome ;
	namespace impl {
		struct Genome: public boost::noncopyable {
			virtual ~Genome() {} ;
			virtual bool compare( Chromosome const& left, Chromosome const& right ) const = 0 ;
		} ;

		struct MapGenome: public Genome {
			MapGenome( char const** strings, std::size_t N ) ;
			~MapGenome() ;
			bool compare( Chromosome const& left, Chromosome const& right ) const ;
		private:
			typedef boost::unordered_map< std::string, int > MapType ;
			MapType m_map ;
		} ;
		
		extern MapGenome human ;
	}
	
	struct Chromosome
	{
	public:
		friend struct impl::MapGenome ;

	public:
		Chromosome( impl::Genome const& genome = impl::human ):
			m_genome( &genome ),
			m_repr()
		{}

		Chromosome( std::string const& chromosome_str, impl::Genome const& genome = impl::human ):
			m_genome( &genome ),
			m_repr( chromosome_str )
		{
		}

		Chromosome( Chromosome const& other ):
			m_genome( other.m_genome ),
			m_repr( other.m_repr )
		{}

		Chromosome& operator=( Chromosome const& other ) {
			m_genome = other.m_genome ;
			m_repr = other.m_repr ;
			return *this ;
		}

		bool operator==( std::string const& other ) const {
			return m_repr == other ;
		}

		bool operator==( Chromosome const& other ) const {
			return m_repr == other.m_repr ;
		}

		bool operator!=( std::string const& other ) const {
			return m_repr != other ;
		}
		
		bool operator!=( Chromosome const& other ) const {
			return m_repr != other.m_repr ;
		}
		
		bool operator<=( Chromosome const& other ) const {
			return !m_genome->compare( other, *this ) ;
		}

		bool operator>=( Chromosome const& other ) const {
			return !m_genome->compare( *this, other ) ;
		}

		bool operator<( Chromosome const& other ) const {
			return m_genome->compare( *this, other ) ;
		}

		bool operator>( Chromosome const& other ) const {
			return m_genome->compare( other, *this ) ;
		}
		
		bool is_sex_determining() const ;
		bool is_autosome() const ;
		bool is_missing() const ;
		
		operator std::string () const ;

	private:
		impl::Genome const* m_genome ;
		boost::optional< std::string > m_repr ;
	} ;
	
	std::ostream& operator<<( std::ostream&, Chromosome const& ) ;
	std::istream& operator>>( std::istream&, Chromosome& ) ;


	struct ChromosomeMismatchError: public std::exception
	{
		ChromosomeMismatchError( Chromosome expected, Chromosome got )
			: m_expected( expected ),
			m_got( got )
		{}
		
		virtual ~ChromosomeMismatchError() throw() {}
		char const* what() const throw() { return "ChromosomeMismatchError" ; }

		Chromosome const& expected() const { return m_expected ; }
		Chromosome const& got() const { return m_got ; }
	private:
		Chromosome m_expected, m_got ;
	} ;
}

#endif
