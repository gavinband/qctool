
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_CHROMOSOME_HPP
#define GENFILE_CHROMOSOME_HPP

#include <string>
#include <iostream>
#include <boost/optional.hpp>

namespace genfile {
	struct Genome {
	} ;

	struct Chromosome
	{
		Chromosome():
			m_repr()
		{}

		Chromosome( std::string const& chromosome_str ):
			m_repr( chromosome_str )
		{}

		Chromosome( Chromosome const& other ):
			m_repr( other.m_repr )
		{}

		Chromosome& operator=( Chromosome const& other ) {
			m_repr = other.m_repr ;
			return *this ;
		}

		Chromosome& operator=( std::string const& other ) {
			m_repr = other ;
			return *this ;
		}
		
		bool operator==( std::string const& other ) const {
			return m_repr == other ;
		}

		bool operator==( Chromosome other ) const {
			return m_repr == other.m_repr ;
		}

		bool operator!=( std::string const& other ) const {
			return m_repr != other ;
		}
		
		bool operator!=( Chromosome other ) const {
			return m_repr != other.m_repr ;
		}
		
		bool operator<=( Chromosome other ) const {
			return m_repr <= other.m_repr ;
		}

		bool operator>=( Chromosome other ) const {
			return m_repr >= other.m_repr ;
		}

		bool operator<( Chromosome other ) const {
			return m_repr < other.m_repr ;
		}

		bool operator>( Chromosome other ) const {
			return m_repr > other.m_repr ;
		}
		
		bool is_sex_determining() const ;
		bool is_autosome() const ;
		bool is_missing() const ;
		
		operator std::string () const ;

	private:
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
