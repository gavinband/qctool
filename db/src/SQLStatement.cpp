
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <string>
#include <stdint.h>
#include "sqlite3/sqlite3.h"
#include "db/SQLStatement.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/string_utils.hpp"

namespace db {
	SQLStatement::~SQLStatement() {}
	
	template<>
	int SQLStatement::get_column< int >( int column_id ) const {
		return this->get_column_int( column_id ) ;
	}

	template<>
	int64_t SQLStatement::get_column< int64_t >( int column_id ) const {
		return this->get_column_int64( column_id ) ;
	}

	template<>
	double SQLStatement::get_column< double >( int column_id ) const {
		return this->get_column_double( column_id ) ;
	}

	template<>
	std::string SQLStatement::get_column< std::string >( int column_id ) const {
		return this->get_column_string( column_id ) ;
	}

	template<>
	char SQLStatement::get_column< char >( int column_id ) const {
		return this->get_column_char( column_id ) ;
	}
	
	SQLStatement& SQLStatement::bind( std::size_t i, char const* value ) {
		bind( i, std::string( value )) ;
		return *this ;
	}
	
	namespace impl {
		struct variant_binder: public boost::static_visitor<>
		{
			variant_binder( std::size_t index, SQLStatement& statement ): m_index( index ), m_statement( statement ) {}
			template< typename T >
			void operator()( T const& value ) {
				m_statement.bind( m_index, value ) ;
			}
			void operator()( genfile::Chromosome const& value ) {
				m_statement.bind( m_index, std::string( value ) ) ;
			}
			void operator()( genfile::MissingValue const& value ) {
				m_statement.bind_NULL( m_index ) ;
			}
			void operator()( genfile::GenomePosition const& value ) {
				m_statement.bind( m_index, genfile::string_utils::to_string( value ) ) ;
			}
		private:
			std::size_t const m_index ;
			SQLStatement& m_statement ;
		} ;
	}

	SQLStatement& SQLStatement::bind( std::size_t index, genfile::VariantEntry const& value ) {
		genfile::apply_visitor( impl::variant_binder( index, *this ), value ) ;
		return *this ;
	}
}
