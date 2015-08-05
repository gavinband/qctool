
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_FILTERING_STATSOURCE_HPP
#define STATFILE_FILTERING_STATSOURCE_HPP

#include <boost/bimap.hpp>
#include "statfile/BuiltInTypeStatSource.hpp"
#include "statfile/BuiltInTypeStatSourceChain.hpp"
#include "genfile/VariantEntry.hpp"

namespace statfile {
	namespace impl {
		// 
		// Homogenise types for value comparison.
		// We do this:
		// * 1. convert any ints to doubles
		// * 2. of one operand is int or double and the other operand is string then we cast the string to a double if possible.
		//
		void cast_types_for_comparison( genfile::VariantEntry& value1, genfile::VariantEntry& value2 ) ;
	}
	
	struct Constraint {
		enum Type { eEqualTo, eNotEqualTo, eLessThan, eLessThanOrEqualTo, eGreaterThan, eGreaterThanOrEqualTo, eBetween } ;

		static Constraint equal_to( genfile::VariantEntry const& value ) ;
		static Constraint not_equal_to( genfile::VariantEntry const& value ) ;
		static Constraint less_than( genfile::VariantEntry const& value ) ;
		static Constraint less_than_or_equal_to( genfile::VariantEntry const& value ) ;
		static Constraint greater_than( genfile::VariantEntry const& value ) ;
		static Constraint greater_than_or_equal_to( genfile::VariantEntry const& value ) ;
		static Constraint between( genfile::VariantEntry const& value1, genfile::VariantEntry const& value2 ) ;
		Constraint( Constraint const& other ) ;
		Constraint& operator=( Constraint const& other ) ;
		bool test( genfile::VariantEntry value ) const ;
		
	private:
		Type m_type ;
		genfile::VariantEntry m_value1 ;
		genfile::VariantEntry m_value2 ;
		
	private:
		Constraint( genfile::VariantEntry const& value1, genfile::VariantEntry const& value2, Type type ) ;
	} ;
	
	struct BoundConstraint {
	public:
		static BoundConstraint parse( std::string const& spec ) ;
	public:
		BoundConstraint( std::string const& variable, Constraint const& constraint ) ;
		BoundConstraint( BoundConstraint const& other ) ;
		BoundConstraint& operator=( BoundConstraint const& other ) ;

	private:
		std::string const m_variable ;
		Constaint m_constraint ;
	} ;
	
	class FilteringStatSource: public BuiltInTypeStatSource {
	public:
		FilteringStatSource(
			BuiltInTypeStatSource::UniquePtr source,
			std::string const column,
			Constraint const& constraint
		) ;
		operator bool() const ;

	private:
		OptionalCount number_of_rows() const { return OptionalCount() ; }
		std::size_t number_of_columns() const { return m_source->number_of_columns() ; }
		std::vector< std::string > column_names() const { return m_source->column_names() ; }
		std::string get_descriptive_text() const { return m_source->get_descriptive_text() ; }
		using Base::read_value ;
		void read_value( int32_t& value ) ;
		void read_value( uint32_t& value ) ;
		void read_value( std::string& value ) ;
		void read_value( double& value ) ;
		void ignore_value() ;
		void ignore_all() ;
		void end_row() ;
		void restart_row() ;
		void move_to_next_row_impl() ;
		void move_to_next_matching_row() ;
		bool row_matches() ;
		
	private:
		BuiltInTypeStatSource::UniquePtr m_source ;
		typedef boost::bimap< std::string, std::size_t > ColumnMap ;
		ColumnMap m_columns ;
		Constraint m_constraint ;
	} ;
}

#endif
