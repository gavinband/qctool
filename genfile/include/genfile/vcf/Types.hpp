
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VCF_TYPES_HPP
#define GENFILE_VCF_TYPES_HPP

#include <limits>
#include <map>
#include <vector>
#include <string>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/string_utils/slice.hpp"

namespace genfile {
	namespace vcf {
		typedef genfile::VariantEntry Entry ;

		struct EntrySetter {
			typedef int64_t Integer ;
			virtual ~EntrySetter() throw() {}
			virtual void operator()( MissingValue const value ) ;
			virtual void operator()( std::string& value ) ;
			virtual void operator()( Integer const value ) ;
			virtual void operator()( double const value ) ;
		} ;

		struct EntriesSetter: public EntrySetter {
			virtual ~EntriesSetter() throw() {}
			virtual void set_number_of_entries( std::size_t n ) = 0 ;
			enum OrderType { eUnorderedList = 0, eOrderedList = 1 } ;
			virtual void set_order_type( OrderType const type ) {} ;
		} ;
		
		struct PerSampleEntriesSetter: public vcf::EntriesSetter, public boost::noncopyable {
			typedef std::auto_ptr< PerSampleEntriesSetter > UniquePtr ;
			virtual ~PerSampleEntriesSetter() throw() {}
			virtual void set_number_of_samples( std::size_t n ) = 0 ;
			virtual void set_sample( std::size_t i ) = 0 ;
		} ;
		
		struct SimpleType: public boost::noncopyable {
			SimpleType() {}
			virtual ~SimpleType() {}
			typedef std::auto_ptr< SimpleType > UniquePtr ;
			typedef boost::shared_ptr< SimpleType > SharedPtr ;
			static UniquePtr create( std::string const& spec, std::string const& scale = "identity" ) ;
			virtual void parse( string_utils::slice const& value, EntrySetter& setter ) const = 0 ;
			Entry parse( string_utils::slice const& value ) const ;
			virtual std::string to_string() const = 0 ;
		} ;
		
		struct StringType: public SimpleType {
			using SimpleType::parse ;
			void parse( string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "String" ; }
		} ;
		
		struct IntegerType: public SimpleType {
			using SimpleType::parse ;
			void parse( string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "Integer" ; }
		} ;

		struct FloatType: public SimpleType {
			using SimpleType::parse ;
			void parse( string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "Float" ; }
		} ;

		struct PhredScaleFloatType: public SimpleType {
			using SimpleType::parse ;
			void parse( string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "PhredScaleFloat" ; }
		} ;

		struct CharacterType: public SimpleType {
			using SimpleType::parse ;
			void parse( string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "Character" ; }
		} ;

		struct FlagType: public SimpleType {
			using SimpleType::parse ;
			void parse( string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "Flag" ; }
		} ;
		
		struct VCFEntryType: public boost::noncopyable
		{
		public:
			typedef std::auto_ptr< VCFEntryType > UniquePtr ;
			typedef std::map< std::string, std::string > Spec ;

			static UniquePtr create( Spec const& spec ) ;
		
		public:
			VCFEntryType( SimpleType::UniquePtr type ) ;
			virtual ~VCFEntryType() {}

			virtual void parse( string_utils::slice const&, std::size_t number_of_alleles, std::size_t ploidy, EntriesSetter& setter ) const ;
			virtual void parse( string_utils::slice const&, std::size_t number_of_alleles, EntriesSetter& setter ) const ;

			// Convenience functions for legacy interface.
			// This interface is deprecated because it involves lots of small memory allocations, which
			// slows things down too much.
			std::vector< Entry > parse( string_utils::slice const&, std::size_t number_of_alleles, std::size_t ploidy ) const ;
			std::vector< Entry > parse( string_utils::slice const&, std::size_t number_of_alleles ) const ;

			virtual void get_missing_value( std::size_t number_of_alleles, std::size_t ploidy, EntriesSetter& setter ) const = 0 ;
			virtual void get_missing_value( std::size_t number_of_alleles, EntriesSetter& setter ) const = 0 ;
			virtual bool check_if_requires_ploidy() const = 0 ;

			virtual SimpleType const& get_type() const { return *m_type ; }
		protected:
			virtual std::vector< string_utils::slice > lex(
				string_utils::slice const& value,
				std::size_t number_of_alleles,
				std::size_t ploidy
			) const = 0 ;

			virtual std::vector< string_utils::slice > lex(
				string_utils::slice const& value,
				std::size_t number_of_alleles
			) const = 0 ;

			void parse_elts( std::vector< string_utils::slice > const& elts, EntriesSetter& setter ) const ;
		protected:
			static std::string m_missing_value ;
		private:
			SimpleType::SharedPtr m_type ;
		} ;
		
		struct ListVCFEntryType: public VCFEntryType {
			ListVCFEntryType( SimpleType::UniquePtr type ): VCFEntryType( type ) {}
			std::vector< string_utils::slice > lex( string_utils::slice const& value, std::size_t number_of_alleles, std::size_t ploidy ) const ;
			std::vector< string_utils::slice > lex( string_utils::slice const& value, std::size_t number_of_alleles ) const ;
			virtual void get_missing_value( std::size_t number_of_alleles, std::size_t ploidy, EntriesSetter& setter ) const ;
			virtual void get_missing_value( std::size_t number_of_alleles, EntriesSetter& setter ) const ;
			typedef std::pair< std::size_t, std::size_t > ValueCountRange ;
			virtual ValueCountRange get_value_count_range( std::size_t number_of_alleles, std::size_t ploidy ) const ;
			virtual ValueCountRange get_value_count_range( std::size_t number_of_alleles ) const = 0 ;
			bool check_if_requires_ploidy() const { return false ; }
		} ;
		
		struct FixedNumberVCFEntryType: public ListVCFEntryType {
			FixedNumberVCFEntryType( std::size_t number, SimpleType::UniquePtr type ) ;
			ValueCountRange get_value_count_range( std::size_t ) const { return ValueCountRange( m_number, m_number ) ; }
		private:
			std::size_t m_number ;
		} ;

		struct DynamicNumberVCFEntryType: public ListVCFEntryType {
			DynamicNumberVCFEntryType( SimpleType::UniquePtr type ): ListVCFEntryType( type ) {} ;
			ValueCountRange get_value_count_range( std::size_t ) const {
				return ValueCountRange( 0, std::numeric_limits< std::size_t >::max() ) ;
			}
		} ;
		
		struct OnePerAlternateAlleleVCFEntryType: public ListVCFEntryType {
			OnePerAlternateAlleleVCFEntryType( SimpleType::UniquePtr type ): ListVCFEntryType( type ) {} ;
			ValueCountRange get_value_count_range( std::size_t number_of_alleles ) const {
				return ValueCountRange( number_of_alleles, number_of_alleles ) ;
			}
		} ;
		
		struct OnePerGenotypeVCFEntryType: public ListVCFEntryType {
			OnePerGenotypeVCFEntryType( SimpleType::UniquePtr type ): ListVCFEntryType( type ) {} ;
			ValueCountRange get_value_count_range( std::size_t number_of_alleles, std::size_t ploidy ) const ;
			ValueCountRange get_value_count_range( std::size_t number_of_alleles ) const ;
			bool check_if_requires_ploidy() const { return true ; }
		} ;
		
		struct GenotypeCallVCFEntryType: public VCFEntryType {
			GenotypeCallVCFEntryType() ;

			std::vector< string_utils::slice > lex( string_utils::slice const& value, std::size_t number_of_alleles, std::size_t ploidy ) const ;
			std::vector< string_utils::slice > lex( string_utils::slice const& value, std::size_t number_of_alleles ) const ;
			void get_missing_value( std::size_t, std::size_t ploidy, EntriesSetter& setter ) const ;
			void get_missing_value( std::size_t number_of_alleles, EntriesSetter& setter ) const ;

			bool check_if_requires_ploidy() const { return true; }

			// A special use of the genotype call is to infer ploidy for the other data.
			// For this use we need to specialise the parse() function.
			using VCFEntryType::parse ;
			void parse( string_utils::slice const& value, std::size_t number_of_alleles, EntriesSetter& setter ) const ;
		} ;

		std::auto_ptr< boost::ptr_map< std::string, VCFEntryType > > get_entry_types(
			std::multimap< std::string, VCFEntryType::Spec > const& spec,
			std::string const& key
		) ;
	}
}

#endif
