
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_VCF_TYPES_HPP
#define GENFILE_VCF_TYPES_HPP

#include <cassert>
#include <limits>
#include <map>
#include <vector>
#include <string>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/types.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/string_utils/slice.hpp"


namespace genfile {
	namespace vcf {
		typedef genfile::VariantEntry Entry ;

		struct EntrySetter {
			typedef int64_t Integer ;
			typedef MissingValue MissingValue ;
			virtual ~EntrySetter() throw() {}
			virtual void set_value( std::size_t, MissingValue const value ) = 0 ;
			virtual void set_value( std::size_t, std::string& value ) ;
			virtual void set_value( std::size_t, Integer const value ) ;
			virtual void set_value( std::size_t, double const value ) ;
			void set_value( std::size_t i, VariantEntry const& value ) {
				if( value.is_missing() ) {
					this->set_value( i, genfile::MissingValue() ) ;
				} else if( value.is_string() ) {
					this->set_value( i, value.as< std::string >() ) ;
				} else if( value.is_int() ) {
					this->set_value( i, value.as< VariantEntry::Integer >() ) ;
				} else if( value.is_double() ) {
					this->set_value( i, value.as< double >() ) ;
				} else {
					assert(0) ;
				}
			}
		} ;

		struct EntriesSetter: public EntrySetter {
			typedef genfile::OrderType OrderType ;
			typedef genfile::ValueType ValueType ;
			virtual ~EntriesSetter() throw() {}
			virtual void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) = 0 ;
		} ;

		struct PerSampleEntriesSetter: public vcf::EntriesSetter, public boost::noncopyable {
			typedef std::auto_ptr< PerSampleEntriesSetter > UniquePtr ;
			virtual ~PerSampleEntriesSetter() throw() {}
			// Prepare to receive values for a variant.
			virtual void initialise( std::size_t nSamples, std::size_t nAlleles ) = 0 ;
			// Return true if we want values for this variant, otherwise false.
			virtual bool set_sample( std::size_t i ) = 0 ;
			// Finalise.
			virtual void finalise() = 0 ;
		} ;
		
		struct SimpleType: public boost::noncopyable {
			SimpleType() {}
			virtual ~SimpleType() {}
			typedef std::auto_ptr< SimpleType > UniquePtr ;
			typedef boost::shared_ptr< SimpleType > SharedPtr ;
			static UniquePtr create( std::string const& spec, std::string const& scale = "identity" ) ;
			virtual void parse( string_utils::slice const& value, EntrySetter& setter ) const = 0 ;
			Entry parse( std::size_t index, string_utils::slice const& value ) const ;
			virtual std::string to_string() const = 0 ;
			virtual ValueType represented_type() const = 0 ;
		} ;
		
		struct StringType: public SimpleType {
			using SimpleType::parse ;
			void parse( std::size_t index, string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "String" ; }
			ValueType represented_type() const { return eUnknownValueType ; }
		} ;
		
		struct IntegerType: public SimpleType {
			IntegerType( ValueType value_type = eUnknownValueType ): m_value_type( value_type ) {}
			using SimpleType::parse ;
			void parse( std::size_t index, string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "Integer" ; }
			ValueType represented_type() const { return m_value_type ; }
		private:
			ValueType const m_value_type ;
		} ;

		struct FloatType: public SimpleType {
			using SimpleType::parse ;
			void parse( std::size_t index, string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "Float" ; }
			ValueType represented_type() const { return eUnknownValueType ; }
		} ;

		struct ProbabilityType: public FloatType {
			std::string to_string() const { return "Probability" ; }
			ValueType represented_type() const { return eProbability ; } ;
		} ;

		struct PhredScaleFloatType: public SimpleType {
			using SimpleType::parse ;
			void parse( std::size_t index, string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "PhredScaleFloat" ; }
			ValueType represented_type() const { return eProbability ; }
		} ;

		struct CharacterType: public SimpleType {
			using SimpleType::parse ;
			void parse( std::size_t index, string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "Character" ; }
			ValueType represented_type() const { return eUnknownValueType ; }
		} ;

		struct FlagType: public SimpleType {
			using SimpleType::parse ;
			void parse( std::size_t index, string_utils::slice const& value, EntrySetter& setter ) const ;
			std::string to_string() const { return "Flag" ; }
			ValueType represented_type() const { return eFlag ; }
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
			// Parse a value returning values into the setter object.
			virtual void parse( string_utils::slice const&, std::size_t number_of_alleles, uint32_t ploidy, EntriesSetter& setter ) const = 0 ;
			// Return a set of missing values into the setter object.
			virtual void get_missing_value( std::size_t number_of_alleles, uint32_t ploidy, EntriesSetter& setter ) const = 0 ;
			// Return true if the entry type requires a valid (non-missing) ploidy value
			// in order to parse.
			virtual bool check_if_requires_ploidy() const = 0 ;
			// Return the underlying value type
			SimpleType const& get_value_type() const { return *m_type ; }
			
			std::string const& missing_value() const { return m_missing_value ; }
		protected:
			static std::string m_missing_value ;
		private:
			SimpleType::SharedPtr m_type ;
		} ;
		
		struct ListVCFEntryType: public VCFEntryType {
			ListVCFEntryType( SimpleType::UniquePtr type ): VCFEntryType( type ) {}
			void parse( string_utils::slice const&, std::size_t number_of_alleles, uint32_t ploidy, EntriesSetter& setter ) const ;
			void get_missing_value( std::size_t number_of_alleles, uint32_t ploidy, EntriesSetter& setter ) const ;
			// Return the range of valid possible value counts.
			typedef std::pair< std::size_t, std::size_t > ValueCountRange ;
			virtual ValueCountRange get_value_count_range( std::size_t number_of_alleles, uint32_t ploidy ) const = 0 ;
		protected:
			// A list entry type has a specific order type, statically known.
			virtual OrderType const get_order_type() const = 0 ;
		private:
			std::vector< string_utils::slice > m_storage ;
		} ;
		
		struct FixedNumberVCFEntryType: public ListVCFEntryType {
			FixedNumberVCFEntryType( std::size_t number, SimpleType::UniquePtr type ) ;
			ValueCountRange get_value_count_range( std::size_t, uint32_t ) const { return ValueCountRange( m_number, m_number ) ; }
			bool check_if_requires_ploidy() const { return false ; }
			OrderType const get_order_type() const { return eOrderedList ; }
		private:
			std::size_t m_number ;
		} ;

		struct DynamicNumberVCFEntryType: public ListVCFEntryType {
			DynamicNumberVCFEntryType( SimpleType::UniquePtr type ): ListVCFEntryType( type ) {} ;
			ValueCountRange get_value_count_range( std::size_t, uint32_t ) const {
				return ValueCountRange( 0, std::numeric_limits< std::size_t >::max() ) ;
			}
			bool check_if_requires_ploidy() const { return false ; }
			OrderType const get_order_type() const { return eOrderedList ; }
		} ;
		
		struct OnePerAlternateAlleleVCFEntryType: public ListVCFEntryType {
			OnePerAlternateAlleleVCFEntryType( SimpleType::UniquePtr type ): ListVCFEntryType( type ) {} ;
			ValueCountRange get_value_count_range( std::size_t number_of_alleles, uint32_t /* ploidy */ ) const {
				assert( number_of_alleles > 0 ) ;
				return ValueCountRange( number_of_alleles - 1, number_of_alleles - 1 ) ;
			}
			bool check_if_requires_ploidy() const { return false ; }
			OrderType const get_order_type() const { return eOrderedList ; }
		} ;
		
		struct OnePerAlleleVCFEntryType: public ListVCFEntryType {
			OnePerAlleleVCFEntryType( SimpleType::UniquePtr type ): ListVCFEntryType( type ) {} ;
			ValueCountRange get_value_count_range( std::size_t number_of_alleles, uint32_t /* ploidy */ ) const {
				return ValueCountRange( number_of_alleles, number_of_alleles ) ;
			}
			bool check_if_requires_ploidy() const { return false ; }
			OrderType const get_order_type() const { return ePerAllele ; }
		} ;
		
		struct OnePerGenotypeVCFEntryType: public ListVCFEntryType {
			OnePerGenotypeVCFEntryType( SimpleType::UniquePtr type ): ListVCFEntryType( type ) {} ;
			ValueCountRange get_value_count_range( std::size_t number_of_alleles, uint32_t ploidy ) const ;
			bool check_if_requires_ploidy() const { return true ; }
			OrderType const get_order_type() const { return ePerUnorderedGenotype ; }
		} ;

		struct GenotypeCallVCFEntryType: public VCFEntryType {
			GenotypeCallVCFEntryType() ;

			// Parse a value returning values into the setter object.
			void parse( string_utils::slice const&, std::size_t number_of_alleles, uint32_t ploidy, EntriesSetter& setter ) const ;
			// Return a set of missing values into the setter object.
			void get_missing_value( std::size_t number_of_alleles, uint32_t ploidy, EntriesSetter& setter ) const ;
			// Return true if the entry type requires a valid (non-missing) ploidy value
			// in order to parse.
			bool check_if_requires_ploidy() const { return false ; }
//			ValueCountRange get_value_count_range( std::size_t number_of_alleles, uint32_t ploidy ) const ;

		private:
			std::string const m_missing_value ;
			void lex( string_utils::slice const& value, std::vector< string_utils::slice >* result ) const ;
			// Parse a value returning values into the setter object.
			void parse_impl( string_utils::slice const&, std::size_t number_of_alleles, EntriesSetter& setter ) const ;
		} ;

		std::auto_ptr< boost::ptr_map< std::string, VCFEntryType > > get_entry_types(
			std::multimap< std::string, VCFEntryType::Spec > const& spec,
			std::string const& key
		) ;
	}
}

#endif
