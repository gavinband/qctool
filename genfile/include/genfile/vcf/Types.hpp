#ifndef GENFILE_VCF_TYPES_HPP
#define GENFILE_VCF_TYPES_HPP

#include <map>
#include <vector>
#include <string>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/noncopyable.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/VariantEntry.hpp"

namespace genfile {
	namespace vcf {
		typedef genfile::VariantEntry Entry ;
		
		struct SimpleType: public boost::noncopyable {
			SimpleType() {}
			virtual ~SimpleType() {}
			typedef std::auto_ptr< SimpleType > UniquePtr ;
			typedef boost::shared_ptr< SimpleType > SharedPtr ;
			static UniquePtr create( std::string const& spec ) ;
			virtual Entry parse( std::string const& value ) const = 0 ;
		} ;
		
		struct StringType: public SimpleType {
			Entry parse( std::string const& value ) const ;
		} ;
		
		struct IntegerType: public SimpleType {
			Entry parse( std::string const& value ) const ;
		} ;

		struct FloatType: public SimpleType {
			Entry parse( std::string const& value ) const ;
		} ;

		struct CharacterType: public SimpleType {
			Entry parse( std::string const& value ) const ;
		} ;

		struct FlagType: public SimpleType {
			Entry parse( std::string const& value ) const ;
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

			std::vector< Entry > parse( std::string const&, std::size_t number_of_alleles, std::size_t ploidy ) const ;
			std::vector< Entry > parse( std::string const&, std::size_t number_of_alleles ) const ;

			virtual std::vector< Entry > get_missing_value( std::size_t number_of_alleles, std::size_t ploidy ) const = 0 ;
			virtual std::vector< Entry > get_missing_value( std::size_t number_of_alleles ) const = 0 ;

		protected:
			virtual std::vector< std::string > lex(
				std::string const& value,
				std::size_t number_of_alleles,
				std::size_t ploidy
			) const = 0 ;

			virtual std::vector< std::string > lex(
				std::string const& value,
				std::size_t number_of_alleles
			) const = 0 ;

			std::vector< Entry > parse_elts( std::vector< std::string > const& elts ) const ;

		private:
			static std::string m_missing_value ;
			SimpleType::SharedPtr m_type ;
		} ;
		
		struct ListVCFEntryType: public VCFEntryType {
			ListVCFEntryType( SimpleType::UniquePtr type ): VCFEntryType( type ) {}
			std::vector< std::string > lex( std::string const& value, std::size_t number_of_alleles, std::size_t ploidy ) const ;
			std::vector< std::string > lex( std::string const& value, std::size_t number_of_alleles ) const ;
			virtual std::vector< Entry > get_missing_value( std::size_t number_of_alleles, std::size_t ploidy ) const ;
			virtual std::vector< Entry > get_missing_value( std::size_t number_of_alleles ) const ;

			typedef std::pair< std::size_t, std::size_t > ValueCountRange ;
			virtual ValueCountRange get_value_count_range( std::size_t number_of_alleles, std::size_t ploidy ) const {
				return get_value_count_range( number_of_alleles ) ;
			} ;
			virtual ValueCountRange get_value_count_range( std::size_t number_of_alleles ) const = 0 ;
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
		} ;
		
		struct GenotypeCallVCFEntryType: public VCFEntryType {
			GenotypeCallVCFEntryType( SimpleType::UniquePtr type ): VCFEntryType( type ) {} ;

			std::vector< std::string > lex( std::string const& value, std::size_t number_of_alleles, std::size_t ploidy ) const ;
			std::vector< std::string > lex( std::string const& value, std::size_t number_of_alleles ) const ;
			std::vector< Entry > get_missing_value( std::size_t, std::size_t ploidy ) const ;
			std::vector< Entry > get_missing_value( std::size_t number_of_alleles ) const ;

			// A special use of the genotype call is to infer ploidy for the other data.
			// To support this we give them their own extra parse function which does not
			// take the ploidy as an argument.
			std::vector< Entry > parse( std::string const& value, std::size_t number_of_alleles ) const ;
		} ;

		std::auto_ptr< boost::ptr_map< std::string, VCFEntryType > > get_entry_types(
			std::multimap< std::string, VCFEntryType::Spec > const& spec,
			std::string const& key
		) ;
	}
}

#endif
