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

			std::vector< Entry > parse( std::string const& ) const ;
			std::vector< Entry > missing() const ;

			// In general, some VCF entries expect a number of values
			// that depends on the number of alleles and/or possible genotypes.
			// The number of possible genotypes can depend both on the number
			// of possible alleles and, e.g. for sex chromosomes, on other information
			// such as the sex of the individual.  To account for this we add two hooks
			// here which must be implemented in subclasses.
			virtual void set_number( std::size_t count ) = 0 ;
			virtual std::size_t get_number() const = 0;
		private:
			static std::string m_missing_value ;
			Spec const m_spec ;
			SimpleType::SharedPtr m_type ;
		} ;
		
		struct FixedNumberVCFEntryType: public VCFEntryType {
			FixedNumberVCFEntryType( std::size_t number, SimpleType::UniquePtr type ) ;
			void set_number( std::size_t ) { assert(0) ; }
			std::size_t get_number() const { return m_number ; }
		private:
			std::size_t m_number ;
		} ;

		struct DynamicNumberVCFEntryType: public VCFEntryType {
			DynamicNumberVCFEntryType( SimpleType::UniquePtr type ): VCFEntryType( type ) {} ;
			void set_number( std::size_t number ) { m_number = number ; }
			std::size_t get_number() const { return m_number ; }
		private:
			std::size_t m_number ;
		} ;

		struct OnePerAlleleVCFEntryType: public DynamicNumberVCFEntryType {
			OnePerAlleleVCFEntryType( SimpleType::UniquePtr type ): DynamicNumberVCFEntryType( type ) {} ;
		} ;

		struct OnePerGenotypeVCFEntryType: public DynamicNumberVCFEntryType {
			OnePerGenotypeVCFEntryType( SimpleType::UniquePtr type ): DynamicNumberVCFEntryType( type ) {} ;
		} ;

		std::auto_ptr< boost::ptr_map< std::string, VCFEntryType > > get_entry_types(
			std::multimap< std::string, VCFEntryType::Spec > const& spec,
			std::string const& key
		) ;
	}
}

#endif
