#include <map>
#include <string>
#include <boost/ptr_container/ptr_map.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/string_utils/slice.hpp"

namespace genfile {
	namespace vcf {
			SimpleType::UniquePtr SimpleType::create( std::string const& spec ) {
				SimpleType::UniquePtr result ;
				if( spec == "String" ) {
					result.reset( new StringType() ) ;
				}
				else if( spec == "Integer" ) {
					result.reset( new IntegerType() ) ;
				}
				else if( spec == "Float" ) {
					result.reset( new FloatType() ) ;
				}
				else if( spec == "Character" ) {
					result.reset( new CharacterType() ) ;
				}
				else if( spec == "Flag" ) {
					result.reset( new FlagType() ) ;
				}
				else {
					throw BadArgumentError( "genfile::vcf::types::SimpleType::create()", "spec = \"" + spec + "\"" ) ;
				}
				return result ;
			}
			
			Entry StringType::parse( string_utils::slice const& value ) const {
				return Entry( value )  ;
			}
			
			Entry IntegerType::parse( string_utils::slice const& value ) const {
				try {
					return Entry( string_utils::to_repr< int >( value ) ) ;
				}
				catch( string_utils::StringConversionError const& ) {
					throw BadArgumentError( "genfile::vcf::IntegerType::parse()", "value = \"" + std::string( value ) + "\"" ) ;
				}
			}

			Entry FloatType::parse( string_utils::slice const& value ) const {
				try {
					return Entry( string_utils::to_repr< double >( value )) ;
				}
				catch( string_utils::StringConversionError const& ) {
					throw BadArgumentError( "genfile::vcf::FloatType::parse()", "value = \"" + std::string( value ) + "\"" ) ;
				}
			}

			Entry CharacterType::parse( string_utils::slice const& value ) const {
				if( value.size() != 1 ) {
					throw BadArgumentError("genfile::vcf::CharacterType::parse()", "value = \"" + std::string( value ) +"\"" ) ;
				}
				return Entry( value ) ;
			}

			Entry FlagType::parse( string_utils::slice const& value ) const {
				assert(0) ;
			}
			
			std::string VCFEntryType::m_missing_value = "." ;
			
			VCFEntryType::UniquePtr VCFEntryType::create( Spec const& spec ) {
				Spec::const_iterator ID = spec.find( "ID" ) ;
				Spec::const_iterator number = spec.find( "Number" ) ;
				Spec::const_iterator type = spec.find( "Type" ) ;
				if( ID == spec.end() || number == spec.end() || type == spec.end() ) {
					throw BadArgumentError( "genfile::vcf::VCFEntryType::create()", "spec" ) ;
				}
				VCFEntryType::UniquePtr result ;

				if( ID->second == "GT" ) {
					if( number->second != "." || type->second != "Integer" ) {
						throw BadArgumentError( "genfile::vcf::VCFEntryType::create()", "spec" ) ;
					}
					result.reset( new GenotypeCallVCFEntryType( SimpleType::create( type->second ))) ;
				}
				else if( number->second == "A" ) {
					result.reset( new OnePerAlternateAlleleVCFEntryType( SimpleType::create( type->second ))) ;
				}
				else if( number->second == "G" ) {
					result.reset( new OnePerGenotypeVCFEntryType( SimpleType::create( type->second ))) ;
				}
				else if( number->second == "." ) {
					new DynamicNumberVCFEntryType(
						SimpleType::create( type->second )
					) ;
				}
				else {
					result.reset(
						new FixedNumberVCFEntryType(
							string_utils::to_repr< std::size_t >( number->second ),
							SimpleType::create( type->second )
						)
					) ;
				}
				return result ;
			}
			
			VCFEntryType::VCFEntryType( SimpleType::UniquePtr type ):
				m_type( type )
			{
			}
			
			std::vector< Entry > VCFEntryType::parse( std::string const& value, std::size_t number_of_alleles, std::size_t ploidy ) const {
				return parse_elts( lex( value, number_of_alleles, ploidy )) ;
			}

			std::vector< Entry > VCFEntryType::parse( std::string const& value, std::size_t number_of_alleles ) const {
				return parse_elts( lex( value, number_of_alleles )) ;
			}
			
			std::vector< Entry > VCFEntryType::parse_elts( std::vector< string_utils::slice > const& elts ) const {
				std::vector< Entry > result( elts.size() ) ;
				for( std::size_t i = 0; i < result.size(); ++i ) {
					if( elts[i] == m_missing_value ) {
						result[i] = Entry( MissingValue() ) ;
					}
					else {
						result[i] = m_type->parse( elts[i] ) ;
					}
				}
				return result ;
			}
			
			std::vector< string_utils::slice > ListVCFEntryType::lex( std::string const& value, std::size_t number_of_alleles, std::size_t ) const {
				return lex( value, number_of_alleles ) ;
			}

			std::vector< string_utils::slice > ListVCFEntryType::lex( std::string const& value, std::size_t number_of_alleles ) const {
				std::vector< string_utils::slice > result = string_utils::slice( value ).split( "," ) ;
				ValueCountRange range = get_value_count_range( number_of_alleles ) ;
				if( result.size() < range.first || result.size() > range.second ) {
					throw BadArgumentError( "genfile::vcf::ListVCFEntryType::lex()", "value = \"" + value + "\"" ) ;
				}
				return result ;
			}
			
			ListVCFEntryType::ValueCountRange ListVCFEntryType::get_value_count_range( std::size_t number_of_alleles, std::size_t ploidy ) const {
				return get_value_count_range( number_of_alleles ) ;
			}
			
			std::vector< Entry > ListVCFEntryType::get_missing_value( std::size_t number_of_alleles, std::size_t ploidy ) const {
				return std::vector< Entry >( get_value_count_range( number_of_alleles, ploidy ).first, MissingValue() ) ;
			}

			std::vector< Entry > ListVCFEntryType::get_missing_value( std::size_t number_of_alleles ) const {
				return std::vector< Entry >( get_value_count_range( number_of_alleles ).first, MissingValue() ) ;
			}
			
			FixedNumberVCFEntryType::FixedNumberVCFEntryType( std::size_t number, SimpleType::UniquePtr type ):
				ListVCFEntryType( type ),
				m_number( number )
			{}
			
			
			namespace impl {
				std::size_t n_choose_k( std::size_t const n, std::size_t const k ) {
					// calculate n choose k, assuming no overflow, using the
					// multiplicative formula given on http://en.wikipedia.org/wiki/Binomial_coefficient
					double result = 0 ;
					for( std::size_t i = 1; i <= k; ++i ) {
						result *= double( n - k - i ) / double( i ) ;
					}
					return std::size_t( result ) ;
				}
			}
			
			ListVCFEntryType::ValueCountRange OnePerGenotypeVCFEntryType::get_value_count_range( std::size_t number_of_alleles, std::size_t ploidy ) const {
				std::size_t N = impl::n_choose_k( number_of_alleles, ploidy ) ;
				return ValueCountRange( N, N ) ;
			}
			
			ListVCFEntryType::ValueCountRange OnePerGenotypeVCFEntryType::get_value_count_range( std::size_t number_of_alleles ) const {
				// OnePerGenotypeVCFEntryType requires the ploidy.
				assert(0) ;
			}
			
			std::vector< string_utils::slice > GenotypeCallVCFEntryType::lex( std::string const& value, std::size_t, std::size_t ploidy ) const {
				std::vector< string_utils::slice > elts = string_utils::slice( value ).split( "|/" ) ;
				if( elts.size() != ploidy ) {
					throw BadArgumentError( "genfile::vcf::GenotypeCallVCFEntryType::lex()", "value = \"" + value + "\"" ) ;
				}
				return elts ;
			}

			std::vector< string_utils::slice > GenotypeCallVCFEntryType::lex( std::string const& value, std::size_t ) const {
				return string_utils::slice( value ).split( "|/" ) ;
			}

			std::vector< Entry > GenotypeCallVCFEntryType::get_missing_value( std::size_t, std::size_t ploidy ) const {
				return std::vector< Entry >( ploidy, MissingValue() ) ;
			}

			std::vector< Entry > GenotypeCallVCFEntryType::get_missing_value( std::size_t number_of_alleles ) const {
				assert(0) ;
			}

			std::vector< Entry > GenotypeCallVCFEntryType::parse(
				std::string const& value,
				std::size_t number_of_alleles
			) const {
				return parse_elts( string_utils::slice( value ).split( "|/" ) ) ;
			}

			std::auto_ptr< boost::ptr_map< std::string, VCFEntryType > > get_entry_types(
				std::multimap< std::string, VCFEntryType::Spec > const& metadata,
				std::string const& key
			) {
				typedef std::multimap< std::string, VCFEntryType::Spec >::const_iterator MetadataIterator ;
				std::pair< MetadataIterator, MetadataIterator > range = metadata.equal_range( key ) ;
				boost::ptr_map< std::string, VCFEntryType > result ;
				for( ; range.first != range.second; ++range.first ) {
					VCFEntryType::Spec::const_iterator where = range.first->second.find( "ID" ) ;
					if( where == range.first->second.end() ) {
						throw BadArgumentError( "genfile::vcf::get_specs()", "metadata" ) ;
					}
					else {
						if( !result.insert( where->second, VCFEntryType::create( range.first->second ) ).second ) {
							assert(0) ;
						}
					}
				}
				return result.release() ;
			}
	}
}
