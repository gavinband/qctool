#include <map>
#include <string>
#include <boost/ptr_container/ptr_map.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/string_utils.hpp"

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
			
			Entry StringType::parse( std::string const& value ) const {
				return Entry( value )  ;
			}
			
			Entry IntegerType::parse( std::string const& value ) const {
				return Entry( string_utils::to_repr< int >( value ) ) ;
			}

			Entry FloatType::parse( std::string const& value ) const {
				return Entry( string_utils::to_repr< double >( value )) ;
			}

			Entry CharacterType::parse( std::string const& value ) const {
				if( !value.size() == 1 ) {
					throw BadArgumentError("genfile::vcf::types::CharacterType::parse()", "value = \"" + value +"\"" ) ;
				}
				return Entry( value ) ;
			}

			Entry FlagType::parse( std::string const& value ) const {
				assert(0) ;
			}
			
			std::string VCFEntryType::m_missing_value = "." ;
			
			VCFEntryType::UniquePtr VCFEntryType::create( Spec const& spec ) {
				Spec::const_iterator number = spec.find( "Number" ) ;
				Spec::const_iterator type = spec.find( "Type" ) ;
				if( number == spec.end() || type == spec.end() ) {
					throw BadArgumentError( "genfile::vcf::VCFEntryType::create()", "spec" ) ;
				}
				VCFEntryType::UniquePtr result ;
				if( number->second == "A" ) {
					result.reset( new OnePerAlleleVCFEntryType( SimpleType::create( type->second ))) ;
				}
				else if( number->second == "G" ) {
					result.reset( new OnePerGenotypeVCFEntryType( SimpleType::create( type->second ))) ;
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
			
			std::vector< Entry > VCFEntryType::parse( std::string const& value ) const {
				std::vector< std::string > elts = string_utils::split( value, "," ) ;
				if( elts.size() != std::size_t( get_number() ) ) {
					throw BadArgumentError( "genfile::vcf::VCFEntryType::parse()", "value = \"" + value + "\"" ) ;
				}
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
			
			std::vector< Entry > VCFEntryType::missing() const {
				return std::vector< Entry >( get_number(), MissingValue() ) ;
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
						result.insert( where->first, VCFEntryType::create( range.first->second ) ) ;
					}
				}
				return result.release() ;
			}
	}
}
