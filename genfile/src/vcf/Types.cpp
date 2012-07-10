
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <map>
#include <string>
#include <iostream>
#include <boost/ptr_container/ptr_map.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/string_utils/strtod.hpp"

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
		
		void EntrySetter::operator()( MissingValue const value ) { assert(0) ; }
		void EntrySetter::operator()( std::string& value ) { assert(0) ; }
		void EntrySetter::operator()( Integer const value ) { assert(0) ; }
		void EntrySetter::operator()( double const value ) { assert(0) ; }
		
		namespace impl {
			struct VariantEntryEntrySetter: public EntrySetter {
				VariantEntryEntrySetter( Entry& entry ): m_entry( entry ) {} 
				virtual void operator()( MissingValue const value ) { m_entry = value ; }
				virtual void operator()( std::string& value ) { m_entry = value ; }
				virtual void operator()( Integer const value ) { m_entry = value ; }
				virtual void operator()( double const value ) { m_entry = value ; }
			private:
				Entry& m_entry ;
			} ;
		}

		Entry SimpleType::parse( string_utils::slice const& value ) const {
			Entry result ;
			impl::VariantEntryEntrySetter setter( result ) ;
			parse( value, setter ) ;
			return result ;
		}

		void StringType::parse( string_utils::slice const& value, EntrySetter& setter ) const {
			std::string string_value = value ;
			setter( string_value ) ;
		}
		
		void IntegerType::parse( string_utils::slice const& value, EntrySetter& setter ) const {
			try {
				setter( string_utils::to_repr< EntrySetter::Integer >( value ) ) ;
			}
			catch( string_utils::StringConversionError const& ) {
				throw BadArgumentError( "genfile::vcf::IntegerType::parse()", "value = \"" + std::string( value ) + "\"" ) ;
			}
		}

		void FloatType::parse( string_utils::slice const& value, EntrySetter& setter ) const {
			try {
				setter( string_utils::strtod( value )) ;
				// return Entry( string_utils::to_repr< double >( value )) ;
			}
			catch( string_utils::StringConversionError const& ) {
				throw BadArgumentError( "genfile::vcf::FloatType::parse()", "value = \"" + std::string( value ) + "\"" ) ;
			}
		}

		void CharacterType::parse( string_utils::slice const& value, EntrySetter& setter ) const {
			if( value.size() != 1 ) {
				throw BadArgumentError("genfile::vcf::CharacterType::parse()", "value = \"" + std::string( value ) +"\"" ) ;
			}
			std::string string_value = value ;
			setter( string_value ) ;
		}

		void FlagType::parse( string_utils::slice const& value, EntrySetter& setter ) const {
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

			if( ID->second == "GT" || type->second == "Genotype" ) {
				if( number->second != "." && number->second != "1" ) {
					throw BadArgumentWithMessageError(
						"\"Number\" of GT field must be encoded as . or 1 in metadata.",
						"genfile::vcf::VCFEntryType::create()",
						"spec"
					) ;
				}
				if( type->second != "String" && type->second != "Genotype" ) {
					throw BadArgumentError( "genfile::vcf::VCFEntryType::create()", "spec" ) ;
				}
				result.reset( new GenotypeCallVCFEntryType()) ;
			}
			else if( number->second == "A" ) {
				result.reset( new OnePerAlternateAlleleVCFEntryType( SimpleType::create( type->second ))) ;
			}
			else if( number->second == "G" ) {
				result.reset( new OnePerGenotypeVCFEntryType( SimpleType::create( type->second ))) ;
			}
			else if( number->second == "." ) {
				result.reset(
					new DynamicNumberVCFEntryType(
						SimpleType::create( type->second )
					)
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
		
		namespace impl {
			struct VariantEntriesEntriesSetter: public EntriesSetter {
				VariantEntriesEntriesSetter( std::vector< Entry >& entries ): m_entries( entries ), entry_i( 0 ) {} 
				void set_number_of_entries( std::size_t n ) { m_entries.resize( n ) ; }
				virtual void operator()( MissingValue const value ) { m_entries[ entry_i++ ] = value ; }
				virtual void operator()( std::string& value ) { m_entries[ entry_i++ ] = value ; }
				virtual void operator()( Integer const value ) { m_entries[ entry_i++ ] = value ; }
				virtual void operator()( double const value ) { m_entries[ entry_i++ ] = value ; }
			private:
				std::vector< Entry >& m_entries ;
				std::size_t entry_i ;
			} ;
		}

		std::vector< Entry > VCFEntryType::parse( string_utils::slice const& value, std::size_t number_of_alleles, std::size_t ploidy ) const {
			std::vector< Entry > result ;
			impl::VariantEntriesEntriesSetter setter( result ) ;
			parse( value, number_of_alleles, ploidy, setter ) ;
			return result ;
		}

		std::vector< Entry > VCFEntryType::parse( string_utils::slice const& value, std::size_t number_of_alleles ) const {
			std::vector< Entry > result ;
			impl::VariantEntriesEntriesSetter setter( result ) ;
			parse( value, number_of_alleles, setter ) ;
			return result ;
		}

		void VCFEntryType::parse(
			string_utils::slice const& value,
			std::size_t number_of_alleles,
			std::size_t ploidy,
			EntriesSetter& setter
		) const {
			parse_elts( lex( value, number_of_alleles, ploidy ), setter ) ;
		}

		void VCFEntryType::parse( string_utils::slice const& value, std::size_t number_of_alleles, EntriesSetter& setter ) const {
			parse_elts( lex( value, number_of_alleles ), setter ) ;
		}
		
		void VCFEntryType::parse_elts( std::vector< string_utils::slice > const& elts, EntriesSetter& setter ) const {
			setter.set_number_of_entries( elts.size() ) ;
			for( std::size_t i = 0; i < elts.size(); ++i ) {
				if( elts[i] == m_missing_value ) {
					setter( MissingValue() ) ;
				}
				else {
					m_type->parse( elts[i], setter ) ;
				}
			}
		}
		
		std::vector< string_utils::slice > ListVCFEntryType::lex( string_utils::slice const& value, std::size_t number_of_alleles, std::size_t ploidy ) const {
			std::vector< string_utils::slice > result ;
			// An empty value is treated as an empty list, (not a list with one empty value)
			if( !value.empty() ) {
				result = string_utils::slice( value ).split( "," ) ;
			}
			ValueCountRange range = get_value_count_range( number_of_alleles, ploidy ) ;
			if( result.size() < range.first || result.size() > range.second ) {
				throw BadArgumentError( "genfile::vcf::ListVCFEntryType::lex()", "value = \"" + std::string( value ) + "\"" ) ;
			}
			return result ;
		}

		std::vector< string_utils::slice > ListVCFEntryType::lex( string_utils::slice const& value, std::size_t number_of_alleles ) const {
			std::vector< string_utils::slice > result ;
			// An empty value is treated as an empty list, (not a list with one empty value)
			if( !value.empty() ) {
				result = string_utils::slice( value ).split( "," ) ;
			}
			ValueCountRange range = get_value_count_range( number_of_alleles ) ;
			if( result.size() < range.first || result.size() > range.second ) {
				throw BadArgumentError( "genfile::vcf::ListVCFEntryType::lex()", "value = \"" + std::string( value ) + "\"" ) ;
			}
			return result ;
		}
		
		ListVCFEntryType::ValueCountRange ListVCFEntryType::get_value_count_range( std::size_t number_of_alleles, std::size_t ploidy ) const {
			return get_value_count_range( number_of_alleles ) ;
		}
		
		void ListVCFEntryType::get_missing_value( std::size_t number_of_alleles, std::size_t ploidy, EntriesSetter& setter ) const {
			std::size_t const N = get_value_count_range( number_of_alleles, ploidy ).first ;
			setter.set_number_of_entries( N ) ;
			for( std::size_t i = 0; i < N; ++i ) {
				setter( MissingValue() ) ;
			}
		}

		void ListVCFEntryType::get_missing_value( std::size_t number_of_alleles, EntriesSetter& setter ) const {
			std::size_t const N = get_value_count_range( number_of_alleles ).first ;
			setter.set_number_of_entries( N ) ;
			for( std::size_t i = 0; i < get_value_count_range( number_of_alleles ).first; ++i ) {
				setter( MissingValue() ) ;
			}
		}
		
		FixedNumberVCFEntryType::FixedNumberVCFEntryType( std::size_t number, SimpleType::UniquePtr type ):
			ListVCFEntryType( type ),
			m_number( number )
		{}
		
		
		namespace impl {
			std::size_t n_choose_k( std::size_t const n, std::size_t const k ) {
				// calculate n choose k, assuming no overflow, using the
				// multiplicative formula given on http://en.wikipedia.org/wiki/Binomial_coefficient
				double result = 1.0 ;
				for( std::size_t i = 1; i <= k; ++i ) {
					result *= double( n - k + i ) / double( i ) ;
				}
				return result ;
			}
			
			std::size_t get_number_of_unphased_genotypes( std::size_t const n_alleles, std::size_t const ploidy ) {
				// The number is equal to the number of ways to fill an n-vector
				// (where n is the number of alleles)
				// with nonnegative integers so that the sum is the ploidy.
				// This has a recursive expression involving filling the first
				// entry and then filling the others.
				// There is one special case: if the ploidy is zero, there are no genotypes at all.
				assert( n_alleles > 0 ) ;
				if( ploidy == 0 || n_alleles == 1 ) {
					return 1 ;
				}
				std::size_t result = 0 ;
				for( std::size_t i = 0; i <= ploidy; ++i ) {
					result += get_number_of_unphased_genotypes( n_alleles - 1, ploidy - i ) ;
				}
				return result ;
			}

			std::size_t get_number_of_phased_genotypes( std::size_t n_alleles, std::size_t ploidy ) {
				// The number is (number of alleles)^(ploidy)
				// except that if the ploidy is zero we report 0 phased genotypes, consistent
				// with 0 unphased genotypes.
				if( ploidy == 0 ) {
					return 0 ;
				}
				std::size_t result = 1 ;
				for( std::size_t i = 0; i < ploidy; ++i ) {
					result *= n_alleles ;
				}
				return result ;
			}
		}
		
		ListVCFEntryType::ValueCountRange OnePerGenotypeVCFEntryType::get_value_count_range( std::size_t number_of_alleles, std::size_t ploidy ) const {
			if( number_of_alleles == 0 && ploidy > 0 ) {
				throw BadArgumentError( "genfile::vcf::OnePerGenotypeVCFEntryType::get_value_count_range()", "number_of_alleles = 0" ) ;
			}
			std::size_t N = ( ploidy == 0 ) ? 0 : impl::get_number_of_unphased_genotypes( number_of_alleles, ploidy ) ;
			//std::cerr << "number_of_alleles = " << number_of_alleles << ", ploidy = " << ploidy << ", N = " << N << ".\n" ;
			return ValueCountRange( N, N ) ;
		}
		
		ListVCFEntryType::ValueCountRange OnePerGenotypeVCFEntryType::get_value_count_range( std::size_t number_of_alleles ) const {
			// OnePerGenotypeVCFEntryType requires the ploidy.
			assert(0) ;
		}
		
		GenotypeCallVCFEntryType::GenotypeCallVCFEntryType():
			VCFEntryType( SimpleType::create( "Integer" ))
		{}
		
		std::vector< string_utils::slice > GenotypeCallVCFEntryType::lex( string_utils::slice const& value, std::size_t number_of_alleles, std::size_t ploidy ) const {
			std::vector< string_utils::slice > elts = lex( value, number_of_alleles ) ;
			if( elts.size() != ploidy ) {
				throw BadArgumentError( "genfile::vcf::GenotypeCallVCFEntryType::lex()", "value = \"" + std::string( value ) + "\"" ) ;
			}
			return elts ;
		}

		std::vector< string_utils::slice > GenotypeCallVCFEntryType::lex( string_utils::slice const& value, std::size_t ) const {
			std::vector< string_utils::slice > elts ;
			// empty value is treated as empty list.
			if( !value.empty() ) {
				elts = string_utils::slice( value ).split( "|/" ) ;
			}
			return elts ;
		}

		void GenotypeCallVCFEntryType::get_missing_value( std::size_t, std::size_t ploidy, EntriesSetter& setter ) const {
			setter.set_number_of_entries( ploidy ) ;
			for( std::size_t i = 0; i < ploidy; ++i ) {
				setter( MissingValue() ) ;
			}
		}

		void GenotypeCallVCFEntryType::get_missing_value( std::size_t number_of_alleles, EntriesSetter& setter ) const {
			assert(0) ;
		}

		namespace impl {
			//
			// 
			struct CheckedGenotypeSetter: public EntriesSetter {
				CheckedGenotypeSetter( EntriesSetter& setter, EntriesSetter::Integer const max_genotype ):
					m_setter( setter ),
					m_max_genotype( max_genotype )
				{}
				
				void set_number_of_entries( std::size_t n ) { m_setter.set_number_of_entries( n ) ; }
				void operator()( MissingValue const value ) { m_setter( value ) ; }
				void operator()( Integer const value ) {
					if( value < 0 || value > m_max_genotype ) {
						throw BadArgumentError( "genfile::vcf::impl::GenotypeSetter::operator()", "value = \"" + string_utils::to_string( value ) + "\"" ) ;
					}
					m_setter( value ) ;
				}
			private:
				EntriesSetter& m_setter ;
				EntriesSetter::Integer m_max_genotype ;
			} ;
		}

		void GenotypeCallVCFEntryType::parse(
			string_utils::slice const& value,
			std::size_t number_of_alleles,
			EntriesSetter& setter
		) const {
			if( number_of_alleles == 0 ) {
				if( value != "" ) {
					throw BadArgumentError( "genfile::vcf::GenotypeCallVCFEntryType::parse()", "value = \"" + std::string( value ) + "\"" ) ;
				}
			}
			int const max = number_of_alleles - 1 ;

			std::vector< Entry > result ;
			
			// Most GT values have one character per allele.  We treat this as a special case.
			bool simple_parse_success = true ;

			if( value.size() == 3 && ( value[1] == '|' || value[1] == '/' ) ) {
				char a, b ;
				a = value[0] - '0' ;
				b = value[2] - '0' ;
				if((( value[0] == m_missing_value[0] ) || ( a >= 0 && a <= max ) ) && ( ( value[2] == m_missing_value[0] ) || ( b >= 0 && b <= max ) ) ) {
					setter.set_number_of_entries( 2 ) ;
					setter.set_order_type( ( value[1] == '|' ) ? EntriesSetter::eOrderedList : EntriesSetter::eUnorderedList ) ;

					if( value[0] == m_missing_value[0] ) {
						setter( MissingValue() ) ;
					} else {
						setter( EntriesSetter::Integer( a ) ) ;
					}
					if( value[2] == m_missing_value[0] ) {
						setter( MissingValue() ) ;
					} else {
						setter( EntriesSetter::Integer( b ) ) ;
					}
				}
				else {
					simple_parse_success = false ;
				}
			}
			else if( value.size() % 2 == 1 ) {
				std::vector< char > simple_values( ( value.size() + 1 ) / 2 ) ;
				for( std::size_t i = 0; i < value.size(); i += 2 ) {
					if( i > 0 && value[i-1] != '|' && value[i-1] != '/' ) {
						simple_parse_success = false ;
						break ;
					}
					if( value[i] == m_missing_value[0] ) {
						simple_values[ i / 2 ] = -1 ;
					}
					else if( value[i] >= '0' && value[i] <= ( '0' + max ) ) {
						simple_values[ i / 2 ] = ( value[i] - '0' ) ;
					}
					else {
						simple_parse_success = false ;
						break ;
					}
				}
				if( simple_parse_success ) {
					setter.set_number_of_entries( simple_values.size() ) ;
					setter.set_order_type( ( value[1] == '|' ) ? EntriesSetter::eOrderedList : EntriesSetter::eUnorderedList ) ;
					for( std::size_t i = 0; i < simple_values.size(); ++i ) {
						if( simple_values[i] == -1 ) {
							setter( MissingValue() ) ;
						} else {
							setter( EntriesSetter::Integer( simple_values[i] )) ;
						}
					}
				}
			}
			else {
				simple_parse_success = false ;
			}

			if( !simple_parse_success ) {
				impl::CheckedGenotypeSetter genotype_setter( setter, number_of_alleles - 1 ) ;
				std::vector< string_utils::slice > elts ;
				if( value.find( '/' ) == std::string::npos ) {
					// If there is a |, or no separator at all (ploidy = 1) the genotypes are phased.
					setter.set_order_type( EntriesSetter::eOrderedList ) ;
				} else {
					// / found, so genotypes are unphased.
					setter.set_order_type( EntriesSetter::eUnorderedList ) ;
				}
				parse_elts( lex( value, number_of_alleles ), genotype_setter ) ;
			}
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
				else if( !result.insert( where->second, VCFEntryType::create( range.first->second ) ).second ) {
					throw DuplicateKeyError( "result of genfile::vcf::get_entry_types", key + " field with ID=\"" + where->second + "\"" ) ;
				}
			}
			return result.release() ;
		}
	}
}
