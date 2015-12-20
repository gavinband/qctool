
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <map>
#include <string>
#include <iostream>
#include <cmath>
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
		SimpleType::UniquePtr SimpleType::create( std::string const& spec, std::string const& scale ) {
			SimpleType::UniquePtr result ;
			if( spec == "String" ) {
				assert( scale == "identity" ) ;
				result.reset( new StringType() ) ;
			}
			else if( spec == "Integer" ) {
				assert( scale == "identity" ) ;
				result.reset( new IntegerType( eUnknownValueType ) ) ;
			}
			else if( spec == "AlleleIndex" ) {
				assert( scale == "identity" ) ;
				result.reset( new IntegerType( eAlleleIndex) ) ;
			}
			else if( spec == "Float" ) {
				if( scale == "identity" ) {
					result.reset( new FloatType() ) ;
				} else {
					result.reset( new PhredScaleFloatType() ) ;
				}
			}
			else if( spec == "Probability" ) {
				result.reset( new ProbabilityType() ) ;
			}
			else if( spec == "Character" ) {
				assert( scale == "identity" ) ;
				result.reset( new CharacterType() ) ;
			}
			else if( spec == "Flag" ) {
				assert( scale == "identity" ) ;
				result.reset( new FlagType() ) ;
			}
			else {
				throw BadArgumentError( "genfile::vcf::types::SimpleType::create()", "spec = \"" + spec + "\"" ) ;
			}
			return result ;
		}
		
		void EntrySetter::set_value( std::size_t, std::string& value ) { assert(0) ; }
		void EntrySetter::set_value( std::size_t, Integer const value ) { assert(0) ; }
		void EntrySetter::set_value( std::size_t, double const value ) { assert(0) ; }
		
		namespace impl {
			struct VariantEntryEntrySetter: public EntrySetter {
				VariantEntryEntrySetter( Entry& entry ): m_entry( entry ) {} 
				virtual void set_value( std::size_t, MissingValue const value ) { m_entry = value ; }
				virtual void set_value( std::size_t, std::string& value ) { m_entry = value ; }
				virtual void set_value( std::size_t, Integer const value ) { m_entry = value ; }
				virtual void set_value( std::size_t, double const value ) { m_entry = value ; }
			private:
				Entry& m_entry ;
			} ;
		}

		Entry SimpleType::parse( std::size_t i, string_utils::slice const& value ) const {
			Entry result ;
			impl::VariantEntryEntrySetter setter( result ) ;
			parse( i, value, setter ) ;
			return result ;
		}

		void StringType::parse( std::size_t i, string_utils::slice const& value, EntrySetter& setter ) const {
			std::string string_value = value ;
			setter.set_value( i, string_value ) ;
		}
		
		void IntegerType::parse( std::size_t i, string_utils::slice const& value, EntrySetter& setter ) const {
			try {
				setter.set_value( i, string_utils::to_repr< EntrySetter::Integer >( value ) ) ;
			}
			catch( string_utils::StringConversionError const& ) {
				throw BadArgumentError( "genfile::vcf::IntegerType::parse()", "value = \"" + std::string( value ) + "\"" ) ;
			}
		}

		void FloatType::parse( std::size_t i, string_utils::slice const& value, EntrySetter& setter ) const {
			try {
				setter.set_value( i, string_utils::strtod( value )) ;
			}
			catch( string_utils::StringConversionError const& ) {
				throw BadArgumentError( "genfile::vcf::FloatType::parse()", "value = \"" + std::string( value ) + "\"" ) ;
			}
		}

		void PhredScaleFloatType::parse( std::size_t i, string_utils::slice const& value, EntrySetter& setter ) const {
			try {
				setter.set_value( i, std::pow( 10, -string_utils::strtod( value ) / 10 ) ) ;
			}
			catch( string_utils::StringConversionError const& ) {
				throw BadArgumentError( "genfile::vcf::PhredScaleFloatType::parse()", "value = \"" + std::string( value ) + "\"" ) ;
			}
		}

		void CharacterType::parse( std::size_t i, string_utils::slice const& value, EntrySetter& setter ) const {
			if( value.size() != 1 ) {
				throw BadArgumentError("genfile::vcf::CharacterType::parse()", "value = \"" + std::string( value ) +"\"" ) ;
			}
			std::string string_value = value ;
			setter.set_value( i, string_value ) ;
		}

		void FlagType::parse( std::size_t i, string_utils::slice const& value, EntrySetter& setter ) const {
			assert(0) ;
		}
		
		std::string VCFEntryType::m_missing_value = "." ;
		
		VCFEntryType::UniquePtr VCFEntryType::create( Spec const& spec ) {
			Spec::const_iterator ID = spec.find( "ID" ) ;
			Spec::const_iterator number = spec.find( "Number" ) ;
			Spec::const_iterator type_i = spec.find( "Type" ) ;
			if( ID == spec.end() || number == spec.end() || type_i == spec.end() ) {
				throw BadArgumentError( "genfile::vcf::VCFEntryType::create()", "spec" ) ;
			}

			VCFEntryType::UniquePtr result ;

			if( ID->second == "GT" || type_i->second == "Genotype" ) {
				if( number->second != "." && number->second != "1" ) {
					throw BadArgumentError(
						"genfile::vcf::VCFEntryType::create()",
						"Number=" + number->second,
						"\"Number\" of GT field must be encoded as . or 1 in metadata."
					) ;
				}
				if( type_i->second != "String" && type_i->second != "Genotype" ) {
					throw BadArgumentError( "genfile::vcf::VCFEntryType::create()", "spec" ) ;
				}
				result.reset( new GenotypeCallVCFEntryType() ) ;
			}
			else if( ID->second == "GP" ) {
				if( number->second != "." && ( number->second == "A" || number->second == "R" ) ) {
					throw BadArgumentError(
						"genfile::vcf::VCFEntryType::create()",
						"Number=" + number->second,
						"\"Number\" of GP field must be encoded as ., G, or a number in metadata (ideally G)."
					) ;
				}
				if( type_i->second != "Float" && type_i->second != "Integer" ) {
					throw BadArgumentError( "genfile::vcf::VCFEntryType::create()", "spec" ) ;
				}

				SimpleType::UniquePtr element_type = SimpleType::create( "Probability", "identity" ) ;
				result.reset( new OnePerGenotypeVCFEntryType( element_type )) ;
			} else {
				SimpleType::UniquePtr element_type = SimpleType::create( type_i->second, "identity" ) ;

				if( number->second == "A" ) {
					result.reset( new OnePerAlternateAlleleVCFEntryType( element_type )) ;
				}
				else if( number->second == "R" ) {
					result.reset( new OnePerAlleleVCFEntryType( element_type )) ;
				}
				else if( number->second == "G" ) {
					result.reset( new OnePerGenotypeVCFEntryType( element_type )) ;
				}
				else if( number->second == "." ) {
					result.reset(
						new DynamicNumberVCFEntryType( element_type )
					) ;
				}
				else {
					result.reset(
						new FixedNumberVCFEntryType(
							string_utils::to_repr< std::size_t >( number->second ),
							element_type
						)
					) ;
				}
			}

			return result ;
		}
		
		namespace impl {
			void lex(
				 string_utils::slice const& value,
				 std::size_t number_of_alleles,
				 char const* separator,
				 std::pair< std::size_t, std::size_t > const& value_count_range,
				 std::vector< string_utils::slice >* result
			) {
				assert( result ) ;
				result->clear() ;
				// An empty value is treated as an empty list, (not a list with one empty value)
				if( !value.empty() ) {
					string_utils::slice( value ).split( separator, result ) ;
				}
				if( result->size() < value_count_range.first || result->size() > value_count_range.second ) {
					throw BadArgumentError( "genfile::vcf::()::lex()", "value = \"" + std::string( value ) + "\"" ) ;
				}
			}

			void parse_elts(
				std::vector< string_utils::slice > const& elts,
				uint32_t ploidy,
				OrderType order_type,
				SimpleType const& value_type,
				std::string const& missing_value,
				EntriesSetter& setter
			) {
				setter.set_number_of_entries( ploidy, elts.size(), order_type, value_type.represented_type() ) ;
				for( std::size_t i = 0; i < elts.size(); ++i ) {
					if( elts[i] == missing_value ) {
						setter.set_value( i, MissingValue() ) ;
					}
					else {
						value_type.parse( i, elts[i], setter ) ;
					}
				}
			}


			void get_missing_value(
				std::size_t number_of_alleles,
				uint32_t ploidy,
				OrderType order_type,
				SimpleType const& value_type,
				std::pair< std::size_t, std::size_t > const& value_count_range,
				EntriesSetter& setter
			) {
				setter.set_number_of_entries( ploidy, value_count_range.first, order_type, value_type.represented_type() ) ;
				for( std::size_t i = 0; i < value_count_range.first; ++i ) {
					setter.set_value( i, MissingValue() ) ;
				}
			}
		}
			
		
		VCFEntryType::VCFEntryType( SimpleType::UniquePtr type ):
			m_type( type )
		{
		}
		
		void ListVCFEntryType::parse( string_utils::slice const& value, std::size_t number_of_alleles, uint32_t ploidy, EntriesSetter& setter ) const {
			std::vector< string_utils::slice > result ;
			impl::lex( value, number_of_alleles, ",", get_value_count_range( number_of_alleles, ploidy ), &result ) ;
			impl::parse_elts( result, ploidy, get_order_type(), get_value_type(), missing_value(), setter ) ;
		}

		void ListVCFEntryType::get_missing_value( std::size_t number_of_alleles, uint32_t ploidy, EntriesSetter& setter ) const {
			impl::get_missing_value(
				number_of_alleles,
				ploidy,
				get_order_type(),
				get_value_type(),
				get_value_count_range( number_of_alleles, ploidy ),
				setter
			) ;
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
			
			uint32_t n_choose_k( uint32_t n, uint32_t k ) {
				if( k == 0 )  {
					return 1 ;
				} else if( k == 1 ) {
					return n ;
				}
				return ( n * n_choose_k(n - 1, k - 1) ) / k ;
			}
			
			std::size_t get_number_of_unphased_genotypes( std::size_t const n_alleles, std::size_t const ploidy ) {
				// The number is equal to the number of ways to fill an n-vector
				// (where n is the number of alleles)
				// with nonnegative integers so that the sum is the ploidy.
				// This has a recursive expression involving filling the first
				// entry and then filling the others.
				// There are two special cases: if the ploidy is zero, there are no genotypes at all,
				// while if there is only one allele then there is only one genotype.
				
				if( ploidy == 0 ) {
					return 0 ;
				} else {
					return n_choose_k( ploidy + n_alleles - 1, n_alleles - 1 ) ;
				}
			}

			std::size_t get_number_of_phased_genotypes( std::size_t n_alleles, uint32_t ploidy ) {
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
		
		ListVCFEntryType::ValueCountRange OnePerGenotypeVCFEntryType::get_value_count_range( std::size_t number_of_alleles, uint32_t ploidy ) const {
			if( ploidy != eUnknownPloidy && number_of_alleles == 0 && ploidy > 0 ) {
				throw BadArgumentError( "genfile::vcf::OnePerGenotypeVCFEntryType::get_value_count_range()", "number_of_alleles = 0" ) ;
			}
			std::size_t N = ( ploidy == 0 ) ? 0 : impl::get_number_of_unphased_genotypes( number_of_alleles, ploidy ) ;
			//std::cerr << "number_of_alleles = " << number_of_alleles << ", ploidy = " << ploidy << ", N = " << N << ".\n" ;
			return ValueCountRange( N, N ) ;
		}
		
		GenotypeCallVCFEntryType::GenotypeCallVCFEntryType():
			VCFEntryType( SimpleType::create( "AlleleIndex" )),
			m_missing_value( "." )
		{}
		
		void GenotypeCallVCFEntryType::lex( string_utils::slice const& value, std::vector< string_utils::slice >* result ) const {
			assert( result ) ;
			std::vector< string_utils::slice > elts ;
			// empty value is treated as empty list.
			if( !value.empty() ) {
				string_utils::slice( value ).split( "|/", result ) ;
			}
		}

		void GenotypeCallVCFEntryType::parse( string_utils::slice const& value, std::size_t number_of_alleles, uint32_t ploidy, EntriesSetter& setter ) const {
			return parse_impl( value, number_of_alleles, setter ) ;
		}

		// Return a set of missing values into the setter object.
		void GenotypeCallVCFEntryType::get_missing_value( std::size_t number_of_alleles, uint32_t ploidy, EntriesSetter& setter ) const {
			assert( ploidy != eUnknownPloidy ) ;
			setter.set_number_of_entries( ploidy, ploidy, ePerUnorderedHaplotype, eAlleleIndex ) ;
			for( std::size_t i = 0; i < ploidy; ++i ) {
				setter.set_value( i, MissingValue() ) ;
			}
		}

		namespace impl {
			//
			// 
			struct RangeCheckedGTSetter: public EntriesSetter {
				RangeCheckedGTSetter( EntriesSetter& setter, EntriesSetter::Integer const max_genotype ):
					m_setter( setter ),
					m_max_genotype( max_genotype )
				{}
				
				void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
					assert( order_type == ePerOrderedHaplotype || order_type == ePerUnorderedHaplotype ) ;
					m_setter.set_number_of_entries( n, n, order_type, value_type ) ;
				}
				void set_value( std::size_t i, MissingValue const value ) { m_setter.set_value( i, value ) ; }
				void set_value( std::size_t i, Integer const value ) {
					if( value < 0 || value > m_max_genotype ) {
						throw BadArgumentError( "genfile::vcf::impl::RangeCheckedGTSetter::set_value", "value = \"" + string_utils::to_string( value ) + "\"" ) ;
					}
					m_setter.set_value( i, value ) ;
				}
			private:
				EntriesSetter& m_setter ;
				EntriesSetter::Integer m_max_genotype ;
			} ;
		}

		void GenotypeCallVCFEntryType::parse_impl(
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

			// Most GT values have one character per allele.
			// We treat this as a special case because it is much quicker
			bool simple_parse_success = true ;

			if( value.size() == 3 && ( value[1] == '|' || value[1] == '/' ) ) {
				char a, b ;
				a = value[0] - '0' ;
				b = value[2] - '0' ;
				if((( value[0] == m_missing_value[0] ) || ( a >= 0 && a <= max ) ) && ( ( value[2] == m_missing_value[0] ) || ( b >= 0 && b <= max ) ) ) {
					setter.set_number_of_entries( 2, 2, ( value[1] == '|' ) ? ePerOrderedHaplotype : ePerUnorderedHaplotype, eAlleleIndex ) ;
					if( value[0] == m_missing_value[0] ) {
						setter.set_value( 0, MissingValue() ) ;
					} else {
						setter.set_value( 0, EntriesSetter::Integer( a ) ) ;
					}
					if( value[2] == m_missing_value[0] ) {
						setter.set_value( 1, MissingValue() ) ;
					} else {
						setter.set_value( 1, EntriesSetter::Integer( b ) ) ;
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
					OrderType const order_type = ( (value.size() == 1) || (value[1] == '|') ) ? ePerOrderedHaplotype : ePerUnorderedHaplotype ;
					setter.set_number_of_entries( simple_values.size(), simple_values.size(), order_type, eAlleleIndex ) ;
					for( std::size_t i = 0; i < simple_values.size(); ++i ) {
						if( simple_values[i] == -1 ) {
							setter.set_value( i, MissingValue() ) ;
						} else {
							setter.set_value( i, EntriesSetter::Integer( simple_values[i] )) ;
						}
					}
				}
			}
			else {
				simple_parse_success = false ;
			}

			if( !simple_parse_success ) {
				impl::RangeCheckedGTSetter checked_genotype_setter( setter, number_of_alleles - 1 ) ;
				std::vector< string_utils::slice > elts ;
				lex( value, &elts ) ;
				OrderType const order_type = ( value.find( '/' ) == std::string::npos )
					? ePerOrderedHaplotype
					: ePerUnorderedHaplotype ; 
				impl::parse_elts( elts, elts.size(), order_type, get_value_type(), m_missing_value, checked_genotype_setter ) ;
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
