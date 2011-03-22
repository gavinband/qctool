#include <vector>
#include <memory>
#include <set>
#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/string/replace.hpp>
#include "genfile/vcf/CallReader.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils.hpp"
#include "test_case.hpp"

using namespace genfile::vcf ;
using namespace genfile::string_utils ;

void test_format_types() ;
void test_malformed_format() ;

struct GenotypeCallChecker
// struct GenotypeCallChecker
// checks the output of CallReader against a specified set of expected calls.
//
// Be careful: this checks both that the values set are the ones expected,
// and, on destruction, that all values have been set.
// You must therefore avoid copying using e.g. boost::cref().
{
	GenotypeCallChecker( std::vector< std::vector< Entry > > const& expected_calls ):
		m_calls( expected_calls )
	{
	}

	GenotypeCallChecker( GenotypeCallChecker const& other ):
		m_calls( other.m_calls ),
		m_indices_of_set_values( other.m_indices_of_set_values )
	{
	}
	
	// Destructor: verify that all expected values were set.
	~GenotypeCallChecker() {
		TEST_ASSERT( m_indices_of_set_values.size() == m_calls.size() ) ;
		for( std::size_t i = 0; i < m_calls.size(); ++i ) {
			TEST_ASSERT( m_indices_of_set_values.find( i ) != m_indices_of_set_values.end() ) ;
		}
	}

	// Operator(): verify that the correct value is passed
	// (and for an in-range individual.)
	void operator()( std::size_t i, std::vector< Entry > const& values ) const {
		TEST_ASSERT( i < m_calls.size() ) ;
		TEST_ASSERT( m_indices_of_set_values.find( i ) == m_indices_of_set_values.end() ) ;
		TEST_ASSERT( values.size() == m_calls[i].size() ) ;
		for( std::size_t j = 0; j < values.size(); ++j ) {
			TEST_ASSERT( values[j] == m_calls[i][j] ) ;
		}
		m_indices_of_set_values.insert( i ) ;
	}

	private:
		std::vector< std::vector< Entry > > const m_calls ;
		mutable std::set< std::size_t > m_indices_of_set_values ;

		GenotypeCallChecker& operator=( GenotypeCallChecker const& other ) ;
} ;

struct Ignore
{
	void operator()( std::size_t, std::vector< Entry > const& ) const {}
} ;


AUTO_TEST_CASE( test_format ) {
	std::cerr << "test_format..." ;
	test_format_types() ;
	test_malformed_format() ;
	std::cerr << "ok.\n" ;
}

void test_format_types() {
	boost::ptr_map< std::string, VCFEntryType > types ;

	// Test with no types

	try {
		CallReader( 2, "", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "HQ", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "GT:HQ", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "HQ:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	// Test with GT type only.
	{
		std::string const GT = "GT" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "Integer" ;
		spec[ "Description" ] = "A test" ;
		types.insert( GT, VCFEntryType::create( spec )) ;
	}

	try {
		CallReader( 2, "", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "GT", "", types ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
		TEST_ASSERT( 0 ) ;
	}

	try {
		CallReader( 2, "GT:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "HQ", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "GT:HQ", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}
	
	try {
		CallReader( 2, "HQ:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	// Test with HQ type only.

	types.clear() ;
	{
		std::string const ID = "HQ" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = ID ;
		spec[ "Number" ] = "5" ;
		spec[ "Type" ] = "String" ;
		spec[ "Description" ] = "A test, HQ" ;
		types.insert( ID, VCFEntryType::create( spec )) ;
	}
	
	try {
		CallReader( 2, "", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
	}

	try {
		CallReader( 2, "GT:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "HQ", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "GT:HQ", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
	}
	
	try {
		CallReader( 2, "HQ:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}
	
	// Test with GT and HQ type
	{
		std::string const GT = "GT" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "Integer" ;
		spec[ "Description" ] = "A test" ;
		types.insert( GT, VCFEntryType::create( spec )) ;
	}
	
	try {
		CallReader( 2, "", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "GT", "", types ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
		TEST_ASSERT( 0 ) ;
	}

	try {
		CallReader( 2, "GT:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "HQ", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "GT:HQ", "", types ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
		TEST_ASSERT( 0 ) ;
	}
	
	try {
		CallReader( 2, "HQ:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}
	
	// Test with three types
	{
		std::string const ID = "another_type" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = ID ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "String" ;
		spec[ "Description" ] = "Another test type" ;
		types.insert( ID, VCFEntryType::create( spec )) ;
	}
	
	try {
		CallReader( 2, "", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "GT", "", types ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
		TEST_ASSERT( 0 ) ;
	}

	try {
		CallReader( 2, "GT:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "HQ", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 2, "GT:HQ", "", types ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
		TEST_ASSERT( 0 ) ;
	}
	
	try {
		CallReader( 2, "HQ:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}
	
	try {
		CallReader( 2, "GT:another_type", "", types ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
		TEST_ASSERT( 0 ) ;
	}

	try {
		CallReader( 2, "GT:HQ:another_type", "", types ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
		TEST_ASSERT( 0 ) ;
	}

	try {
		CallReader( 2, "GT:another_type:HQ", "", types ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
		TEST_ASSERT( 0 ) ;
	}

	try {
		CallReader( 2, "another_type:GT:HQ", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
	}

	try {
		CallReader( 2, "GT:another_type:HQ:another_type", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {
	}
}

AUTO_TEST_CASE( test_setters ) {
	std::cerr << "test_setters()..." ;

	using namespace genfile ;
	
	boost::ptr_map< std::string, VCFEntryType > types ;
	{
		std::string const GT = "GT" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "Integer" ;
		spec[ "Description" ] = "A test" ;
		types.insert( GT, VCFEntryType::create( spec )) ;
	}
	{
		std::string const ID = "HQ" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = ID ;
		spec[ "Number" ] = "5" ;
		spec[ "Type" ] = "String" ;
		spec[ "Description" ] = "A test, HQ" ;
		types.insert( ID, VCFEntryType::create( spec )) ;
	}
	{
		std::string const ID = "another_type" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = ID ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "String" ;
		spec[ "Description" ] = "Another test type" ;
		types.insert( ID, VCFEntryType::create( spec )) ;
	}
	
	try {
		CallReader( 2, "GT", "1|1", types ) ;
	}
	catch( BadArgumentError const& ) {
		TEST_ASSERT(0) ;
	}
	catch( MalformedInputError const& ) {
		TEST_ASSERT(0) ;
	}

	try {
		std::vector< std::vector< Entry > > expected ;
		expected.push_back( std::vector< Entry >( 2, Entry( 1 ))) ;
		CallReader( 2, "GT", "1|1", types )
			( "GT", boost::cref( GenotypeCallChecker( expected ))) ;
	}
	catch( BadArgumentError const& ) {
		TEST_ASSERT(0) ;
	}
	catch( MalformedInputError const& ) {
		TEST_ASSERT(0) ;
	}

	try {
		std::vector< std::vector< Entry > > expected ;
		expected.push_back( std::vector< Entry >( 5, MissingValue() )) ;
		CallReader( 2, "GT", "1|1", types )
			( "HQ", boost::cref( GenotypeCallChecker( expected ))) ;
	}
	catch( BadArgumentError const& ) {
		TEST_ASSERT(0) ;
	}
	catch( MalformedInputError const& ) {
		TEST_ASSERT(0) ;
	}

	try {
		std::vector< std::vector< Entry > > expected ;
		expected.push_back( std::vector< Entry >( 5, MissingValue() )) ;
		CallReader( 2, "GT:HQ", "1|1", types )
			( "HQ", boost::cref( GenotypeCallChecker( expected ))) ;
	}
	catch( BadArgumentError const& ) {
		TEST_ASSERT(0) ;
	}
	catch( MalformedInputError const& ) {
		TEST_ASSERT(0) ;
	}

	try {
		std::vector< std::vector< Entry > > expected ;
		expected.push_back( std::vector< Entry >( 5, MissingValue() )) ;
		CallReader( 2, "GT", "1|1", types )
			( "GT", Ignore() )
			( "HQ", boost::cref( GenotypeCallChecker( expected ))) ;
	}
	catch( BadArgumentError const& ) {
		TEST_ASSERT(0) ;
	}
	catch( MalformedInputError const& ) {
		TEST_ASSERT(0) ;
	}

	try {
		std::vector< std::vector< Entry > > expected ;
		expected.push_back( std::vector< Entry >( 1, MissingValue() )) ;
		{
			CallReader( 2, "GT:HQ", "1|1", types )
				( "another_type", boost::cref( GenotypeCallChecker( expected ))) ;
		}
		expected.push_back( std::vector< Entry >( 1, MissingValue() )) ;
		{
			CallReader( 2, "GT:HQ", "1|1\t1|1", types )
				( "another_type", boost::cref( GenotypeCallChecker( expected ))) ;
		}
	}
	catch( BadArgumentError const& ) {
		TEST_ASSERT(0) ;
	}
	catch( MalformedInputError const& ) {
		TEST_ASSERT(0) ;
	}
	std::cerr << "ok.\n" ;
}

void test_malformed_format() {
	boost::ptr_map< std::string, VCFEntryType > types ;
	{
		std::string const GT = "GT" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "Integer" ;
		spec[ "Description" ] = "A test" ;
		types.insert( GT, VCFEntryType::create( spec )) ;
	}
	{
		std::string const ID = "another_type" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = ID ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "String" ;
		spec[ "Description" ] = "Another test type" ;
		types.insert( ID, VCFEntryType::create( spec )) ;
	}
	{
		std::string const ID = "another_type" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = ID ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "String" ;
		spec[ "Description" ] = "Another test type" ;
		types.insert( ID, VCFEntryType::create( spec )) ;
	}

	std::string format = "GT:HQ:another_type" ;
	for( int c = 0; c < std::numeric_limits< char >::max(); ++c ) {
		for( std::size_t i = 0; i < format.size(); ++i ) {
			std::string malformed_format = format.substr( 0, i ) + std::string( 1, c ) + format.substr( i + 1, format.size() ) ;
			try {
				CallReader( 2, "GT:HQ:another_type", "", types ) ;
				TEST_ASSERT( 0 ) ;
			}
			catch( genfile::BadArgumentError const& e ) {
			}
		}
	}
}

AUTO_TEST_CASE( test_alleles ) {
	std::cerr << "test_alleles..." ;
	boost::ptr_map< std::string, VCFEntryType > types ;
	{
		std::string const GT = "GT" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "Integer" ;
		spec[ "Description" ] = "A test" ;
		types.insert( GT, VCFEntryType::create( spec )) ;
	}
	{
		std::string const AT = "AT" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = AT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "String" ;
		spec[ "Description" ] = "Another test type" ;
		types.insert( AT, VCFEntryType::create( spec )) ;
	}
	
	for( std::size_t n = 0; n < 100; ++n ) {
		try {
			CallReader( n, "GT:AT", "", types ) ;
			if( n == 0 ) {
				TEST_ASSERT( 0 ) ;
			}
		}
		catch( genfile::BadArgumentError const& e ) {
			if( n > 0 ) {
				TEST_ASSERT( 0 ) ;
			}
		}

		try {
			CallReader( n, "GT", "", types ) ;
			if( n == 0 ) {
				TEST_ASSERT( 0 ) ;
			}
		}
		catch( genfile::BadArgumentError const& e ) {
			if( n > 0 ) {
				TEST_ASSERT( 0 ) ;
			}
		}
	}
	std::cerr << "ok.\n" ;
}

std::pair< std::string, std::vector< std::vector< Entry > > > make_call_data(
	std::size_t const number_of_individuals,
	std::size_t const number_of_alleles,
	std::size_t const ploidy,
	std::string const& sep
) {
	std::vector< std::vector< Entry > > calls( number_of_individuals ) ;
	std::ostringstream ostr ;
	for( std::size_t i = 0; i < number_of_individuals; ++i ) {
		if( i > 0 ) {
			ostr << "\t" ;
		}
		for( std::size_t p = 0; p < ploidy; ++p ) {
			if( p > 0 ) {
				ostr << sep ;
			}
			// Make calls with separation equal to i (the individual id)
			// And first equal to the number of individuals modulo the number of alleles.
			calls[i].push_back( int( ( number_of_individuals + ( p * i ) ) % number_of_alleles ) ) ;
			ostr << calls[i].back() ;
		}
		ostr << ":." ;
	}
	return std::make_pair( ostr.str() , calls ) ;
}

std::auto_ptr< boost::ptr_map< std::string, VCFEntryType > > make_some_types() {
	boost::ptr_map< std::string, VCFEntryType > types ;
	{
		std::string const GT = "GT" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "Integer" ;
		spec[ "Description" ] = "A test" ;
		types.insert( GT, VCFEntryType::create( spec )) ;
	}
	{
		std::string const AT = "AT" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = AT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "String" ;
		spec[ "Description" ] = "Another test type" ;
		types.insert( AT, VCFEntryType::create( spec )) ;
	}
	return types.release() ;
}

AUTO_TEST_CASE( test_simple_gt_values ) {
	std::cerr << "test_simple_gt_values..." ;
	boost::ptr_map< std::string, VCFEntryType > types( make_some_types() ) ;

	// A simple example.
	for( int c1 = 0; c1 < std::numeric_limits< char >::max(); ++c1 ) {
		for( int c2 = 0; c2 < std::numeric_limits< char >::max(); ++c2 ) {
			if( c1 == '\t' ) ++c1 ;
			if( c2 == '\t' ) ++c2 ;
			bool bad_c1 = ( ( c1 != '.' ) && ( c1 < '0' || c1 > '4' ) ) ;
			bool bad_c2 = ( ( c2 != '.' ) && ( c2 < '0' || c2 > '4' ) ) ;

			try {
				std::string data = ".|." ;
				data[0] = c1 ;
				data[2] = c2 ;
				
				std::vector< std::vector< Entry > > expected_calls( 1, std::vector< Entry >( 2 )) ;

				if( bad_c1 || bad_c2 ) {
					// no calls should be set.
					expected_calls.clear() ;
				}
				else {
					// missing calls
					if( c1 != '.' ) {
						expected_calls[0][0] = c1 - '0' ;
					}
					if( c2 != '.' ) {
						expected_calls[0][1] = c2 - '0' ;
					}
				}

				{
					// std::cerr << data << ".\n" ;
					GenotypeCallChecker checker( expected_calls ) ;
					CallReader( 5, "GT", data, types )
						( "GT", boost::ref( checker ) ) ;
				}
				
				if( bad_c1 || bad_c2 ) {
					// std::cerr << "c1 = " << c1 << "('" << char(c1) << "')" << ", c2 = " << c2 << "('" << char(c2) << "').\n" ;
					TEST_ASSERT(0) ;
				}
			}
			catch( genfile::MalformedInputError const& ) {
				if( !bad_c1 && !bad_c2 ) {
					// std::cerr << "c1 = " << c1 << "('" << char(c1) << "')" << ", c2 = " << c2 << "('" << char(c2) << "').\n" ;
					TEST_ASSERT(0) ;
				}
			}
		}
	}

	// A more complex example.
	for( int c1 = 0; c1 < std::numeric_limits< char >::max(); ++c1 ) {
		for( int c2 = 0; c2 < std::numeric_limits< char >::max(); ++c2 ) {
			if( c1 == '\t' ) ++c1 ;
			if( c2 == '\t' ) ++c2 ;
			
			bool bad_c1 = ( ( c1 != '.' ) && ( c1 < '0' || c1 > '4' ) ) ;
			bool bad_c2 = ( ( c2 != '.' ) && ( c2 < '0' || c2 > '4' ) ) ;

			try {
				std::string data = "1|1\t1|" + std::string( 1, c2 ) + "\t" + std::string( 1, c1 ) + "|" + std::string( 1, c2 ) ;
				CallReader( 5, "GT", data, types ) ;
				if( bad_c1 || bad_c2 ) {
					// std::cerr << "c1 = " << c1 << "('" << char(c1) << "')" << ", c2 = " << c2 << "('" << char(c2) << "').\n" ;
					TEST_ASSERT(0) ;
				}
			}
			catch( genfile::MalformedInputError const& ) {
				if( !bad_c1 && !bad_c2 ) {
					TEST_ASSERT(0) ;
				}
			}
		}
	}

	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_omitted_values ) {
	std::cerr << "test_omitted_values..." ;
	boost::ptr_map< std::string, VCFEntryType > types( make_some_types() ) ;
	// types:
	// GT
	// AT - Number=1, Type=String
	// 
	using namespace genfile ;
	
	try {
		// CallReader will not parse the field unless asked to
		CallReader( 2, "GT:AT", "1|1:", types ) ;
	}
	catch( genfile::MalformedInputError const& ) {
		TEST_ASSERT( 0 ) ;
	}

	try {
		CallReader( 2, "GT:AT", "1|1:", types )
			( "AT", Ignore() ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::MalformedInputError const& ) {
	}

	try {
		CallReader( 2, "GT:AT", "1|1", types ) ;
		
		std::vector< std::vector< Entry > > expected_values ;
		expected_values.push_back( std::vector< Entry >( 1, MissingValue() )) ;
		CallReader( 2, "GT:AT", "1|1", types )
			( "AT", boost::cref( GenotypeCallChecker( expected_values ))) ;
	}
	catch( genfile::MalformedInputError const& ) {
		TEST_ASSERT(0) ;
	}

	try {
		CallReader( 2, "GT:AT", "1|1:hello\t0|1", types ) ;

		std::vector< std::vector< Entry > > expected_values ;
		expected_values.push_back( std::vector< Entry >( 1, std::string( "hello" ) )) ;
		expected_values.push_back( std::vector< Entry >( 1, MissingValue() )) ;
		CallReader( 2, "GT:AT", "1|1:hello\t0|1", types )
			( "AT", boost::cref( GenotypeCallChecker( expected_values ))) ;
	}
	catch( genfile::MalformedInputError const& ) {
		TEST_ASSERT(0) ;
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_gt_delimiter ) {
	std::cerr << "test_gt_delimiter..." ;
	boost::ptr_map< std::string, VCFEntryType > types( make_some_types() ) ;

	// test the delimiter between genotypes.

	for( std::size_t number_of_individuals = 0; number_of_individuals < 10; ++number_of_individuals ) {
		
		std::string data ;
		std::vector< std::vector< Entry > > expected_calls ;
		boost::tie( data, expected_calls ) = make_call_data( number_of_individuals, 5, 2, "|" ) ;
		for( int c = 0; c < std::numeric_limits< char >::max(); ++c ) {
			while( char( c ) == '|' || char( c ) == '/' || char( c ) == '\t' || char( c ) == '0' ) {
				++c ;
			}

			try {
				std::string used_data = data ;
				for( std::size_t i = 0; i < number_of_individuals; ++i ) {
					used_data[ (i*6) + 1 ] = char( c ) ;
				}

				CallReader( 5, "GT:AT", used_data, types ) ;

				if( number_of_individuals > 0 ) {
					// std::cerr << "data = \"" << data << "\".\n" ;
					// std::cerr << "char = " << c << "('" << char( c ) << "'), data = \"" << used_data << "\".\n" ;

					TEST_ASSERT(0) ;
				}
			}
			catch( genfile::MalformedInputError const& ) {
				if( number_of_individuals == 0 ) {
					TEST_ASSERT(0) ;
				}
			}
			catch( genfile::BadArgumentError const& ) {
				TEST_ASSERT(0) ;
			}
		}
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_complex_gt_values ) {
	std::cerr << "test_complex_gt_values..." ;
	boost::ptr_map< std::string, VCFEntryType > types( make_some_types() ) ;

	for( std::size_t number_of_individuals = 0; number_of_individuals < 25; ++number_of_individuals ) {
		// std::cerr << number_of_individuals << " " ;
		for( std::size_t number_of_alleles = 1; number_of_alleles < 10; ++number_of_alleles ) {
			for( std::size_t ploidy = 0; ploidy < 10; ++ploidy ) {
				std::vector< std::vector< Entry > > expected_calls ;
				std::string phased_data, unphased_data ;
				boost::tie( phased_data, expected_calls ) = make_call_data( number_of_individuals, number_of_alleles, ploidy, "|" ) ;
				boost::tie( unphased_data, expected_calls ) = make_call_data( number_of_individuals, number_of_alleles, ploidy, "/" ) ;

				try {
					{
						// std::cerr << "here:" << number_of_individuals << " " << number_of_alleles << " " << ploidy << ": \"" << phased_data << "\"\n" ;
						CallReader( number_of_alleles, "GT:AT", phased_data, types ) ;
						GenotypeCallChecker checker( expected_calls ) ;
						CallReader( number_of_alleles, "GT:AT", phased_data, types )
						 	( "GT", boost::ref( checker )) ;
					}
					{
						CallReader( number_of_alleles, "GT:AT", unphased_data, types ) ;
						GenotypeCallChecker checker( expected_calls ) ;
						CallReader( number_of_alleles, "GT:AT", unphased_data, types )
						 	( "GT", boost::ref( checker )) ;
					}
				}
				catch( genfile::MalformedInputError const& e ) {
					TEST_ASSERT(0) ;
				}
				catch( genfile::BadArgumentError const& e ) {
					std::cerr << e.what() << ":" << e.function() << ":" << e.arguments() << ".\n" ;
					TEST_ASSERT(0) ;
				}
			}
		}
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_gt_genotype_bounds ) {
	std::cerr << "test_gt_genotype_bounds..." ;
	boost::ptr_map< std::string, VCFEntryType > types( make_some_types() ) ;

	for( int n_alleles = 1; n_alleles < 100; ++n_alleles ) {
		for( int call = 0; call < 110; ++call ) {
			bool bad_call = ( call >= n_alleles ) ;
			std::vector< std::vector< Entry > > expected_set_calls ;
			if( !bad_call ) {
				expected_set_calls.push_back( std::vector< Entry >( 1, call )) ;
			}
			try {
				GenotypeCallChecker checker( expected_set_calls ) ;
				CallReader(
					n_alleles,
					"GT",
					::to_string( call ),
					types
				)
				( "GT", boost::ref( checker )) ;
				
				if( bad_call ) {
					TEST_ASSERT(0) ;
				}
			}
			catch( genfile::MalformedInputError const& ) {
				if( !bad_call ) {
					TEST_ASSERT(0) ;
				}
			}
			catch( genfile::BadArgumentError const& ) {
				TEST_ASSERT(0) ;
			}
		}
	}

	for( int n_alleles = 1; n_alleles < 100; ++n_alleles ) {
		for( int call = 0; call < 110; ++call ) {
			std::string data = to_string( std::min( call, n_alleles - 1 ) ) + "\t" + to_string( call ) ;
			bool bad_call = ( call >= n_alleles ) ;
			std::vector< std::vector< Entry > > expected_set_calls ;
			expected_set_calls.push_back( std::vector< Entry >( 1, std::min( call, n_alleles - 1 ) ) ) ;
			if( !bad_call ) {
				expected_set_calls.push_back( std::vector< Entry >( 1, call )) ;
			}
			try {
				GenotypeCallChecker checker( expected_set_calls ) ;
				CallReader(
					n_alleles,
					"GT",
					data,
					types
				)
				( "GT", boost::ref( checker )) ;
				
				if( bad_call ) {
					TEST_ASSERT(0) ;
				}
			}
			catch( genfile::MalformedInputError const& ) {
				if( !bad_call ) {
					TEST_ASSERT(0) ;
				}
			}
			catch( genfile::BadArgumentError const& ) {
				TEST_ASSERT(0) ;
			}
		}
	}

	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_float_fields ) {
	std::cerr << "test_float_fields..." ;
	
	// verify that we can read float fields out of the data correctly.

	boost::ptr_map< std::string, VCFEntryType > types( make_some_types() ) ;
	{
		std::string const ID = "F1" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = ID ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "Float" ;
		spec[ "Description" ] = "" ;
		types.insert( ID, VCFEntryType::create( spec )) ;
	}
	{
		std::string const ID = "F2" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = ID ;
		spec[ "Number" ] = "3" ;
		spec[ "Type" ] = "Float" ;
		spec[ "Description" ] = "" ;
		types.insert( ID, VCFEntryType::create( spec )) ;
	}

	{
		try {
			std::vector< std::vector< Entry > > expected_result ;
			CallReader( 2, "GT:F1", "", types )
				( "F1", boost::cref( GenotypeCallChecker( expected_result )) ) ;
		}
		catch( genfile::MalformedInputError const& ) {
			TEST_ASSERT(0) ;
		}

		try {
			std::vector< std::vector< Entry > > expected_result ;
			expected_result.push_back( std::vector< Entry >( 1, 0.1 )) ;
			CallReader( 2, "GT:F1", "1|1:0.1", types )( "F1", boost::cref( GenotypeCallChecker( expected_result ) )) ;
		}
		catch( genfile::MalformedInputError const& ) {
			TEST_ASSERT(0) ;
		}

		try {
			std::vector< std::vector< Entry > > expected_result ;
			expected_result.push_back( std::vector< Entry >( 1, 0.1 )) ;
			expected_result.push_back( std::vector< Entry >( 1, 1.1 )) ;
			expected_result.push_back( std::vector< Entry >( 1, 2.1 )) ;
			CallReader( 2, "GT:F1", "1|1:0.1\t0|1:1.1\t0|0:2.1", types )
			 	( "F1", boost::cref( GenotypeCallChecker( expected_result ) )) ;
		}
		catch( genfile::MalformedInputError const& ) {
			TEST_ASSERT(0) ;
		}

		try {
			std::vector< std::vector< Entry > > expected_result ;
			expected_result.push_back( std::vector< Entry >( 1, 0.1 )) ;
			expected_result.push_back( std::vector< Entry >( 1, genfile::MissingValue() )) ;
			expected_result.push_back( std::vector< Entry >( 1, 2.1 )) ;

			CallReader( 2, "GT:F1", "1|1:0.1\t0|1:.\t0|0:2.1", types )
			 	( "F1", boost::cref( GenotypeCallChecker( expected_result ) )) ;
		}
		catch( genfile::MalformedInputError const& ) {
			TEST_ASSERT(0) ;
		}

		try {
			CallReader( 2, "GT:F2", "", types )
			 	( "F2", boost::cref( GenotypeCallChecker( std::vector< std::vector< Entry > >() ) )) ;
		}
		catch( genfile::MalformedInputError const& ) {
			TEST_ASSERT(0) ;
		}

		try {
			std::vector< std::vector< Entry > > expected_result ;
			expected_result.push_back( std::vector< Entry >()) ;
			expected_result.back().push_back( 0.1 ) ;
			expected_result.back().push_back( -1.5568 ) ;
			expected_result.back().push_back( -1e-08 ) ;
			
			CallReader( 2, "GT:F2", "1|1:0.1,-1.5568,-1e-08", types )
			 	( "F2", boost::cref( GenotypeCallChecker( expected_result )) ) ;
		}
		catch( genfile::MalformedInputError const& ) {
			TEST_ASSERT(0) ;
		}

		try {
			std::vector< std::vector< Entry > > expected_result ;
			expected_result.push_back( std::vector< Entry >() ) ;
			expected_result.back().push_back( 0.1 ) ;
			expected_result.back().push_back( -1.5568 ) ;
			expected_result.back().push_back( -1e-08 ) ;
			expected_result.push_back( std::vector< Entry >() ) ;
			expected_result.back().push_back( genfile::MissingValue() ) ;
			expected_result.back().push_back( -1000.556 ) ;
			expected_result.back().push_back( double( 10 ) ) ;
			
			CallReader( 2, "GT:F2", "1|1:0.1,-1.5568,-1e-08\t0/1:.,-1000.556,10", types )
			 	( "F2", boost::cref( GenotypeCallChecker( expected_result ) )) ;
		}
		catch( genfile::MalformedInputError const& ) {
			TEST_ASSERT(0) ;
		}

		try {
			std::vector< std::vector< Entry > > expected_result ;
			expected_result.push_back( std::vector< Entry >() ) ;
			expected_result.back().push_back( 0.1 ) ;
			expected_result.back().push_back( -1.5568 ) ;
			expected_result.back().push_back( -1e-08 ) ;
			expected_result.push_back( std::vector< Entry >() ) ;
			expected_result.back().push_back( genfile::MissingValue() ) ;
			expected_result.back().push_back( -1000.556 ) ;
			expected_result.back().push_back( double( 10 ) ) ;
			expected_result.push_back( std::vector< Entry >() ) ;
			expected_result.back().push_back( double( 1000 ) ) ;
			expected_result.back().push_back( double( 10000 ) ) ;
			expected_result.back().push_back( double( 100000 ) ) ;
			
			CallReader( 2, "GT:F2", "1|1:0.1,-1.5568,-1e-08\t0/1:.,-1000.556,10\t0|0:1000,10000,100000", types )
			 	( "F2", boost::cref( GenotypeCallChecker( expected_result ) )) ;
		}
		catch( genfile::MalformedInputError const& ) {
			TEST_ASSERT(0) ;
		}

		try {
			std::vector< std::vector< Entry > > expected_result1 ;
			expected_result1.push_back( std::vector< Entry >() ) ;
			expected_result1.back().push_back( 15.05 ) ;
			expected_result1.push_back( std::vector< Entry >() ) ;
			expected_result1.back().push_back( 110.110e02 ) ;
			expected_result1.push_back( std::vector< Entry >() ) ;
			expected_result1.back().push_back( 16.2 ) ;

			std::vector< std::vector< Entry > > expected_result2 ;
			expected_result2.push_back( std::vector< Entry >() ) ;
			expected_result2.back().push_back( 0.1 ) ;
			expected_result2.back().push_back( -1.5568 ) ;
			expected_result2.back().push_back( -1e-08 ) ;
			expected_result2.push_back( std::vector< Entry >() ) ;
			expected_result2.back().push_back( genfile::MissingValue() ) ;
			expected_result2.back().push_back( -1000.556 ) ;
			expected_result2.back().push_back( double( 10 ) ) ;
			expected_result2.push_back( std::vector< Entry >() ) ;
			expected_result2.back().push_back( double( 1000 ) ) ;
			expected_result2.back().push_back( double( 10000 ) ) ;
			expected_result2.back().push_back( double( 100000 ) ) ;
			
			CallReader( 2, "GT:F1:F2", "1|1:15.05:0.1,-1.5568,-1e-08\t0/1:110.110e02:.,-1000.556,10\t0|0:16.2:1000,10000,100000", types )
			 	( "F1", boost::cref( GenotypeCallChecker( expected_result1 ) ))
				( "F2", boost::cref( GenotypeCallChecker( expected_result2 ) )) ;

			CallReader( 2, "GT:F1:F2", "1|1:15.05:0.1,-1.5568,-1e-08\t0/1:110.110e02:.,-1000.556,10\t0|0:16.2:1000,10000,100000", types )
			 	( "F1", boost::cref( GenotypeCallChecker( expected_result1 ) ))
			 	( "F2", boost::cref( GenotypeCallChecker( expected_result2 ) ))
			 	( "F2", boost::cref( GenotypeCallChecker( expected_result2 ) )) ;

			CallReader( 2, "GT:F2:F1", "1|1:0.1,-1.5568,-1e-08:15.05\t0/1:.,-1000.556,10:110.110e02\t0|0:1000,10000,100000:16.2", types )
			 	( "F1", boost::cref( GenotypeCallChecker( expected_result1 ) ))
				( "F2", boost::cref( GenotypeCallChecker( expected_result2 ) )) ;
		}
		catch( genfile::MalformedInputError const& ) {
			TEST_ASSERT(0) ;
		}
	}

	std::cerr << "ok.\n" ;
}

AUTO_TEST_MAIN {
	test_format() ;
	test_alleles() ;
	test_omitted_values() ;
	test_setters() ;
	test_float_fields() ;
	test_simple_gt_values() ;
	test_gt_delimiter() ;
	test_complex_gt_values() ;
	test_gt_genotype_bounds() ;
}
