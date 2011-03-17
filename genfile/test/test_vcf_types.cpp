#include "test_case.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/Error.hpp"

using namespace genfile::vcf ;
using namespace std ;

AUTO_TEST_CASE( test_string ) {
	std::cerr << "test_string()..." ;
	TEST_ASSERT( StringType().parse( string( "" ) ).as< string >() == "" ) ;
	TEST_ASSERT( StringType().parse( string( "Hello. there" ) ).as< string >() == "Hello. there" ) ;
	TEST_ASSERT( StringType().parse( string( "\t\n\rb" ) ).as< string >() == "\t\n\rb" ) ;
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_character ) {
	std::cerr << "test_character()..." ;
	try {
		CharacterType().parse( string( "" ) ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}

	try {
		CharacterType().parse( string( "Hello" ) ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}

	for( int c = 0; c < 256; ++c ) {
		TEST_ASSERT( CharacterType().parse( string( 1, char(c) ) ).as< string >() == string( 1, char(c) ) ) ;
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_integer ) {
	std::cerr << "test_integer()..." ;
	TEST_ASSERT( IntegerType().parse( string( "0" )).as< int >() == 0 ) ;
	TEST_ASSERT( IntegerType().parse( string( "-0" )).as< int >() == 0 ) ;
	TEST_ASSERT( IntegerType().parse( string( "1000" )).as< int >() == 1000 ) ;
	TEST_ASSERT( IntegerType().parse( string( "-10000" )).as< int >() == -10000 ) ;
	
	try {
		IntegerType().parse( string( "" )) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}

	try {
		IntegerType().parse( string( "f" )) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}

	try {
		IntegerType().parse( string( "5 2" )) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_float ) {
	std::cerr << "test_float()..." ;
	TEST_ASSERT( FloatType().parse( string( "0" )).as< double >() == 0.0 ) ;
	TEST_ASSERT( FloatType().parse( string( "-0" )).as< double >() == -0.0 ) ;
	TEST_ASSERT( FloatType().parse( string( "0.0" )).as< double >() == 0.0 ) ;
	TEST_ASSERT( FloatType().parse( string( "-0.0" )).as< double >() == -0.0 ) ;
	TEST_ASSERT( FloatType().parse( string( "1000" )).as< double >() == 1000 ) ;
	TEST_ASSERT( FloatType().parse( string( "-10000" )).as< double >() == -10000 ) ;
	TEST_ASSERT( FloatType().parse( string( "5E02" )).as< double >() == 500.0 ) ;
	TEST_ASSERT( FloatType().parse( string( "1.1E03" )).as< double >() == 1100.0 ) ;
	TEST_ASSERT( FloatType().parse( string( "1.2345E03" )).as< double >() == 1234.5 ) ;
	TEST_ASSERT( FloatType().parse( string( "-1.2345E03" )).as< double >() == -1234.5 ) ;
	
	try {
		FloatType().parse( string( "" )) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}

	try {
		FloatType().parse( string( "f" )) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}

	try {
		FloatType().parse( string( "5 2" )) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}

	try {
		FloatType().parse( string( " -1.2345E03" )) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}

	try {
		FloatType().parse( string( "\n-1.2345E03" )) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}

	try {
		FloatType().parse( string( "-1.2345E03g" )) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {
		// ok.
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_MAIN {
	test_string() ;
	test_character() ;
	test_integer() ;
	test_float() ;
}