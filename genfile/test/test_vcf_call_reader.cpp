#include "genfile/vcf/CallReader.hpp"
#include "genfile/Error.hpp"
#include "test_case.hpp"

using namespace genfile::vcf ;

AUTO_TEST_CASE( test_format ) {
	std::cerr << "test_format..." ;

	boost::ptr_map< std::string, VCFEntryType > types ;

	// Format invalid, missing or bad GT field.
	try {
		CallReader( 0, "", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 0, "HQ", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}
	
	try {
		CallReader( 0, "HQ:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	// Format invalid, GT field not in types.
	try {
		CallReader( 0, "GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

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
		CallReader( 0, "GT", "", types ) ;
	}
	catch( genfile::BadArgumentError const& e ) {}

	try {
		CallReader( 0, "HQ:GT", "", types ) ;
		TEST_ASSERT( 0 ) ;
		
	}
	catch( genfile::BadArgumentError const& e ) {
		TEST_ASSERT( e.function().substr( 0, 26 ) == "genfile::vcf::CallReader::" ) ;
	}
}

AUTO_TEST_MAIN {
	test_format() ;
}