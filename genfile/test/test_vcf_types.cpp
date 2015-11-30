
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "test_case.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

using namespace genfile::vcf ;
using namespace std ;

BOOST_AUTO_TEST_SUITE( test_vcf_types )

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

namespace {
	void insert_elt(
		std::vector< Entry >* destination,
		Entry const& value
	) {
		destination->push_back( value ) ;
	}

	struct VariantEntriesEntriesSetter: public EntriesSetter {
		VariantEntriesEntriesSetter( std::vector< Entry >& entries ): m_entries( entries ), entry_i( 0 ) {} 
		void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
			m_entries.resize( n ) ;
		}
		virtual void set_value( genfile::MissingValue const value ) { m_entries[ entry_i++ ] = value ; }
		virtual void set_value( std::string& value ) { m_entries[ entry_i++ ] = value ; }
		virtual void set_value( Integer const value ) { m_entries[ entry_i++ ] = value ; }
		virtual void set_value( double const value ) { m_entries[ entry_i++ ] = value ; }
	private:
		std::vector< Entry >& m_entries ;
		std::size_t entry_i ;
	} ;

	template <typename Type >
	std::vector< Entry > parse(
		Type const& type,
		std::string const& value,
		std::size_t n_alleles,
		uint32_t ploidy
	) {
		std::vector< Entry > result ;
		VariantEntriesEntriesSetter setter( result ) ;
		type.parse( value, n_alleles, ploidy, setter ) ;
		return result ;
	}
	
	
}
AUTO_TEST_CASE( test_fixed_number_entry_type ) {
	std::cerr << "test_fixed_number_entry_type()..." ;

	using genfile::BadArgumentError ;
	using genfile::MissingValue ;
	
	// First test some simple parseable values.	
	{
		typedef std::vector< Entry > Result ;
		for( std::size_t ploidy = 0; ploidy < 10; ++ploidy ) {
			for( std::size_t n_alleles = 0; n_alleles < 10; ++n_alleles ) {
				Result r, rp ;
				{
					FixedNumberVCFEntryType type( 0, SimpleType::create( "Integer" ) ) ;
					r = parse( type, std::string( "" ), n_alleles, ploidy ) ;
					rp = parse( type, std::string( "" ), n_alleles, ploidy ) ;
					TEST_ASSERT( r == rp ) ;
					TEST_ASSERT( r.size() == 0 ) ;
				}
				{
					FixedNumberVCFEntryType type( 1, SimpleType::create( "Integer" ) ) ;
					r = parse( type, std::string( "1" ), n_alleles, ploidy ) ;
					rp = parse( type, std::string( "1" ), n_alleles, ploidy ) ;
					TEST_ASSERT( r == rp ) ;
					TEST_ASSERT( r.size() == 1 ) ;
					TEST_ASSERT( r[0] == Entry( 1 )) ;

					r = parse( type, std::string( "5" ), n_alleles, ploidy ) ;
					rp = parse( type, std::string( "5" ), n_alleles, ploidy ) ;
					TEST_ASSERT( r == rp ) ;
					TEST_ASSERT( r.size() == 1 ) ;
					TEST_ASSERT( r[0] == Entry( 5 )) ;

					r = parse( type, std::string( "-1000" ), n_alleles, ploidy ) ;
					rp = parse( type, std::string( "-1000" ), n_alleles, ploidy ) ;
					TEST_ASSERT( r == rp ) ;
					TEST_ASSERT( r.size() == 1 ) ;
					TEST_ASSERT( r[0] == Entry( -1000 )) ;
				}

				{
					FixedNumberVCFEntryType type( 5, SimpleType::create( "Integer" ) ) ;
					r = parse( type, std::string( "1,2,3,4,5" ), n_alleles, ploidy ) ;
					rp = parse( type, std::string( "1,2,3,4,5" ), n_alleles, ploidy ) ;
					TEST_ASSERT( r == rp ) ;
					TEST_ASSERT( r.size() == 5 ) ;
					TEST_ASSERT( r[0] == Entry( 1 )) ;
					TEST_ASSERT( r[1] == Entry( 2 )) ;
					TEST_ASSERT( r[2] == Entry( 3 )) ;
					TEST_ASSERT( r[3] == Entry( 4 )) ;
					TEST_ASSERT( r[4] == Entry( 5 )) ;

					r = parse( type, std::string( "5,4,.,0,-1" ), n_alleles, ploidy ) ;
					rp = parse( type, std::string( "5,4,.,0,-1" ), n_alleles, ploidy ) ;
					TEST_ASSERT( r == rp ) ;
					TEST_ASSERT( r.size() == 5 ) ;
					TEST_ASSERT( r[0] == Entry( 5 )) ;
					TEST_ASSERT( r[1] == Entry( 4 )) ;
					TEST_ASSERT( r[2] == Entry( MissingValue() )) ;
					TEST_ASSERT( r[3] == Entry( 0 )) ;
					TEST_ASSERT( r[4] == Entry( -1 )) ;

					r = parse( type, std::string( ".,.,.,.,." ), n_alleles, ploidy ) ;
					rp = parse( type, std::string( ".,.,.,.,." ), n_alleles, ploidy ) ;
					TEST_ASSERT( r == rp ) ;
					TEST_ASSERT( r.size() == 5 ) ;
					TEST_ASSERT( r[0] == Entry( MissingValue() )) ;
					TEST_ASSERT( r[1] == Entry( MissingValue() )) ;
					TEST_ASSERT( r[2] == Entry( MissingValue() )) ;
					TEST_ASSERT( r[3] == Entry( MissingValue() )) ;
					TEST_ASSERT( r[4] == Entry( MissingValue() )) ;
				}
			}
		}
	}

	uint32_t const ploidy = genfile::eUnknownPloidy ;
	
	// Now check lots of parseable and non-parseable values.
	for( std::size_t n_alleles = 0; n_alleles < 100; ++n_alleles ) {
		FixedNumberVCFEntryType type( 0, SimpleType::create( "Integer" ) ) ;
		FixedNumberVCFEntryType type2( 0, SimpleType::create( "String" ) ) ;
		try { parse(  type, std::string( "" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "-1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( ".,." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2,3" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}

		try { parse(  type, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
	}

	for( std::size_t n_alleles = 0; n_alleles < 100; ++n_alleles ) {
		FixedNumberVCFEntryType type( 1, SimpleType::create( "Integer" ) ) ;
		FixedNumberVCFEntryType type2( 1, SimpleType::create( "String" ) ) ;
		try { parse(  type, std::string( "" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "2" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "-1" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "." ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "1,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( ".,1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( ".,." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "-1,-1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2,3" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "-1,.,-2000" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}

		try { parse(  type, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "3.5" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
	}

	for( std::size_t n_alleles = 0; n_alleles < 100; ++n_alleles ) {
		FixedNumberVCFEntryType type( 2, SimpleType::create( "Integer" ) ) ;
		FixedNumberVCFEntryType type2( 2, SimpleType::create( "String" ) ) ;
		try { parse(  type, std::string( "" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "-1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "1,2" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( ".,." ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "-1,-1" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "1,2,3" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "-1,.,-2000" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}

		try { parse(  type, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type2, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
	}

	for( std::size_t n_alleles = 0; n_alleles < 100; ++n_alleles ) {
		FixedNumberVCFEntryType type( 3, SimpleType::create( "Integer" ) ) ;
		FixedNumberVCFEntryType type2( 3, SimpleType::create( "String" ) ) ;
		try { parse(  type, std::string( "" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "-1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "1,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( ".,." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ;  } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "-1,-1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2,3" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "1,.,2" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "-1,.,-2000" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "A,.,2" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type2, std::string( "1,2-3,2" ), n_alleles, ploidy ) ;  } catch( BadArgumentError const& ) { TEST_ASSERT(0) ;}
		try { parse(  type2, std::string( "1,2," ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }

		try { parse(  type, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2.," ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type2, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
	}

	for( std::size_t n_alleles = 0; n_alleles < 100; ++n_alleles ) {
		FixedNumberVCFEntryType type( 10, SimpleType::create( "Integer" ) ) ;
		FixedNumberVCFEntryType type2( 10, SimpleType::create( "String" ) ) ;
		try { parse(  type, std::string( "" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "-1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "1,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( ".,." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ;  } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "-1,-1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2,3" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "1,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {  }
		try { parse(  type, std::string( "-1,.,-2000" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "1,.,2,3,3,50,50,.,76,." ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}

		try { parse(  type, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type2, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
	}

	for( std::size_t n_alleles = 0; n_alleles < 100; ++n_alleles ) {
		FixedNumberVCFEntryType type( 300, SimpleType::create( "Integer" ) ) ;
		FixedNumberVCFEntryType type2( 300, SimpleType::create( "String" ) ) ;
		
		try { parse(  type, std::string( "" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "-1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "1,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( ".,." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ;  } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "-1,-1" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2,3" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "1,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {  }
		try { parse(  type, std::string( "-1,.,-2000" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "1,.,2,3,3,50,50,.,76,." ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}

		try { parse(  type, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "1,2.," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
		try { parse(  type2, std::string( "3.5" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) { }
	}

	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_dynamic_number_entry_type ) {
	std::cerr << "test_dynamic_number_entry_type()..." ;

	using namespace genfile ;

	// test some simple parseable values
	{
		typedef std::vector< Entry > Result ;
		DynamicNumberVCFEntryType type( SimpleType::create( "Integer" ) ) ;
		
		for( std::size_t ploidy = 0; ploidy < 10; ++ploidy ) {
			for( std::size_t n_alleles = 0; n_alleles < 10; ++n_alleles ) {
				Result r, rp ;
				r = parse(  type, std::string( "" ), n_alleles, ploidy ) ;
				rp = parse( type, std::string( "" ), n_alleles, ploidy ) ;
				TEST_ASSERT( r == rp ) ;
				TEST_ASSERT( r.size() == 0 ) ;

				r = parse(  type, std::string( "1" ), n_alleles, ploidy ) ;
				rp = parse( type, std::string( "1" ), n_alleles, ploidy ) ;
				TEST_ASSERT( r == rp ) ;
				TEST_ASSERT( r.size() == 1 ) ;
				TEST_ASSERT( r[0] == Entry( 1 )) ;

				r = parse(  type, std::string( "5" ), n_alleles, ploidy ) ;
				rp = parse( type, std::string( "5" ), n_alleles, ploidy ) ;
				TEST_ASSERT( r == rp ) ;
				TEST_ASSERT( r.size() == 1 ) ;
				TEST_ASSERT( r[0] == Entry( 5 )) ;

				r = parse(  type, std::string( "-1000" ), n_alleles, ploidy ) ;
				rp = parse( type, std::string( "-1000" ), n_alleles, ploidy ) ;
				TEST_ASSERT( r == rp ) ;
				TEST_ASSERT( r.size() == 1 ) ;
				TEST_ASSERT( r[0] == Entry( -1000 )) ;

				r = parse(  type, std::string( "1,2,3,4,5" ), n_alleles, ploidy ) ;
				rp = parse( type, std::string( "1,2,3,4,5" ), n_alleles, ploidy ) ;
				TEST_ASSERT( r == rp ) ;
				TEST_ASSERT( r.size() == 5 ) ;
				TEST_ASSERT( r[0] == Entry( 1 )) ;
				TEST_ASSERT( r[1] == Entry( 2 )) ;
				TEST_ASSERT( r[2] == Entry( 3 )) ;
				TEST_ASSERT( r[3] == Entry( 4 )) ;
				TEST_ASSERT( r[4] == Entry( 5 )) ;

				r = parse(  type, std::string( "5,4,.,0,-1" ), n_alleles, ploidy ) ;
				rp = parse( type, std::string( "5,4,.,0,-1" ), n_alleles, ploidy ) ;
				TEST_ASSERT( r == rp ) ;
				TEST_ASSERT( r.size() == 5 ) ;
				TEST_ASSERT( r[0] == Entry( 5 )) ;
				TEST_ASSERT( r[1] == Entry( 4 )) ;
				TEST_ASSERT( r[2] == Entry( MissingValue() )) ;
				TEST_ASSERT( r[3] == Entry( 0 )) ;
				TEST_ASSERT( r[4] == Entry( -1 )) ;

				r = parse(  type, std::string( ".,.,.,.,." ), n_alleles, ploidy ) ;
				rp = parse( type, std::string( ".,.,.,.,." ), n_alleles, ploidy ) ;
				TEST_ASSERT( r == rp ) ;
				TEST_ASSERT( r.size() == 5 ) ;
				TEST_ASSERT( r[0] == Entry( MissingValue() )) ;
				TEST_ASSERT( r[1] == Entry( MissingValue() )) ;
				TEST_ASSERT( r[2] == Entry( MissingValue() )) ;
				TEST_ASSERT( r[3] == Entry( MissingValue() )) ;
				TEST_ASSERT( r[4] == Entry( MissingValue() )) ;
			}
		}
	}

	uint32_t const ploidy = genfile::eUnknownPloidy ;
	// test the boundaries of parseability
	for( std::size_t n_alleles = 0; n_alleles < 100; ++n_alleles ) {
		DynamicNumberVCFEntryType type( SimpleType::create( "Integer" ) ) ;
		DynamicNumberVCFEntryType type2( SimpleType::create( "String" ) ) ;
		try { parse(  type, std::string( "" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "2" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "-1" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "." ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "1,2" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( ".,." ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "-1,-1" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "1,2,3" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "1,.,2" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "-1,.,-2000" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "1,.,2,3,3,50,50,.,76,." ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "A,.,2" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type2, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
		try { parse(  type2, std::string( "1,2," ), n_alleles, ploidy ) ; } catch( BadArgumentError const& ) { TEST_ASSERT(0) ; }
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_one_per_allele_entry_type ) {
	std::cerr << "test_one_per_alternate_allele_entry_type()..." ;

	using namespace genfile ;

	{
		typedef std::vector< Entry > Result ;
		OnePerAlleleVCFEntryType type( SimpleType::create( "Integer" ) ) ;
		Result r, rp ;
		for( std::size_t ploidy = 0; ploidy < 100; ++ploidy ) {
			r = parse(  type, std::string( "" ), 0, ploidy ) ;
			rp = parse( type, std::string( "" ), 0, ploidy ) ;
			TEST_ASSERT( r == rp ) ;
			TEST_ASSERT( r.size() == 0 ) ;

			r = parse(  type, std::string( "5,6" ), 2, ploidy ) ;
			rp = parse( type, std::string( "5,6" ), 2, ploidy ) ;
			TEST_ASSERT( r == rp ) ;
			TEST_ASSERT( r.size() == 2 ) ;
			TEST_ASSERT( r[0] == Entry( 5 )) ;
			TEST_ASSERT( r[1] == Entry( 6 )) ;

			r = parse(  type, std::string( "-1,.,5,." ), 4, ploidy ) ;
			rp = parse( type, std::string( "-1,.,5,." ), 4, ploidy ) ;
			TEST_ASSERT( r == rp ) ;
			TEST_ASSERT( r.size() == 4 ) ;
			TEST_ASSERT( r[0] == Entry( -1 )) ;
			TEST_ASSERT( r[1] == Entry( MissingValue() )) ;
			TEST_ASSERT( r[2] == Entry( 5 )) ;
			TEST_ASSERT( r[3] == Entry( MissingValue() )) ;
		}
	}
	
	uint32_t const ploidy = genfile::eUnknownPloidy ;
	for( std::size_t n_alleles = 0; n_alleles < 100; ++n_alleles ) {
		OnePerAlleleVCFEntryType type( SimpleType::create( "Integer" ) ) ;
		OnePerAlleleVCFEntryType type2( SimpleType::create( "String" ) ) ;
		try { parse(  type, std::string( "" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 0 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 0 ) ; }
		try { parse(  type, std::string( "," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 1 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 1 ) ; }
		try { parse(  type, std::string( "2" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 1 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 1 ) ; }
		try { parse(  type, std::string( "-1" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 1 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 1 ) ; }
		try { parse(  type, std::string( "." ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 1 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 1 ) ; }
		try { parse(  type, std::string( "1,2" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 2 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 2 ) ; }
		try { parse(  type, std::string( ".,." ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 2 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 2 ) ; }
		try { parse(  type, std::string( "-1,-1" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 2 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 2 ) ; }
		try { parse(  type, std::string( "1,2,3" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 3 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 3 ) ; }
		try { parse(  type, std::string( "1,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 3 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 3 ) ; }
		try { parse(  type, std::string( "-1,.,-2000" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 3 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 3 ) ; }
		try { parse(  type, std::string( "1,.,2,3,3,50,50,.,76,." ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 10 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 10 ) ; }
		try { parse(  type, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT(0) ; } catch( BadArgumentError const& ) {}
		try { parse(  type2, std::string( "A,.,2" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 3 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 3 ) ; }
		try { parse(  type2, std::string( "1,2-3,2" ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 3 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 3 ) ; }
		try { parse(  type2, std::string( "1,2," ), n_alleles, ploidy ) ; TEST_ASSERT( n_alleles == 3 ) ; } catch( BadArgumentError const& ) { TEST_ASSERT( n_alleles != 3 ) ; }
	}
	std::cerr << "ok.\n" ;
}

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
}

AUTO_TEST_CASE( test_one_per_genotype_entry_type ) {
	std::cerr << "test_one_per_genotype_entry_type()..." ;

	using genfile::BadArgumentError ;
	using genfile::string_utils::to_string ;

	typedef std::vector< Entry > Result ;
	OnePerGenotypeVCFEntryType type( SimpleType::create( "Integer" ) ) ;
	{
		Result r ;

		r = parse( type, std::string( "" ), 1, 0 ) ;
		TEST_ASSERT( r.size() == 0 ) ;

		// if number of alleles is zero, so should be the ploidy.
		try {
			r = parse( type, std::string( "" ), 0, 0 ) ;
			TEST_ASSERT( r.size() == 0 ) ;
		}
		catch( BadArgumentError const& ) {
			TEST_ASSERT(0) ;
		}

		try {
			r = parse( type, std::string( "" ), 0, 1 ) ;
			TEST_ASSERT( 0 ) ;
		}
		catch( BadArgumentError const& ) {}

		try {
			r = parse( type, std::string( "" ), 0, 2 ) ;
			TEST_ASSERT( 0 ) ;
		}
		catch( BadArgumentError const& ) {}

		// Test one allele with various ploidies.
		// (There is only one genotype.)
		for( std::size_t ploidy = 0; ploidy < 10; ++ploidy ) {
			try {
				r = parse( type, std::string( "0" ), 1, ploidy ) ;
				TEST_ASSERT( ploidy > 0 ) ;
				TEST_ASSERT( r.size() == 1 ) ;
				TEST_ASSERT( r[0] == Entry( 0 )) ;
			} catch( BadArgumentError const& ) {
				TEST_ASSERT( ploidy == 0 ) ;
			}
		}

		// Test two alleles, ploidy one, 2 genotypes.
		try {
			r = parse( type, std::string( "0,1" ), 2, 1 ) ;
			TEST_ASSERT( r.size() == 2 ) ;
			TEST_ASSERT( r[0] == Entry( 0 )) ;
			TEST_ASSERT( r[1] == Entry( 1 )) ;
		}
		catch( BadArgumentError const& ) {
			TEST_ASSERT(0) ;
		}
	}
	
	for( std::size_t ploidy = 0; ploidy < 4; ++ploidy ) {
		for( std::size_t number_of_alleles = 0; number_of_alleles < 10; ++number_of_alleles ) {
			std::size_t N = 0 ;
			// number of entries is the number of ways to fill an n-vector
			// with nonnegative integers adding up to the ploidy.
			if( ploidy == 1 ) {
				N = number_of_alleles ;
			}
			else if( ploidy == 2 ) {
				// one 2, or two ones.
				N = number_of_alleles ;
				N += ( number_of_alleles < 2 ) ? 0 : impl::n_choose_k( number_of_alleles, 2 ) ;
			}
			else if( ploidy == 3 ) {
				// one 3, one 1 and one 2, or three ones.
				N = number_of_alleles ;
				N += ( number_of_alleles < 2 ) ? 0 : ( number_of_alleles * ( number_of_alleles - 1 ) ) ;
				N += ( number_of_alleles < 3 ) ? 0 : impl::n_choose_k( number_of_alleles, 3 ) ;
			}

			for( std::size_t i = 0; i < 20; ++i ) {
				std::string data ;
				for( std::size_t j = 0; j < i; ++j ) {
					if( j > 0 ) {
						data += "," ;
					}
					data += genfile::string_utils::to_string( int(j * i) - 1 ) ;
				}
				try {
					// std::cerr << "i = " << i << ", number_of_alleles = " << number_of_alleles << ", ploidy = " << ploidy << ", N = " << N << ".\n" ;
					Result r = parse( type, data, number_of_alleles, ploidy ) ;

					TEST_ASSERT( ( number_of_alleles != 0 || ploidy == 0 ) && i == N ) ;

					TEST_ASSERT( r.size() == N ) ;
					for( std::size_t j = 0; j < N; ++j ) {
						TEST_ASSERT( r[j] == Entry( int(j*i) - 1 ) ) ;
					}
				}
				catch( BadArgumentError const& ) {
					TEST_ASSERT( ( number_of_alleles == 0 && ploidy > 0 ) || i != N ) ;
				}
			}
		}
	}

	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_genotype_call_entry_type ) {
	std::cerr << "test_genotype_call_entry_type()..." ;
	using namespace genfile ;
	using genfile::string_utils::to_string ;

	GenotypeCallVCFEntryType type ;
	typedef std::vector< Entry > Result ;

	uint32_t const ploidy = genfile::eUnknownPloidy ;

	for( std::size_t number_of_alleles = 0; number_of_alleles < 20; ++number_of_alleles ) {
		for( int c = 0; c < std::numeric_limits< char >::max(); ++c ) {
			while( c >= '0' && c <= '9' ) {
				++c ;
			}
			for( std::size_t n = 0; n < 20; ++n ) {
				std::string data ;
				int max_call = -1, min_call = 100000 ;
				for( std::size_t i = 0; i < n; ++i ) {
					if( i > 0 ) {
						data += std::string( 1, char(c) ) ;
					}
					int v = int(i+n) - 5 ;
					data += to_string( v ) ;
					max_call = std::max( v, max_call ) ;
					min_call = std::min( v, min_call ) ;
				}
				
				bool bad_call = ( max_call >= int( number_of_alleles ) ) || ( min_call < 0 ) ;
				bool bad_char = n > 1 && c != '|' && c != '/' ;
				
				try {
					// std::cerr << "number_of_alleles = " << number_of_alleles << ", c = " << c << "('" << char(c) << "'), n = " << n << ", data = \"" << data << "\".\n" ;
					// std::cerr << "bad_char = " << ( bad_char ? "true" : "false" ) << ", bad_call = " << ( bad_call ? "true" : "false" ) << ".\n" ;
					Result r = parse(  type, data, number_of_alleles, ploidy ) ;
					TEST_ASSERT( (!bad_char) && (!bad_call) ) ;
					TEST_ASSERT( r.size() == n ) ;
					for( std::size_t i = 0; i < n; ++i ) {
						TEST_ASSERT( r[i].as<int>() == ( int( i+n ) - 5 ) ) ;
					}
				}
				catch( BadArgumentError const& ) {
					TEST_ASSERT( bad_char || bad_call ) ;
				}
				
				if( (!bad_char) && (!bad_call) ) {
					for( std::size_t ploidy = 0; ploidy < 10; ++ploidy ) {
						// std::cerr << "ploidy = " << ploidy << ".\n" ;
						try {
							Result r = parse( type, data, number_of_alleles, ploidy ) ;
							TEST_ASSERT( ploidy == n ) ;
							TEST_ASSERT( r.size() == n ) ;
							for( std::size_t i = 0; i < n; ++i ) {
								TEST_ASSERT( r[i].as<int>() == int( ( i+n ) - 5 ) ) ;
							}
						}
						catch( BadArgumentError const& ) {
							TEST_ASSERT( ploidy != n ) ;
						}
					}
				}
			}
		}
	}
	
	std::cerr << "ok.\n" ;
}

void test_gt_spec() {

	std::string const GT = "GT" ;

	try {
		std::string const GT = "GT" ;
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "Integer" ;
		spec[ "Description" ] = "A test" ;
		VCFEntryType::create( spec ) ;
	}
	catch( genfile::BadArgumentError const& ) {
		TEST_ASSERT(0) ;
	}

	try {
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "." ;
		spec[ "Type" ] = "Integer" ;
		spec[ "Description" ] = "A test" ;
		VCFEntryType::create( spec ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {}

	try {
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "2" ;
		spec[ "Type" ] = "Integer" ;
		spec[ "Description" ] = "A test" ;
		VCFEntryType::create( spec ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {}

	try {
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "String" ;
		spec[ "Description" ] = "A test" ;
		VCFEntryType::create( spec ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {}

	try {
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "Float" ;
		spec[ "Description" ] = "A test" ;
		VCFEntryType::create( spec ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {}

	try {
		VCFEntryType::Spec spec ;
		spec[ "ID" ] = GT ;
		spec[ "Number" ] = "1" ;
		spec[ "Type" ] = "Flag" ;
		spec[ "Description" ] = "A test" ;
		VCFEntryType::create( spec ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::BadArgumentError const& ) {}	
}

BOOST_AUTO_TEST_SUITE_END()
