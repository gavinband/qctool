
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/GenomePosition.hpp"
#include "test_case.hpp"

BOOST_AUTO_TEST_SUITE( test_variant_identifying_data )

BOOST_AUTO_TEST_CASE( test_constructors ) {
	using namespace genfile ;
	{
		VariantIdentifyingData data ;
		TEST_ASSERT( data.number_of_alleles() == 0 ) ;
		TEST_ASSERT( data.number_of_identifiers() == 1 ) ;
		TEST_ASSERT( data.get_rsid() == "(unknown variant)" ) ;
	}

	{
		VariantIdentifyingData data( "rs1", GenomePosition( Chromosome( "01" ), 1000 ), "A", "G" ) ;
		TEST_ASSERT( data.number_of_alleles() == 2 ) ;
		TEST_ASSERT( data.number_of_identifiers() == 1 ) ;
		TEST_ASSERT( data.get_rsid() == "rs1" ) ;
		TEST_ASSERT( data.get_alternate_identifiers_as_string() == "" ) ;
	}

	{
		VariantIdentifyingData data( "SNP1", "rs1", GenomePosition( Chromosome( "01" ), 1000 ), "A", "G" ) ;
		TEST_ASSERT( data.number_of_alleles() == 2 ) ;
		TEST_ASSERT( data.number_of_identifiers() == 2 ) ;
		TEST_ASSERT( data.get_rsid() == "rs1" ) ;
		TEST_ASSERT( data.get_alternate_identifiers_as_string() == "SNP1" ) ;
	}
}

BOOST_AUTO_TEST_CASE( test_identifiers ) {
	using namespace genfile ;
	{
		VariantIdentifyingData data( "rs0" ) ;
		data.set_primary_id( "rs1" ) ;
		TEST_ASSERT( data.number_of_alleles() == 0 ) ;
		TEST_ASSERT( data.number_of_identifiers() == 1 ) ;
		TEST_ASSERT( data.get_rsid() == "rs1" ) ;
		TEST_ASSERT( data.get_alternate_identifiers_as_string() == "" ) ;

		data.add_identifier( "rs1" ) ;
		TEST_ASSERT( data.number_of_identifiers() == 1 ) ;
		TEST_ASSERT( data.get_rsid() == "rs1" ) ;
		TEST_ASSERT( data.get_alternate_identifiers_as_string() == "" ) ;

		data.add_identifier( "rsx" ) ;
		TEST_ASSERT( data.number_of_identifiers() == 2 ) ;
		TEST_ASSERT( data.get_rsid() == "rs1" ) ;
		TEST_ASSERT( data.get_alternate_identifiers_as_string() == "rsx" ) ;

		data.add_identifier( "rsx" ) ;
		TEST_ASSERT( data.number_of_identifiers() == 2 ) ;
		TEST_ASSERT( data.get_rsid() == "rs1" ) ;
		TEST_ASSERT( data.get_alternate_identifiers_as_string() == "rsx" ) ;

		data.add_identifier( "rsy" ) ;
		TEST_ASSERT( data.number_of_identifiers() == 3 ) ;
		TEST_ASSERT( data.get_rsid() == "rs1" ) ;
		TEST_ASSERT( data.get_alternate_identifiers_as_string() == "rsx,rsy" ) ;
	}
	
	{
		genfile::VariantIdentifyingData snp( "rs1" ) ;
		snp.set_primary_id( std::string( "rs1_changed" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_rsid(), "rs1_changed" ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 0 ) ;
	}

	{
		genfile::GenomePosition pos( genfile::Chromosome( "01" ), 1000 ) ;
		genfile::VariantIdentifyingData snp( "rs1", pos, "A", "G" ) ;
		snp.add_identifier( std::string( "ID" ) ) ;
		snp.set_primary_id( std::string( "rs1_changed" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_rsid(), "rs1_changed" ) ;
		BOOST_CHECK_EQUAL( snp.get_allele(0), "A" ) ;
		BOOST_CHECK_EQUAL( snp.get_allele(1), "G" ) ;
		BOOST_CHECK_EQUAL( snp.get_position(), pos ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers().size(), 1 ) ;
		BOOST_CHECK_EQUAL( snp.get_alternative_identifiers()[0], "ID" ) ;
	}
}

BOOST_AUTO_TEST_CASE( test_alleles ) {
	using namespace genfile ;
	using genfile::string_utils::slice ;
	{
		VariantIdentifyingData data( "rs0" ) ;
		TEST_ASSERT( data.number_of_alleles() == 0 ) ;
		data.add_allele( "A" ) ;
		TEST_ASSERT( data.number_of_alleles() == 1 ) ;
		TEST_ASSERT( data.get_allele(0) == "A" ) ;
		data.add_allele( "A" ) ;
		TEST_ASSERT( data.number_of_alleles() == 2 ) ;
		TEST_ASSERT( data.get_allele(0) == "A" ) ;
		TEST_ASSERT( data.get_allele(1) == "A" ) ;
		data.add_allele( "G" ) ;
		TEST_ASSERT( data.number_of_alleles() == 3 ) ;
		TEST_ASSERT( data.get_allele(0) == "A" ) ;
		TEST_ASSERT( data.get_allele(1) == "A" ) ;
		TEST_ASSERT( data.get_allele(2) == "G" ) ;
		data.add_allele( "CATTAC" ) ;
		TEST_ASSERT( data.number_of_alleles() == 4 ) ;
		TEST_ASSERT( data.get_allele(0) == "A" ) ;
		TEST_ASSERT( data.get_allele(1) == "A" ) ;
		TEST_ASSERT( data.get_allele(2) == "G" ) ;
		TEST_ASSERT( data.get_allele(3) == "CATTAC" ) ;
	}

	{
		VariantIdentifyingData data( "rs0" ) ;
		TEST_ASSERT( data.number_of_alleles() == 0 ) ;
		data.add_allele( "A" ) ;
		data.set_allele( 0, "G" ) ;
		TEST_ASSERT( data.number_of_alleles() == 1 ) ;
		TEST_ASSERT( data.get_allele(0) == "G" ) ;
		data.add_allele( "A" ) ;
		data.set_allele( 0, "C" ) ;
		TEST_ASSERT( data.number_of_alleles() == 2 ) ;
		TEST_ASSERT( data.get_allele(0) == "C" ) ;
	}

	{
		VariantIdentifyingData data( "rs0", genfile::GenomePosition( genfile::Chromosome( "01" ), 1000 ), "A", "G" ) ;
		TEST_ASSERT( data.number_of_alleles() == 2 ) ;
		data.set_allele( 0, "CC" ) ;
		TEST_ASSERT( data.number_of_alleles() == 2 ) ;
		TEST_ASSERT( data.get_allele(0) == "CC" ) ;
		TEST_ASSERT( data.get_allele(1) == "G" ) ;
		data.set_allele( 1, "GG" ) ;
		TEST_ASSERT( data.number_of_alleles() == 2 ) ;
		TEST_ASSERT( data.get_allele(0) == "CC" ) ;
		TEST_ASSERT( data.get_allele(1) == "GG" ) ;
		data.add_allele( "TT" ) ;
		TEST_ASSERT( data.number_of_alleles() == 3 ) ;
		TEST_ASSERT( data.get_allele(0) == "CC" ) ;
		TEST_ASSERT( data.get_allele(1) == "GG" ) ;
		TEST_ASSERT( data.get_allele(2) == "TT" ) ;
		data.set_allele( 0, "TTT" ) ;
		TEST_ASSERT( data.number_of_alleles() == 3 ) ;
		TEST_ASSERT( data.get_allele(0) == "TTT" ) ;
		TEST_ASSERT( data.get_allele(1) == "GG" ) ;
		TEST_ASSERT( data.get_allele(2) == "TT" ) ;
	}
}


BOOST_AUTO_TEST_SUITE_END() ;


