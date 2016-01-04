
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
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 0 ) ;
		BOOST_CHECK_EQUAL( data.number_of_identifiers(), 1 ) ;
		BOOST_CHECK_EQUAL( data.get_rsid(), "(unknown variant)" ) ;
	}

	{
		VariantIdentifyingData data( "rs1", GenomePosition( Chromosome( "01" ), 1000 ), "A", "G" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 2 ) ;
		BOOST_CHECK_EQUAL( data.number_of_identifiers(), 1 ) ;
		BOOST_CHECK_EQUAL( data.get_rsid(), "rs1" ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers_as_string( "," ), "rs1" ) ;
	}

	{
		VariantIdentifyingData data( "SNP1", "rs1", GenomePosition( Chromosome( "01" ), 1000 ), "A", "G" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 2 ) ;
		BOOST_CHECK_EQUAL( data.number_of_identifiers(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_rsid(), "rs1" ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers_as_string( "," ), "rs1,SNP1" ) ;
	}
}

BOOST_AUTO_TEST_CASE( test_identifiers ) {
	using namespace genfile ;
	{
		VariantIdentifyingData data( "rs0" ) ;
		data.set_primary_id( "rs1" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 0 ) ;
		BOOST_CHECK_EQUAL( data.number_of_identifiers(), 1 ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers().size(), 1 ) ;
		BOOST_CHECK_EQUAL( data.get_rsid(), "rs1" ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers_as_string( "," ), "rs1" ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers_as_string( ",", 1 ), "" ) ;

		data.add_identifier( "rs1" ) ;
		BOOST_CHECK_EQUAL( data.number_of_identifiers(), 1 ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers().size(), 1 ) ;
		BOOST_CHECK_EQUAL( data.get_rsid(), "rs1" ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers_as_string( "," ), "rs1" ) ;

		data.add_identifier( "rsx" ) ;
		BOOST_CHECK_EQUAL( data.number_of_identifiers(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers().size(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_rsid(), "rs1" ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers_as_string( "," ), "rs1,rsx" ) ;

		data.add_identifier( "rsx" ) ;
		BOOST_CHECK_EQUAL( data.number_of_identifiers(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers().size(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_rsid(), "rs1" ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers_as_string( "," ), "rs1,rsx" ) ;

		data.add_identifier( "rs1" ) ;
		BOOST_CHECK_EQUAL( data.number_of_identifiers(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers().size(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_rsid(), "rs1" ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers_as_string( "," ), "rs1,rsx" ) ;

		data.add_identifier( "rsy" ) ;
		BOOST_CHECK_EQUAL( data.number_of_identifiers(), 3 ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers().size(), 3 ) ;
		BOOST_CHECK_EQUAL( data.get_rsid(), "rs1" ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers_as_string( "," ), "rs1,rsx,rsy" ) ;
		BOOST_CHECK_EQUAL( data.get_identifiers_as_string( ",", 1 ), "rsx,rsy" ) ;
	}
	
	{
		genfile::VariantIdentifyingData snp( "RSID_1", genfile::GenomePosition( genfile::Chromosome(), 0 ), "A", "G" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 1 ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 1 ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 2 ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "SNPID_1" ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 2 ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "SNPID_1" ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 2 ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "SNPID_1" ) ;
		snp.add_identifier( std::string( "AnotherID" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 3 ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[2], "AnotherID" ) ;
		snp.add_identifier( std::string( "RSID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 3 ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[2], "AnotherID" ) ;
		snp.add_identifier( std::string( "SNPID_1" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 3 ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "RSID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "SNPID_1" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[2], "AnotherID" ) ;
	}
	
	{
		genfile::VariantIdentifyingData snp( "rs1" ) ;
		snp.set_primary_id( std::string( "rs1_changed" ) ) ;
		BOOST_CHECK_EQUAL( snp.get_rsid(), "rs1_changed" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 1 ) ;
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
		BOOST_CHECK_EQUAL( snp.get_identifiers().size(), 2 ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[0], "rs1_changed" ) ;
		BOOST_CHECK_EQUAL( snp.get_identifiers()[1], "ID" ) ;
	}
}

BOOST_AUTO_TEST_CASE( test_alleles ) {
	using namespace genfile ;
	using genfile::string_utils::slice ;
	{
		VariantIdentifyingData data( "rs0" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 0 ) ;
		data.add_allele( "A" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 1 ) ;
		BOOST_CHECK_EQUAL( data.get_alleles().size(), 1 ) ;
		BOOST_CHECK_EQUAL( data.get_allele(0), "A" ) ;
		data.add_allele( "A" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_alleles().size(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_allele(0), "A" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(1), "A" ) ;
		data.add_allele( "G" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 3 ) ;
		BOOST_CHECK_EQUAL( data.get_alleles().size(), 3 ) ;
		BOOST_CHECK_EQUAL( data.get_allele(0), "A" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(1), "A" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(2), "G" ) ;
		data.add_allele( "CATTAC" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 4 ) ;
		BOOST_CHECK_EQUAL( data.get_alleles().size(), 4 ) ;
		BOOST_CHECK_EQUAL( data.get_allele(0), "A" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(1), "A" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(2), "G" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(3), "CATTAC" ) ;
	}

	{
		VariantIdentifyingData data( "rs0" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 0 ) ;
		BOOST_CHECK_EQUAL( data.get_alleles().size(), 0 ) ;
		data.add_allele( "A" ) ;
		data.set_allele( 0, "G" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 1 ) ;
		BOOST_CHECK_EQUAL( data.get_alleles().size(), 1 ) ;
		BOOST_CHECK_EQUAL( data.get_allele(0), "G" ) ;
		data.add_allele( "A" ) ;
		data.set_allele( 0, "C" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_alleles().size(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_allele(0), "C" ) ;
	}

	{
		VariantIdentifyingData data( "rs0", genfile::GenomePosition( genfile::Chromosome( "01" ), 1000 ), "A", "G" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 2 ) ;
		data.set_allele( 0, "CC" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_allele(0), "CC" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(1), "G" ) ;
		data.set_allele( 1, "GG" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 2 ) ;
		BOOST_CHECK_EQUAL( data.get_allele(0), "CC" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(1), "GG" ) ;
		data.add_allele( "TT" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 3 ) ;
		BOOST_CHECK_EQUAL( data.get_allele(0), "CC" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(1), "GG" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(2), "TT" ) ;
		data.set_allele( 0, "TTT" ) ;
		BOOST_CHECK_EQUAL( data.number_of_alleles(), 3 ) ;
		BOOST_CHECK_EQUAL( data.get_allele(0), "TTT" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(1), "GG" ) ;
		BOOST_CHECK_EQUAL( data.get_allele(2), "TT" ) ;
	}
}


BOOST_AUTO_TEST_SUITE_END()


