#include <sstream>
#include <vector>
#include <string>
#include <set>
#include "genfile/vcf/MetadataParser.hpp"
#include "genfile/Error.hpp"
#include "test_case.hpp"

namespace data {
	// This example is from http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
	std::string const ex1 =
		"##fileformat=VCFv4.1\n"
		"##fileDate=20090805\n"
		"##source=myImputationProgramV3.1\n"
		"##reference=file:/seq/references/1000GenomesPilot-NCBI36.fasta\n"
		"##contig=<name=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=Homo sapiens,taxonomy=x>\n"
		"##phasing=partial\n" ;

	std::string const ex1b =
		"##fileformat=VCFv4.0\n" ;

	std::string const ex1c =
		"##fileformat=VCFv3.9\n" ;
	
	std::string const ex2 =
		"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
		"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
		"##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"
		"##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">\n"
		"##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">\n" ;

	std::string const ex2b =
		"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
		"##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n"
		"##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"
		"##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">\n"
		"##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">\n" ;
	
	std::string const ex3 =
		"##FILTER=<ID=q10,Description=\"Quality below 10\">\n"
		"##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">\n" ;
		
	std::string const ex4 =
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n"
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
		"##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">\n" ;
		
	std::string const ex5 =
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003\n" ;

	std::string const ex6 = 
		"20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.\n"
		"20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3\n"
		"20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4\n"
		"20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2\n"
		"20\t1234567\tmicrosat1\tGTC\tG,GTCTC\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3\n" ;
}

AUTO_TEST_CASE( test_spec_example ) {
	using namespace data ;
	std::cerr << "test_spec_example()..." ;
	std::istringstream istr( ex1 + ex2 + ex3 + ex4 + ex5 + ex6 ) ;
	genfile::vcf::MetadataParser parser( "VCF_4_1_example_file", istr ) ;
	std::string s ;
	istr >> s ;
	TEST_ASSERT( s == "#CHROM" ) ;

	// Test the metadata is correct.
	// Actually, these tests might fail because multimap is not guaranteed to preserve
	// the order of the elements inseretd in a range of elements comparing equal.
	// For now I'll just assume it does though.
	genfile::vcf::MetadataParser::Metadata metadata = parser.get_metadata() ;
	TEST_ASSERT( std::distance( metadata.equal_range( "fileformat" ).first, metadata.equal_range( "fileformat" ).second ) == 1 )  ;
	TEST_ASSERT( std::distance( metadata.equal_range( "fileDate" ).first, metadata.equal_range( "fileDate" ).second ) == 1 )  ;
	TEST_ASSERT( std::distance( metadata.equal_range( "source" ).first, metadata.equal_range( "source" ).second ) == 1 )  ;
	TEST_ASSERT( std::distance( metadata.equal_range( "reference" ).first, metadata.equal_range( "reference" ).second ) == 1 )  ;
	TEST_ASSERT( std::distance( metadata.equal_range( "contig" ).first, metadata.equal_range( "contig" ).second ) == 1 )  ;
	TEST_ASSERT( std::distance( metadata.equal_range( "phasing" ).first, metadata.equal_range( "phasing" ).second ) == 1 )  ;
	TEST_ASSERT( std::distance( metadata.equal_range( "INFO" ).first, metadata.equal_range( "INFO" ).second ) == 6 )  ;
	TEST_ASSERT( std::distance( metadata.equal_range( "FILTER" ).first, metadata.equal_range( "FILTER" ).second ) == 2 )  ;
	TEST_ASSERT( std::distance( metadata.equal_range( "FORMAT" ).first, metadata.equal_range( "FORMAT" ).second ) == 4 )  ;
	TEST_ASSERT( metadata.equal_range( "DUMMY" ).first == metadata.equal_range( "DUMMY" ).second ) ;
	TEST_ASSERT( metadata.equal_range( "INFO2" ).first == metadata.equal_range( "INFO2" ).second ) ;

	// Check data for some rows
	genfile::vcf::MetadataParser::Metadata::const_iterator where ;
	where = metadata.equal_range( "fileformat" ).first ;
	TEST_ASSERT( where->second.find( "version" ) != where->second.end() ) ;
	TEST_ASSERT( where->second.find( "version" )->second == "4.1" ) ;

	// INFO
	where = metadata.equal_range( "INFO" ).first ;
	TEST_ASSERT( where->second.find( "ID" ) != where->second.end() ) ;
	TEST_ASSERT( where->second.find( "ID" )->second == "NS" ) ;
	TEST_ASSERT( where->second.find( "Number" ) != where->second.end() ) ;
	TEST_ASSERT( where->second.find( "Number" )->second == "1" ) ;
	TEST_ASSERT( where->second.find( "Type" ) != where->second.end() ) ;
	TEST_ASSERT( where->second.find( "Type" )->second == "Integer" ) ;
	TEST_ASSERT( where->second.find( "Description" ) != where->second.end() ) ;
	TEST_ASSERT( where->second.find( "Description" )->second == "\"Number of Samples With Data\"" ) ;
	++++where ;
	TEST_ASSERT( where->second.find( "ID" ) != where->second.end() ) ;
	TEST_ASSERT( where->second.find( "ID" )->second == "AF" ) ;
	TEST_ASSERT( where->second.find( "Number" ) != where->second.end() ) ;
	TEST_ASSERT( where->second.find( "Number" )->second == "A" ) ;
	TEST_ASSERT( where->second.find( "Type" ) != where->second.end() ) ;
	TEST_ASSERT( where->second.find( "Type" )->second == "Float" ) ;
	TEST_ASSERT( where->second.find( "Description" ) != where->second.end() ) ;
	TEST_ASSERT( where->second.find( "Description" )->second == "\"Allele Frequency\"" ) ;
	
	// FILTER
	where = metadata.equal_range( "FILTER" ).first ;
	
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_no_fileformat ) {
	using namespace data ;
	std::cerr << "test_no_fileformat()..." ;
	std::string content = ex1 + ex2 + ex3 + ex4 + ex5 + ex6 ;
	// turn it into "gileformat"
	content[3] = 'g' ;
	std::istringstream istr( content ) ;
	try {
		genfile::vcf::MetadataParser( "VCF_4_1_example_file", istr ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::MalformedInputError const& e ) {
		// ok.
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_malformed_file ) {
	using namespace data ;
	std::cerr << "test_malformed_file()..." ;
	// ex1 must be present.  Metadata ends with EOF or a a final line not starting with ##.
	try {
		std::istringstream istr( ex1 + ex2 + ex3 + ex4 ) ;
		genfile::vcf::MetadataParser( "VCF_4_1_example_file", istr ) ;
	}
	catch( genfile::MalformedInputError const& e ) {
		TEST_ASSERT(0) ;
	}

	try {
		std::istringstream istr( ex2 + ex3 + ex4 ) ;
		genfile::vcf::MetadataParser( "VCF_4_1_example_file", istr ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::MalformedInputError const& e ) {
		// ok.
	}

	try {
		std::istringstream istr( ex2 + ex3 + ex4 + ex5 ) ;
		genfile::vcf::MetadataParser( "VCF_4_1_example_file", istr ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::MalformedInputError const& e ) {
		// ok.
	}
	std::cerr << "ok.\n" ;
}	

void test_info_or_format_field( std::string const& field ) {
	TEST_ASSERT( field == "INFO" || field == "FORMAT" ) ;
	std::cerr << "test_info_or_format_field(" << field << ")..." ;
	std::vector< std::string > types ;
	types.push_back( "Integer" ) ;
	types.push_back( "Float" ) ;
	types.push_back( "String" ) ;
	types.push_back( "Character" ) ;
	types.push_back( "Flag" ) ;
	types.push_back( "." ) ;
	types.push_back( "N" ) ;

	std::set< std::size_t > bad_types ;
	bad_types.insert( 5 ) ;
	bad_types.insert( 6 ) ;
	if( field == "FORMAT" ) {
		bad_types.insert( 4 ) ;
	}

	std::vector< std::string > numbers ;
	numbers.push_back( "0" ) ;
	numbers.push_back( "1" ) ;
	numbers.push_back( "100" ) ;
	numbers.push_back( "A" ) ;
	numbers.push_back( "G" ) ;
	numbers.push_back( "X" ) ;
	numbers.push_back( "-1" ) ;
	numbers.push_back( "-1000" ) ;
	
	std::set< std::size_t > bad_numbers ;
	bad_numbers.insert( 5 ) ;
	bad_numbers.insert( 6 ) ;
	bad_numbers.insert( 7 ) ;
	
	std::vector< std::string > descriptions ;
	descriptions.push_back( "" ) ;
	descriptions.push_back( "Hi there" ) ;
	descriptions.push_back( "\"\"" ) ;
	descriptions.push_back( "\"Hi there\"" ) ;
	std::set< std::size_t > bad_descriptions ;
	bad_descriptions.insert( 0 ) ;
	bad_descriptions.insert( 1 ) ;

	std::vector< std::string > IDs ;
	IDs.push_back( "AnID" ) ;
	IDs.push_back( "HG" ) ;
	std::set< std::size_t > bad_IDs ;
	
	using namespace data ;

	for( std::size_t id_i = 0; id_i < IDs.size(); ++id_i ) {
		std::string id = IDs[id_i] ;
		for( std::size_t number_i = 0; number_i < numbers.size(); ++number_i ) {
			std::string number = numbers[ number_i ] ;
			for( std::size_t type_i = 0; type_i < types.size(); ++type_i ) {
				std::string type = types[ type_i ] ;
				for( std::size_t description_i = 0; description_i < descriptions.size(); ++description_i ) {
					std::string description = descriptions[ description_i ] ;
					std::string content = ex1 ;
					content = content + "##" + field + "=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=" + description + ">\n" ;
					content = content + "##" + field + "=<Number=" + number + ",ID=" + id + ",Type=" + type + ",Description=" + description + ">\n" ;
					content = content + "END" ;
					std::istringstream istr( content ) ;
					try {
						genfile::vcf::MetadataParser( "VCF_4_1_example_file", istr ) ;
						TEST_ASSERT(
							bad_IDs.count( id_i ) == 0
							&& bad_numbers.count( number_i ) == 0
							&& bad_types.count( type_i ) == 0
							&& bad_descriptions.count( description_i ) == 0
						) ;
						std::string s ;
						istr >> s ;
						TEST_ASSERT( s == "END" ) ;
					}
					catch( genfile::MalformedInputError const& ) {
						TEST_ASSERT(
							bad_IDs.count( id_i ) > 0
							|| bad_numbers.count( number_i ) > 0
							|| bad_types.count( type_i ) > 0
							|| bad_descriptions.count( description_i ) > 0
						) ;
					}
				}
			}
		}
	}

	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_info_and_format_fields ) {
	test_info_or_format_field( "INFO" ) ;
	test_info_or_format_field( "FORMAT" ) ;
}

AUTO_TEST_CASE( test_filter_fields ) {
	std::cerr << "test_filter_fields()..." ;
	std::vector< std::string > IDs ;
	IDs.push_back( "AnID" ) ;
	IDs.push_back( "HG" ) ;
	IDs.push_back( "222 222" ) ;
	std::set< std::size_t > bad_IDs ;
	
	std::vector< std::string > descriptions ;
	descriptions.push_back( "" ) ;
	descriptions.push_back( "Hi there" ) ;
	descriptions.push_back( "\"\"" ) ;
	descriptions.push_back( "\"Hi there\"" ) ;
	std::set< std::size_t > bad_descriptions ;
	bad_descriptions.insert( 0 ) ;
	bad_descriptions.insert( 1 ) ;

	using namespace data ;

	for( std::size_t id_i = 0; id_i < IDs.size(); ++id_i ) {
		std::string id = IDs[id_i] ;
		for( std::size_t description_i = 0; description_i < descriptions.size(); ++description_i ) {
			std::string description = descriptions[ description_i ] ;
			std::string content = ex1 ;
			content = content + "##FILTER=<ID=" + id + ",Description=" + description + ">\n" ;
			content = content + "##FILTER=<Description=" + description + ",ID=" + id + ">\n" ;
			content = content + "END" ;
			std::istringstream istr( content ) ;
			try {
				genfile::vcf::MetadataParser( "VCF_4_1_example_file", istr ) ;
				TEST_ASSERT(
					bad_IDs.count( id_i ) == 0
					&& bad_descriptions.count( description_i ) == 0
				) ;
				std::string s ;
				istr >> s ;
				TEST_ASSERT( s == "END" ) ;
			}
			catch( genfile::MalformedInputError const& ) {
				TEST_ASSERT(
					bad_IDs.count( id_i ) > 0
					|| bad_descriptions.count( description_i ) > 0
				) ;
			}
		}
	}

	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_tab_characters ) {
	std::cerr << "test_tab_characters()..." ;
	// no tab characters are allowed anywhere in the metadata.
	using namespace data ;
	std::string base_content = ex1 + ex2 + ex3 + ex4 ;
	for( std::size_t i = 0; i < base_content.size(); ++i ) {
		std::string content = base_content.substr(0, i) + "\t" + base_content.substr( i, base_content.size() ) + "END" ;
		// std::cerr << i << " " ;
		try {
			std::istringstream istr( content ) ;
			genfile::vcf::MetadataParser( "VCF_4_1_example_file", istr ) ;
			TEST_ASSERT( 0 ) ;
		}
		catch( genfile::MalformedInputError const& ) {
			// ok.
		}
	}
	std::cerr << "ok.\n" ;
}

AUTO_TEST_CASE( test_v4_0 ) {
	using namespace data ;
	std::cerr << "test_v4_0()..." ;
	try {
		std::istringstream istr( ex1b + ex2 + ex3 + ex4 + ex5 + ex6 ) ;
		genfile::vcf::MetadataParser parser( "VCF_4_0_example_file", istr ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::MalformedInputError const& ) {
		// ok.
	}

	try {
		std::istringstream istr( ex1b + ex2b + ex3 + ex4 + ex5 + ex6 ) ;
		genfile::vcf::MetadataParser parser( "VCF_4_0_example_file", istr ) ;
		// ok.
	}
	catch( genfile::MalformedInputError const& ) {
		TEST_ASSERT(0) ;
	}

	try {
		std::istringstream istr( ex1c + ex2 + ex3 + ex4 + ex5 + ex6 ) ;
		genfile::vcf::MetadataParser parser( "VCF_4_0_example_file", istr ) ;
		TEST_ASSERT(0) ;
	}
	catch( genfile::FormatUnsupportedError const& ) {
		// ok.
	}
	std::cerr << "ok.\n" ;
}

