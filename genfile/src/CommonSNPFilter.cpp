#include <cassert>
#include <iterator>

#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/CommonSNPFilter.hpp"
#include "genfile/RSIDInListTest.hpp"
#include "genfile/SNPIDInListTest.hpp"
#include "genfile/SNPIDFieldsInListTest.hpp"
#include "genfile/ChromosomeInSetTest.hpp"
#include "genfile/SNPIDMatchesTest.hpp"
#include "genfile/SNPPositionInRangeTest.hpp"
#include "genfile/snp_data_utils.hpp"

namespace genfile {
	CommonSNPFilter::CommonSNPFilter() {
	}

	bool CommonSNPFilter::operator()(
		std::string SNPID,
		std::string RSID,
		GenomePosition position,
		std::string first_allele,
		std::string second_allele
	) const {
		return m_filter( SNPID, RSID, position, first_allele, second_allele ) ;
	}
		
	std::string CommonSNPFilter::display() const { return m_filter.display() ; }
	
	CommonSNPFilter& CommonSNPFilter::exclude_snps_in_file( std::string const& filename, int fields ) {
		return exclude_snps_in_set( read_strings_from_file( filename ), fields ) ;
	}

	CommonSNPFilter& CommonSNPFilter::include_snps_in_file( std::string const& filename, int fields ) {
		return include_snps_in_set( read_strings_from_file( filename ), fields ) ;
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps_not_in_file( std::string const& filename, int fields ) {
		return exclude_snps_not_in_set( read_strings_from_file( filename ), fields ) ;
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps_in_set( std::set< std::string > const& set, int fields ) {
		SNPIdentifyingDataTest::UniquePtr test = construct_snp_exclusion_test( set, fields ) ;
		test.reset( new SNPIdentifyingDataTestNegation( test )) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::include_snps_in_set( std::set< std::string > const& set, int fields ) {
		SNPIdentifyingDataTest::UniquePtr test = construct_snp_exclusion_test( set, fields ) ;
		add_inclusion_filter_if_necessary( "id" ) ;
		m_inclusion_filters[ "id" ]->add_subtest( test ) ;
		return *this ;
	}

	void CommonSNPFilter::add_inclusion_filter_if_necessary( std::string const& name ) {
		if( m_inclusion_filters.find( "name" ) == m_inclusion_filters.end() ) {
			SNPIdentifyingDataTestDisjunction::UniquePtr test( new SNPIdentifyingDataTestDisjunction() ) ;
			m_inclusion_filters[ name ] = test.get() ;
			m_filter.add_subtest( SNPIdentifyingDataTest::UniquePtr( test.release() )) ;
		}
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps_not_in_set( std::set< std::string > const& set, int fields ) {
		SNPIdentifyingDataTest::UniquePtr test = construct_snp_exclusion_test( set, fields ) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

	std::set< std::string > CommonSNPFilter::read_strings_from_file( std::string const& filename ) {
		std::set< std::string > result ;
		std::auto_ptr< std::istream > file = open_text_file_for_input( filename, "no_compression" ) ;
		result.insert( std::istream_iterator< std::string >( *file ), std::istream_iterator< std::string >() ) ;
		return result ;
	}
	
	SNPIdentifyingDataTest::UniquePtr CommonSNPFilter::construct_snp_exclusion_test( std::set< std::string > const& set, int fields ) {
		SNPIdentifyingDataTest::UniquePtr test ;
		switch( fields ) {
			case RSIDs:
				test.reset( new RSIDInListTest( set ) ) ;
				break ;
			case SNPIDs:
				test.reset( new SNPIDInListTest( set ) ) ;
				break ;
			case RSIDs | SNPIDs:
				test.reset( new SNPIDFieldsInListTest( set ) ) ;
				break ;
			default:
				assert(0) ;
		}
		return test ;
	}

	CommonSNPFilter& CommonSNPFilter::exclude_chromosomes_in_set( std::set< genfile::Chromosome > const& set ) {
		SNPIdentifyingDataTest::UniquePtr test( new ChromosomeInSetTest( set ) ) ;
		test.reset( new SNPIdentifyingDataTestNegation( test )) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::include_chromosomes_in_set( std::set< genfile::Chromosome > const& set ) {
		SNPIdentifyingDataTest::UniquePtr test( new ChromosomeInSetTest( set ) ) ;
		add_inclusion_filter_if_necessary( "chromosome" ) ;
		m_inclusion_filters[ "chromosome" ]->add_subtest( test ) ;
		return *this ;
	}
	
	CommonSNPFilter& CommonSNPFilter::exclude_chromosomes_not_in_set( std::set< genfile::Chromosome > const& set ) {
		SNPIdentifyingDataTest::UniquePtr test( new ChromosomeInSetTest( set ) ) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps_matching( std::string const& expression ) {
		SNPIdentifyingDataTest::UniquePtr test( new SNPIDMatchesTest( expression ) ) ;
		test.reset( new SNPIdentifyingDataTestNegation( test )) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::include_snps_matching( std::string const& expression ) {
		SNPIdentifyingDataTest::UniquePtr test( new SNPIDMatchesTest( expression ) ) ;
		add_inclusion_filter_if_necessary( "match" ) ;
		m_inclusion_filters[ "match" ]->add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps_not_matching( std::string const& expression ) {
		SNPIdentifyingDataTest::UniquePtr test( new SNPIDMatchesTest( expression ) ) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps_in_range( genfile::GenomePositionRange const& range ) {
		SNPIdentifyingDataTest::UniquePtr test( new SNPPositionInRangeTest( range ) ) ;
		test.reset( new SNPIdentifyingDataTestNegation( test )) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::include_snps_in_range( genfile::GenomePositionRange const& range ) {
		SNPIdentifyingDataTest::UniquePtr test( new SNPPositionInRangeTest( range ) ) ;
		add_inclusion_filter_if_necessary( "range" ) ;
		m_inclusion_filters[ "range" ]->add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps_not_in_range( genfile::GenomePositionRange const& range ) {
		SNPIdentifyingDataTest::UniquePtr test( new SNPPositionInRangeTest( range ) ) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

}
