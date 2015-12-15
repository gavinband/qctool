
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <iterator>

#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/CommonSNPFilter.hpp"
#include "genfile/RSIDInListTest.hpp"
#include "genfile/SNPIDInListTest.hpp"
#include "genfile/SNPIDFieldsInListTest.hpp"
#include "genfile/ChromosomeInSetTest.hpp"
#include "genfile/SNPIDMatchesTest.hpp"
#include "genfile/SNPInSetTest.hpp"
#include "genfile/SNPPositionInRangeTest.hpp"
#include "genfile/PositionInListTest.hpp"
#include "genfile/snp_data_utils.hpp"

namespace genfile {
	CommonSNPFilter::CommonSNPFilter() {
	}

	bool CommonSNPFilter::operator()( VariantIdentifyingData const& data ) const {
		return m_filter( data ) ;
	}
	
	std::string CommonSNPFilter::display() const { return m_filter.display() ; }
	
	CommonSNPFilter& CommonSNPFilter::exclude_snps_in_file( std::string const& filename, int fields ) {
		exclude_snps_in_set( read_strings_from_file( filename ), fields ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::include_snps_in_file( std::string const& filename, int fields ) {
		include_snps_in_set( read_strings_from_file( filename ), fields ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps_not_in_file( std::string const& filename, int fields ) {
		return exclude_snps_not_in_set( read_strings_from_file( filename ), fields ) ;
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps_in_set( std::set< std::string > const& set, int fields ) {
		SNPIdentifyingDataTest::UniquePtr test = construct_snp_inclusion_test( set, fields ) ;
		test.reset( new SNPIdentifyingDataTestNegation( test )) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::include_snps_in_set( std::set< std::string > const& set, int fields ) {
		SNPIdentifyingDataTest::UniquePtr test = construct_snp_inclusion_test( set, fields ) ;
		add_inclusion_filter_if_necessary( "id" ) ;
		m_inclusion_filters[ "id" ]->add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps( std::vector< VariantIdentifyingData > const& snps, VariantIdentifyingData::CompareFields const& comparer ) {
		SNPIdentifyingDataTest::UniquePtr test( new SNPInSetTest( std::set< VariantIdentifyingData >( snps.begin(), snps.end() ), comparer ) ) ;
		test.reset( new SNPIdentifyingDataTestNegation( test )) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

	CommonSNPFilter& CommonSNPFilter::include_snps( std::vector< VariantIdentifyingData > const& snps, VariantIdentifyingData::CompareFields const& comparer ) {
		SNPIdentifyingDataTest::UniquePtr test( new SNPInSetTest( std::set< VariantIdentifyingData >( snps.begin(), snps.end() ), comparer ) ) ;
		add_inclusion_filter_if_necessary( "snps" ) ;
		m_inclusion_filters[ "snps" ]->add_subtest( test ) ;
		return *this ;
	}

	void CommonSNPFilter::add_inclusion_filter_if_necessary( std::string const& name ) {
		if( m_inclusion_filters.find( name ) == m_inclusion_filters.end() ) {
			SNPIdentifyingDataTestDisjunction::UniquePtr test( new SNPIdentifyingDataTestDisjunction() ) ;
			m_inclusion_filters[ name ] = test.get() ;
			m_filter.add_subtest( SNPIdentifyingDataTest::UniquePtr( test.release() )) ;
		}
	}

	CommonSNPFilter& CommonSNPFilter::exclude_snps_not_in_set( std::set< std::string > const& set, int fields ) {
		SNPIdentifyingDataTest::UniquePtr test = construct_snp_inclusion_test( set, fields ) ;
		m_filter.add_subtest( test ) ;
		return *this ;
	}

	std::set< std::string > CommonSNPFilter::read_strings_from_file( std::string const& filename ) {
		std::set< std::string > result ;
		std::auto_ptr< std::istream > file = open_text_file_for_input( filename, "no_compression" ) ;
		result.insert( std::istream_iterator< std::string >( *file ), std::istream_iterator< std::string >() ) ;
		return result ;
	}
	
	SNPIdentifyingDataTest::UniquePtr CommonSNPFilter::construct_snp_inclusion_test( std::set< std::string > const& set, int fields ) {
		SNPIdentifyingDataTest::UniquePtr test ;
		std::set< GenomePosition > position_set ;
		std::set< std::string >::const_iterator i = set.begin() ;
		std::set< std::string >::const_iterator const end_i = set.end() ;
		
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
			case Positions:
				for( ; i != end_i; ++i ) {
					position_set.insert( GenomePosition( *i ) ) ;
				}
				test.reset( new PositionInListTest( position_set )) ;
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
