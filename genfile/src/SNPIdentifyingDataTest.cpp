
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include <string>
#include <sstream>
#include <memory>
#include "unistd.h"
#include "genfile/GenomePosition.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {

	namespace {
		void append_with_comma( std::string* target, std::string const& value ) {
			(*target) += "," + value ;
		}
	}

	std::ostream& operator<<( std::ostream& out, SNPIdentifyingDataTest const& test ) {
		return out << test.display() ;
	}
	
	std::vector< std::size_t > SNPIdentifyingDataTest::get_indices_of_filtered_in_snps( std::vector< VariantIdentifyingData> const& snps ) const {
		std::vector< std::size_t > result ;
		for( std::size_t i = 0; i < snps.size(); ++i ) {
			if( this->operator()( snps[i] )) {
				result.push_back( i ) ;
			}
		}
		return result ;
	}

	std::vector< std::size_t > SNPIdentifyingDataTest::get_indices_of_filtered_in_snps( std::size_t number_of_snps, SNPGetter snps ) const {
		std::vector< std::size_t > result ;
		for( std::size_t i = 0; i < number_of_snps; ++i ) {
			if( this->operator()( snps(i) )) {
				result.push_back( i ) ;
			}
		}
		return result ;
	}
	
	SNPIdentifyingDataTestNegation::SNPIdentifyingDataTestNegation( SNPIdentifyingDataTest::UniquePtr subtest ):
	 	m_subtest( subtest )
	{
		assert( m_subtest.get() ) ;
	}

	bool SNPIdentifyingDataTestNegation::operator()( VariantIdentifyingData const& data ) const {
		return !m_subtest->operator()( data ) ;
	}
	
	std::string SNPIdentifyingDataTestNegation::display() const {
		return "NOT( " + m_subtest->display() + " )" ;
	}
	
	CompoundSNPIdentifyingDataTest::~CompoundSNPIdentifyingDataTest() {
		for( std::size_t i = 0; i < m_subtests.size(); ++i ) {
			delete m_subtests[i] ;
		}
	}

	void CompoundSNPIdentifyingDataTest::add_subtest( SNPIdentifyingDataTest::UniquePtr subtest ) {
		m_subtests.push_back( subtest.release() ) ;
	}

	std::size_t CompoundSNPIdentifyingDataTest::get_number_of_subtests() const {
		return m_subtests.size() ;
	}

	SNPIdentifyingDataTest const& CompoundSNPIdentifyingDataTest::get_subtest( std::size_t index ) const {
		assert( index < m_subtests.size() ) ;
		return *m_subtests[ index ] ;
	}

	bool SNPIdentifyingDataTestConjunction::operator()( VariantIdentifyingData const& data ) const {
		for( std::size_t i = 0; i < get_number_of_subtests(); ++i ) {
			if( !get_subtest(i)( data ) ) {
				return false ;
			}
		}
		return true ;
	}

	std::string SNPIdentifyingDataTestConjunction::display() const {
		if( get_number_of_subtests() == 0 ) {
			return "true" ;
		}
		else if( get_number_of_subtests() == 1 ) {
			return get_subtest( 0 ).display() ;
		}
		else {
			std::ostringstream result ;
			for( std::size_t i = 0; i < get_number_of_subtests(); ++i ) {
				if( i > 0 ) {
					result << " AND " ;
				}
				result << get_subtest( i ).display() ;
			}
			return result.str() ;
		}
	}

	bool SNPIdentifyingDataTestDisjunction::operator()( VariantIdentifyingData const& data ) const {
		for( std::size_t i = 0; i < get_number_of_subtests(); ++i ) {
			if( get_subtest( i )( data ) ) {
				return true ;
			}
		}
		return false ;
	}
	
	std::string SNPIdentifyingDataTestDisjunction::display() const {
		if( get_number_of_subtests() == 0 ) {
			return "false" ;
		}
		else if( get_number_of_subtests() == 1 ) {
			return get_subtest( 0 ).display() ;
		}
		else {
			std::ostringstream result ;
			for( std::size_t i = 0; i < get_number_of_subtests(); ++i ) {
				if( i > 0 ) {
					result << " OR " ;
				}
				result << get_subtest( i ).display() ;
			}
			return result.str() ;
		}
	}
	
}
