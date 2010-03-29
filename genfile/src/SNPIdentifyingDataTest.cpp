#include <cassert>
#include <string>
#include <memory>
#include "unistd.h"
#include "genfile/GenomePosition.hpp"
#include "genfile/Chromosome.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {
	SNPIdentifyingDataTestNegation::SNPIdentifyingDataTestNegation( std::auto_ptr< SNPIdentifyingDataTest > subtest ):
	 	m_subtest( subtest )
	{
		assert( m_subtest.get() ) ;
	}

	bool SNPIdentifyingDataTestNegation::operator()(
		std::string SNPID,
		std::string RSID,
		GenomePosition position,
		char first_allele,
		char second_allele
	) const {
		return !m_subtest->operator()( SNPID, RSID, position, first_allele, second_allele ) ;
	}
	
	CompoundSNPIdentifyingDataTest::~CompoundSNPIdentifyingDataTest() {
		for( std::size_t i = 0; i < m_subtests.size(); ++i ) {
			delete m_subtests[i] ;
		}
	}

	void CompoundSNPIdentifyingDataTest::add_subtest( std::auto_ptr< SNPIdentifyingDataTest > subtest ) {
		m_subtests.push_back( subtest.release() ) ;
	}

	std::size_t CompoundSNPIdentifyingDataTest::get_number_of_subtests() const {
		return m_subtests.size() ;
	}

	SNPIdentifyingDataTest const& CompoundSNPIdentifyingDataTest::get_subtest( std::size_t index ) const {
		assert( index < m_subtests.size() ) ;
		return *m_subtests[ index ] ;
	}

	bool SNPIdentifyingDataTestConjunction::operator()(
		std::string SNPID,
		std::string RSID,
		GenomePosition position,
		char first_allele,
		char second_allele
	) const {
		for( std::size_t i = 0; i < get_number_of_subtests(); ++i ) {
			if( !get_subtest(i)(
				SNPID,
				RSID,
				position,
				first_allele,
				second_allele
				)
			) {
				return false ;
			}
		}
		return true ;
	}
}
