#ifndef SNPIDENTIFYINGDATATEST_HPP
#define SNPIDENTIFYINGDATATEST_HPP

#include <string>
#include <memory>
#include <vector>
#include "unistd.h"
#include "genfile/GenomePosition.hpp"
#include "genfile/Chromosome.hpp"

namespace genfile {
	class SNPIdentifyingDataTest
	{
	public:
		virtual ~SNPIdentifyingDataTest() {} ;
		virtual bool operator()(
			std::string SNPID,
			std::string RSID,
			GenomePosition position,
			char first_allele,
			char second_allele
		) const = 0 ;
	protected:
		SNPIdentifyingDataTest() {} ;
	private:
		// forbid copying, assignment.
		SNPIdentifyingDataTest( SNPIdentifyingDataTest const& ) ;
		SNPIdentifyingDataTest& operator=( SNPIdentifyingDataTest const& ) ;
	} ;
	
	struct SNPIdentifyingDataTestNegation: public SNPIdentifyingDataTest
	{
		SNPIdentifyingDataTestNegation( std::auto_ptr< SNPIdentifyingDataTest > ) ;
		bool operator()(
			std::string SNPID,
			std::string RSID,
			GenomePosition position,
			char first_allele,
			char second_allele
		) const ;
	private:
		std::auto_ptr< SNPIdentifyingDataTest > m_subtest ;
	} ;
	
	struct CompoundSNPIdentifyingDataTest: public SNPIdentifyingDataTest
	{
		~CompoundSNPIdentifyingDataTest() ;
		void add_subtest( std::auto_ptr< SNPIdentifyingDataTest > subtest ) ;
		std::size_t get_number_of_subtests() const ;
		SNPIdentifyingDataTest const& get_subtest( std::size_t index ) const ;
	private:
		std::vector< SNPIdentifyingDataTest* > m_subtests ;
	} ;

	struct SNPIdentifyingDataTestConjunction: public CompoundSNPIdentifyingDataTest
	{
		bool operator()(
			std::string SNPID,
			std::string RSID,
			GenomePosition position,
			char first_allele,
			char second_allele
		) const ;
	} ;
}

#endif
