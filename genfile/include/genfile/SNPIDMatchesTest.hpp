#ifndef GENFILE_SNPIDMATCHES_TEST_HPP
#define GENFILE_SNPIDMATCHES_TEST_HPP

#include <string>
#include "SNPIdentifyingDataTest.hpp"

namespace genfile {
	struct SNPIDMatchesTest: public SNPIdentifyingDataTest
	{
	public:
		SNPIDMatchesTest( std::string const& expression ) ;

		bool operator()(
			std::string SNPID,
			std::string RSID,
			GenomePosition position,
			std::string first_allele,
			std::string second_allele
		) const ;
		
		std::string display() const ;

	private:
		void setup( std::string const& ) ;
		bool match( std::string const& ) const ;
	private:
		enum Type { eSNPID = 1, eRSID = 2, eEITHER = 3 } ;	
		Type m_type ;
		std::string m_expression, m_prefix, m_suffix ;
		char m_wildcard_char ;
		bool m_have_wildcard ;
	} ;
}

#endif

