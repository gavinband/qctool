#ifndef GENFILE_COMMON_SNP_FILTER_HPP
#define GENFILE_COMMON_SNP_FILTER_HPP

#include <string>
#include <memory>
#include <set>
#include <string>
#include "genfile/GenomePosition.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"

namespace genfile {
	class CommonSNPFilter: public SNPIdentifyingDataTest
		/*
		* This class acts as an easier-to-use frontend for a conjunction
		* of common SNPIdentifyingDataTest types.
		*/
	{
	public:
		typedef std::auto_ptr< CommonSNPFilter > UniquePtr ;
		
		using SNPIdentifyingDataTest::operator() ;
		bool operator()(
			std::string SNPID,
			std::string RSID,
			GenomePosition position,
			char first_allele,
			char second_allele
		) const ;

		std::string display() const ;

	public:
		CommonSNPFilter() ;
		
		enum Field { RSIDs = 1, SNPIDs = 2 } ;
		
		CommonSNPFilter& exclude_snps_in_file( std::string const& filename, int fields ) ;
		CommonSNPFilter& exclude_snps_not_in_file( std::string const& filename, int fields ) ;
		
		CommonSNPFilter& exclude_snps_in_set( std::set< std::string > const& set, int fields ) ;
		CommonSNPFilter& exclude_snps_not_in_set( std::set< std::string > const& set, int fields ) ;
		
		CommonSNPFilter& exclude_chromosomes_in_set( std::set< genfile::Chromosome > const& set ) ;
		CommonSNPFilter& exclude_chromosomes_not_in_set( std::set< genfile::Chromosome > const& set ) ;
		
		CommonSNPFilter& exclude_snps_matching( std::string const& expression ) ;
		CommonSNPFilter& exclude_snps_not_matching( std::string const& expression ) ;

	private:
		SNPIdentifyingDataTestConjunction m_filter ;

	private:
		SNPIdentifyingDataTest::UniquePtr construct_test( std::set< std::string > const& set, int fields ) ;
		std::set< std::string > read_strings_from_file( std::string const& filename ) ;
	} ;
}

#endif
