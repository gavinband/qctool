
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_COMMON_SNP_FILTER_HPP
#define GENFILE_COMMON_SNP_FILTER_HPP

#include <string>
#include <memory>
#include <set>
#include <string>
#include "genfile/GenomePosition.hpp"
#include "genfile/GenomePositionRange.hpp"
#include "genfile/SNPDataSource.hpp"
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
			std::string first_allele,
			std::string second_allele
		) const ;

		std::string display() const ;

	public:
		CommonSNPFilter() ;
		
		enum Field { RSIDs = 1, SNPIDs = 2, Positions = 4 } ;
		
		CommonSNPFilter& exclude_snps_in_file( std::string const& filename, int fields ) ;
		CommonSNPFilter& exclude_snps_not_in_file( std::string const& filename, int fields ) ;
		CommonSNPFilter& include_snps_in_file( std::string const& filename, int fields ) ;
		
		CommonSNPFilter& exclude_snps_in_set( std::set< std::string > const& set, int fields ) ;
		CommonSNPFilter& exclude_snps_not_in_set( std::set< std::string > const& set, int fields ) ;
		CommonSNPFilter& include_snps_in_set( std::set< std::string > const& set, int fields ) ;

		CommonSNPFilter& exclude_snps( std::vector< SNPIdentifyingData > const& snps ) ;
		CommonSNPFilter& include_snps( std::vector< SNPIdentifyingData > const& snps ) ;
		
		CommonSNPFilter& exclude_chromosomes_in_set( std::set< genfile::Chromosome > const& set ) ;
		CommonSNPFilter& exclude_chromosomes_not_in_set( std::set< genfile::Chromosome > const& set ) ;
		CommonSNPFilter& include_chromosomes_in_set( std::set< genfile::Chromosome > const& set ) ;
		
		CommonSNPFilter& exclude_snps_matching( std::string const& expression ) ;
		CommonSNPFilter& exclude_snps_not_matching( std::string const& expression ) ;
		CommonSNPFilter& include_snps_matching( std::string const& expression ) ;

		CommonSNPFilter& exclude_snps_in_range( genfile::GenomePositionRange const& range ) ;
		CommonSNPFilter& exclude_snps_not_in_range( genfile::GenomePositionRange const& range ) ;
		CommonSNPFilter& include_snps_in_range( genfile::GenomePositionRange const& range ) ;

	private:
		SNPIdentifyingDataTestConjunction m_filter ;
		std::map< std::string, CompoundSNPIdentifyingDataTest* > m_inclusion_filters ;
	private:
		SNPIdentifyingDataTest::UniquePtr construct_snp_inclusion_test( std::set< std::string > const& set, int fields ) ;
		void add_inclusion_filter_if_necessary( std::string const& name ) ;
		std::set< std::string > read_strings_from_file( std::string const& filename ) ;
	} ;
}

#endif
