#ifndef __GTOOL_GENROW_HPP__
#define __GTOOL_GENROW_HPP__

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "GenotypeProportions.hpp"

class GenRow
{
	public:

		std::size_t number_of_columns() const ;
 		std::size_t number_of_samples() const ;
		std::string SNPID() const { return m_SNPID ; } 
		std::string RSID() const { return m_RSID ; } 
		int SNP_position() const { return m_SNP_position ; }
		char first_allele() const { return m_1st_allele; }
		char second_allele() const { return m_2nd_allele; }
		GenotypeProportions const& genotype_proportions_for_sample( std::size_t ) const ;

		typedef std::vector< GenotypeProportions >::const_iterator genotype_proportion_iterator ;
		
        genotype_proportion_iterator begin_genotype_proportions() const { return m_genotype_proportions.begin() ; }
        genotype_proportion_iterator end_genotype_proportions() const { return m_genotype_proportions.end() ; }

		void reserveSpaceForNSamples( std::size_t ) ;

	public:
		bool operator==( GenRow const& right ) const ;
		
	private:

		std::string m_SNPID ;
		std::string m_RSID ;
		int m_SNP_position ;
		char m_1st_allele, m_2nd_allele ;
		std::vector< GenotypeProportions > m_genotype_proportions ;

	    friend std::istream& operator>>( std::istream&, GenRow& ) ;
} ;

std::ostream& operator<<( std::ostream&, GenRow const& ) ;


#endif

