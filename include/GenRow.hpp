#ifndef __GTOOL_GENROW_HPP__
#define __GTOOL_GENROW_HPP__

#include <vector>
#include <string>
#include <iostream>
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

		std::ostream& write_to_text_stream( std::ostream& aStream ) const ;
		std::istream& read_from_text_stream( std::istream& aStream ) ;
		std::ostream& write_to_binary_stream( std::ostream& aStream ) const ;
		std::istream& read_from_binary_stream( std::istream& aStream ) ;

	private:

		// Data accessors, needed for genbin support.
		void set_number_of_samples( std::size_t n ) { m_genotype_proportions.resize( n ) ; }
		void set_SNPID( std::string const& str ) { m_SNPID = str ; }
		void set_RSID( std::string const& str ) { m_RSID = str ; }
		void set_SNP_position( int pos ) { m_SNP_position = pos ; }
		void set_alleles( char a, char b ) { m_1st_allele = a; m_2nd_allele = b ; }
		void set_genotype_probabilities( std::size_t i, double aa, double ab, double bb ) { m_genotype_proportions[i] = GenotypeProportions( aa, ab, bb ) ; }
		double get_AA_probability( std::size_t i ) const { return m_genotype_proportions[i].AA() ; }
		double get_AB_probability( std::size_t i ) const { return m_genotype_proportions[i].AB() ; }
		double get_BB_probability( std::size_t i ) const { return m_genotype_proportions[i].BB() ; }

		// Data fields
		std::string m_SNPID ;
		std::string m_RSID ;
		int m_SNP_position ;
		char m_1st_allele, m_2nd_allele ;
		std::vector< GenotypeProportions > m_genotype_proportions ;

	    friend std::istream& operator>>( std::istream&, GenRow& ) ;
} ;

std::ostream& operator<<( std::ostream&, GenRow const& ) ;


#endif

