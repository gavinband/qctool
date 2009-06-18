#ifndef __GTOOL_GENOTYPEASSAYBASICSTATISTICS__
#define __GTOOL_GENOTYPEASSAYBASICSTATISTICS__

#include <vector>
#include <iostream>
#include <map>
#include "GenotypeProportions.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"

// This class holds data representing an assay of one or more SNPs.
// It has methods to return the amounts of AA, AB and BB genotypes in the sample,
// as well as genotype proportions and allele proportions. 
struct GenotypeAssayBasicStatistics
{
	public:
		GenotypeAssayBasicStatistics() {} ;
		virtual ~GenotypeAssayBasicStatistics() {} ;

		template< typename Iterator >
		void process( Iterator begin, Iterator const& end ) {
			m_number_of_samples = std::distance( begin, end ) ;
			GenotypeProportions zero_amounts( 0.0, 0.0, 0.0 ) ;
			m_genotype_amounts = std::accumulate( begin, end, zero_amounts ) ;
		}
		
		std::size_t number_of_samples() const { return m_number_of_samples ; }
		double amount_of_missing_data() const { return number_of_samples() - m_genotype_amounts.sum(); }

		GenotypeProportions const& get_genotype_amounts() const { return m_genotype_amounts ; }
		GenotypeProportions get_mean_genotype_proportions() const { return m_genotype_amounts / m_genotype_amounts.sum() ; }
		AlleleProportions get_allele_amounts() const {
			return AlleleProportions( 	2.0 * m_genotype_amounts.AA() + m_genotype_amounts.AB(),
										2.0 * m_genotype_amounts.BB() + m_genotype_amounts.AB() ) ;
		}
		AlleleProportions get_mean_allele_proportions() const { return get_allele_amounts() / (2.0 * m_genotype_amounts.sum()) ; }

		// Methods to manipulate stored genotype count.
		void floor_genotype_amounts() ;
		void round_genotype_amounts() ;

		// Methods for formatted output.
		virtual std::ostream& format_column_headers( std::ostream& ) const ;
		virtual std::ostream& format_statistic_values( std::ostream& aStream ) const ;

	private:
		std::size_t m_number_of_samples ;
		GenotypeProportions m_genotype_amounts ;
} ;

std::ostream& operator<<( std::ostream& aStream, GenotypeAssayBasicStatistics const& statistics ) ;

#endif

