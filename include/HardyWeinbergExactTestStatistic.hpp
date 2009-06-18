#ifndef __HARDY_WEINBERG_EXACT_TEST__
#define __HARDY_WEINBERG_EXACT_TEST__

#include "GenotypeAssayStatistics.hpp"

// Implementation of the "exact" test for Hardy-Weinberg equilibrium
// as described in 
// (*) Wigginton et al "A Note on Exact Tests of Hardy-Weinberg Equilibrium", Am. J. Hum. Genet. (2005).
//
// This class implements the test in a direct way, i.e. directly using formula (1) of the paper, rather than
// via the recurrence relations used in their implementation.
//
// Note: notation in this class comes from that in the paper (*) above.  Thus "a" always denotes the minor allele.
class HardyWeinbergExactTestStatistic: public GenotypeAssayStatistic
{
	public:
	
		double calculate_value( GenotypeAssayStatistics const& statistics ) const ;
	
	protected:
	
		void calculate_allele_amounts() const ;
		
		// The next two functions implement formula (1) of the Wigginton paper.
		// To allow efficient calculation, we split up the formula into the nab-invariant
		// and nab-dependent parts
		double calculate_log_of_invariant_part_of_formula1( double const observed_na, double const observed_nb ) const ;
		double calculate_probability_of_nab_heterozygotes( double const observed_na, double const observed_nb, double const nab, double const log_of_nab_invariant_part ) const ;

		// Implementation of P_{HWE}
		double calculate_probability_of_HWE( double const observed_na, double const observed_nb, double const observed_nab ) const ;
	
	private:
	
		double m_observed_na ;
		double m_observed_nb ;
		double m_observed_N ;
		double m_observed_nab ;
} ;

#endif
