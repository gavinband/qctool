#include "InformationStatistic.hpp"

// Calculate Information statistic as per J.Marchini's email to me
// on 28/09/2009.
double InformationStatistic::calculate_value( GenRowStatistics const& stats ) const {

	double result ;

	GenRow const& row = stats.row() ;

	if( row.number_of_samples() == 0 ) {
		return 0.0 ;
	}

	double T1 = stats.get_genotype_amounts().sum() ;
	double T2 = 0.0 ;
	std::vector< double > e( row.number_of_samples() ) ;
	std::vector< double > f( row.number_of_samples() ) ;

	GenRow::genotype_proportion_const_iterator
		i = row.begin_genotype_proportions(),
		end_i = row.end_genotype_proportions() ;

	for( std::size_t index = 0; i != end_i ; ++i, ++index ) {
		e[index] = i->AB() + (2.0 * i->BB()) ;
		f[index] = i->AB() + (4.0 * i->BB()) ;
		T2 += e[index] ;
	}

	if( T1 == 0.0 ) {
		return (T2 == 0.0) ? 1.0 : 0.0 ;
	}

	double E = T2 / (2.0 * T1) ;

	if( E == 0 || E == 1 ) {
		result = 1.0 ;
	}
	else {
		double E_times_1_minus_E = E * ( 1.0 - E ) ;
		double I1 = (2.0 * T1) / E_times_1_minus_E ;
		double V = 0.0 ;
		for( std::size_t i = 0; i < row.number_of_samples(); ++i ) {
			V += f[i] - (e[i]*e[i]) ;
		}
		V /= E_times_1_minus_E * E_times_1_minus_E ;
		result =  (I1 - V) / I1 ;
	}
	
	if( result < 0.0 ) {
		result = 0.0 ;
	}
	
	return result ;
}

