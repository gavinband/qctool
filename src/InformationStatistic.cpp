#include "InformationStatistic.hpp"
#include "fill_genotype_probabilities.hpp"
#include "rescale_genotype_probabilities.hpp"

// Calculate Information statistic as per J.Marchini's email to me
// on 28/09/2009.
// Update 26/10/2010: rescale the genotype probabilities per person to avoid problem
// with info when there is missing data.  (To wit, flipping the alleles can give
// drastically different info values.)
double InformationStatistic::calculate_value( GenRow const& row ) const {

	double result ;

	if( row.number_of_samples() == 0 ) {
		return 0.0 ;
	}

	double T1 = 0.0 ;
	double T2 = 0.0 ;
	std::vector< double > e( row.number_of_samples() ) ;
	std::vector< double > f( row.number_of_samples() ) ;

	GenRow::genotype_proportion_const_iterator
		i = row.begin_genotype_proportions(),
		end_i = row.end_genotype_proportions() ;

	for( std::size_t index = 0; i != end_i ; ++i, ++index ) {
		e[index] = i->AB() + (2.0 * i->BB()) ;
		f[index] = i->AB() + (4.0 * i->BB()) ;
		T1 += i->sum() ;
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

double PlainInformationStatistic::calculate_value( GenRowStatistics const& stats ) const {
	return InformationStatistic::calculate_value( stats.row() ) ;
}

double FillingInformationStatistic::calculate_value( GenRowStatistics const& stats ) const {
	InternalStorageGenRow row( stats.row() ) ;
	fill_genotype_probabilities( &row, 0.25, 0.5, 0.25 ) ;
	return InformationStatistic::calculate_value( row ) ;
}

double ScalingInformationStatistic::calculate_value( GenRowStatistics const& stats ) const {
	InternalStorageGenRow row( stats.row() ) ;
	rescale_genotype_probabilities( &row, 0.1 ) ;
	return InformationStatistic::calculate_value( row ) ;
}

