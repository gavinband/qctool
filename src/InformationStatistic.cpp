#include <cmath>

#include "InformationStatistic.hpp"
#include "fill_genotype_probabilities.hpp"
#include "rescale_genotype_probabilities.hpp"

namespace impl {
	double add_genotype_sums( double value, GenotypeProportions const& genotypes ) {
		return value + genotypes.sum() ;
	}
}

// Calculate information statistic, updated to deal with missing genotype calls, as per J.Marchini's email to me
// on 26 October 2010.
double InformationStatistic::calculate_value( GenRow const& row ) const {
	if( row.number_of_samples() == 0 ) {
		return 0.0 ;
	}

	double non_missingness = 0.0 ;
	double theta_mle = 0.0 ;
	double v[3] = { 0.0, 0.0, 0.0 } ;
	double c[3] = { 0.0, 0.0, 0.0 } ;

	GenRow::genotype_proportion_const_iterator
		begin_genotypes = row.begin_genotype_proportions(),
		end_genotypes = row.end_genotype_proportions() ;
	
	for(
		GenRow::genotype_proportion_const_iterator i = begin_genotypes;
		i != end_genotypes;
		++i
	) {
		non_missingness += i->sum() ;
		theta_mle += i->AB() + 2.0 * i->BB() ;
		for( std::size_t g = 0; g < 3; ++g ) {
			v[g] += (*i)[g] * ( 1 - (*i)[g] ) ;
			c[g] -= (*i)[g] * (*i)[ (g+1) % 3 ] ;
		}
	}

	if( non_missingness == 0.0 ) {
		return 0.0 ;
	}
	else {
		theta_mle /= 2.0 * non_missingness ;
	}

	//
	//      F_1 + 2F_2     2F_0 + F_1
	// U =  ----------  -  ----------
	//        \theta       1 - \theta
	//
	// or
	//
	//     F_1 + 2 F_2 - 2 ( F_0 + F_1 + F_2 ) \theta
	// U = ------------------------------------------
	//            \theta ( 1 - \theta  )
	//
	// We calculated the variance as sum of variances of terms
	// plus twice the sum of covariances between terms.
	//
	// std::cerr << "jonathans_information: v0 = " << v[0] << ", v1 = " << v[1]<< ", v2 = " << v[2] << ", c01 = " << c[0] << ", c12 = " << c[1] << ", c02 = " << c[2] << "...\n" ;
	
	double vU
		// variance terms
		= 4.0 * std::pow( theta_mle, 2.0 ) * v[0]
		+ std::pow( 1.0 - 2.0 * theta_mle, 2 ) * v[1]
		+ 4.0 * std::pow( 1.0 - theta_mle, 2) * v[2]
		// covariance terms
		- ( 4.0 * theta_mle * ( 1 - 2.0 * theta_mle ) * c[0] )
		- ( 8.0 * theta_mle * ( 1 - theta_mle ) * c[2] )
		+ ( 4.0 * ( 1.0 - theta_mle ) * (1.0 - 2 * theta_mle ) * c[1] ) ;

	double const I = 1.0 - ( vU / (2 * non_missingness * theta_mle * ( 1- theta_mle )) ) ;
	// std::cerr << "jonathans_information: theta_mle = " << theta_mle << ", eI = " << expected_I << ", vU = " << vU << "...\n" ;
	return I ;
}

double InformationStatistic::calculate_value( GenRowStatistics const& stats ) const {
	return InformationStatistic::calculate_value( stats.row() ) ;
}

// Calculate original information statistic as per J.Marchini's email to me
// on 28/09/2009.
double OldInformationStatistic::calculate_value( GenRow const& row ) const {

	double result ;

	if( row.number_of_samples() == 0 ) {
		return 0.0 ;
	}

	double non_missingness = 0.0 ;
	double theta_mle = 0.0 ;
	std::vector< double > e( row.number_of_samples() ) ;
	std::vector< double > f( row.number_of_samples() ) ;
	GenRow::genotype_proportion_const_iterator
		i = row.begin_genotype_proportions(),
		end_i = row.end_genotype_proportions() ;

	for( std::size_t index = 0; i != end_i ; ++i, ++index ) {
		e[index] = i->AB() + (2.0 * i->BB()) ;
		f[index] = i->AB() + (4.0 * i->BB()) ;
		theta_mle += e[index] ;
		non_missingness += i->sum() ;
	}
	
	if( non_missingness == 0.0 ) {
		return (theta_mle == 0.0) ? 1.0 : 0.0 ;
	}

	theta_mle /= (2.0 * non_missingness) ;

	if( theta_mle == 0.0 || theta_mle == 1.0 ) {
		result = 1.0 ;
	}
	else {
		double V = 0.0 ;
		for( std::size_t i = 0; i < row.number_of_samples(); ++i ) {
			V += f[i] - (e[i]*e[i]) ;
		}
		// result =  1.0 - ( V / (2.0 * non_missingness ) * theta_mle * ( 1.0 - theta_mle ) ) ;
		result =  1.0 - ( V / (2.0 * row.number_of_samples() * theta_mle * ( 1.0 - theta_mle )) ) ;
	}
	
/*	if( result < 0.0 ) {
		result = 0.0 ;
	}
*/	
	return result ;
}

double OldInformationStatistic::calculate_value( GenRowStatistics const& stats ) const {
	return OldInformationStatistic::calculate_value( stats.row() ) ;
}

// Calculate gavin's information statistic, obtained as for the OldInformationStatistic
// but adjusted to take account of missingness in genotype calls.  See "info_measure.lyx".
double GavinsInformationStatistic::calculate_value( GenRow const& row ) const {
	std::vector< double > e( row.number_of_samples() ) ;
	std::vector< double > f( row.number_of_samples() ) ;
	std::vector< double > adjustment1( row.number_of_samples() ) ;
	std::vector< double > adjustment2( row.number_of_samples() ) ;
	
	double theta_mle = 0.0 ;
	double non_missingness = 0.0 ;
	
	GenRow::genotype_proportion_const_iterator
		i = row.begin_genotype_proportions(),
		end_i = row.end_genotype_proportions() ;

	for( std::size_t index = 0; i != end_i ; ++i, ++index ) {
		e[index] = i->AB() + (2.0 * i->BB()) ;
		f[index] = i->AB() + (4.0 * i->BB()) ;
		theta_mle += e[index] ;
		non_missingness += i->sum() ;
		adjustment1[index] = (1.0 - i->sum()) * e[index] ; // probability of null call times expected.
		adjustment2[index] = (1.0 - i->sum()) * ( 1.0 - i->sum() ) ; // square of null call times.
	}
	
	if( non_missingness == 0.0 ) {
		return 0.0 ;
	}
	
	theta_mle /= 2.0 * non_missingness ;

	if( theta_mle == 0.0 || theta_mle == 1.0 ) {
		return 1.0 ;
	}

	double variance = 0.0 ;
	for( std::size_t i = 0; i < e.size(); ++i ) {
		variance += f[ i ] - ( e[ i ] * e[ i ] ) ;
		// So far as with the old statistic.  Now do the adjustments
		variance -= 4.0 * theta_mle * adjustment1[ i ] ;
		variance -= 4.0 * theta_mle * theta_mle * adjustment2[ i ] ;
		// std::cerr << "i = " << i << ", adjustment1 is " << adjustment1[i] << ", adjustment2 is " << adjustment2[ i ] << ".\n" ;
	}
	double const missingness = row.number_of_samples() - non_missingness ;

	variance += missingness * 2.0 * theta_mle * ( 1.0 + theta_mle ) ;

	double denominator = 2.0 * row.number_of_samples() * theta_mle * ( 1.0 - theta_mle ) ;
	double result = 1.0 - ( variance / denominator ) ;
	
	// std::cerr << "missingness is " << missingness << ", variance is " << variance << "`...\n" ;

	return result ;
}

double GavinsInformationStatistic::calculate_value( GenRowStatistics const& stats ) const {
	return GavinsInformationStatistic::calculate_value( stats.row() ) ;
}
