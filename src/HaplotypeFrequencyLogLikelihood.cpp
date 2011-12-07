#include <string>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/Error.hpp"
#include "HaplotypeFrequencyComponent.hpp"

HaplotypeFrequencyLogLikelihood::HaplotypeFrequencyLogLikelihood( Matrix const& genotype_table ):
	m_genotype_table( genotype_table ),
	m_D_ll( 3 ),
	m_DDt_ll( 3, 3 )
{
	// store derivatives of elements of the parameters pi with respect to pi.
	m_dpi.push_back( std::vector< RowVector >( 2, RowVector::Zero( 3 ) )) ;
	m_dpi.push_back( std::vector< RowVector >( 2, RowVector::Zero( 3 ) )) ;
	m_dpi[0][0] = RowVector::Constant( 3, -1 ) ;
	m_dpi[0][1]( 0 ) = 1 ;
	m_dpi[1][0]( 1 ) = 1 ;
	m_dpi[1][1]( 2 ) = 1 ;
}

HaplotypeFrequencyLogLikelihood::Vector HaplotypeFrequencyLogLikelihood::get_MLE_by_EM() const {
	// Only missing data is resolution of middle cell
	// Call the alleles at the first SNP A and a and at the second SNP B and b.
	Matrix const& G = m_genotype_table ;
	double AB_ab = G(1,1) / 2 ;
	Vector pi = estimate_parameters( AB_ab ) ;
	Vector old_pi ;
	std::size_t count = 0 ;
	std::size_t const max_count = 100 ;
	double const tolerance = 0.0000001 ;
	do {
		old_pi = pi ;
		double pi00 = 1.0 - pi(0) - pi(1) - pi(2) ;
		AB_ab = G(1,1) * ( pi00 * pi( 2 ) ) / ( pi00 * pi(2) + pi( 0 ) * pi( 1 )) ;
		pi = estimate_parameters( AB_ab ) ;
	}
	while( ( pi - old_pi ).squaredNorm() > tolerance && ++count < max_count ) ;
	if( count == max_count ) {
		throw genfile::OperationFailedError(
			"HaplotypeFrequencyLogLikelihood::maximise_by_EM()",
			"object of type HaplotypeFrequencyLogLikelihood",
			"convergence"
		) ;
	}
	return pi ;
}

HaplotypeFrequencyLogLikelihood::Vector HaplotypeFrequencyLogLikelihood::estimate_parameters(
	double const AB_ab
) const {
	Matrix const& G = m_genotype_table ;
	double const Ab_aB = G( 1, 1 ) - AB_ab ;
	Vector result( 3 ) ;
	result <<
		G( 0, 1 ) + 2 * G( 0, 2 ) + G( 1, 2 ) + Ab_aB,
		G( 2, 1 ) + 2 * G( 2, 0 ) + G( 1, 0 ) + Ab_aB,
		G( 0, 1 ) + 2 * G( 2, 2 ) + G( 1, 2 ) + AB_ab
	;
	result /= 2.0 * G.sum() ;
	return result ;
}

void HaplotypeFrequencyLogLikelihood::evaluate_at( Vector const& pi ) {
	double const pi00 = 1.0 - pi.sum() ;
	double const& pi01 = pi(0) ;
	double const& pi10 = pi(1) ;
	double const& pi11 = pi(2) ;
	
	using std::log ;
	double const lpi00 = log( pi00 ) ;
	double const lpi01 = log( pi01 ) ;
	double const lpi10 = log( pi10 ) ;
	double const lpi11 = log( pi11 ) ;
	
	double het_probability = pi00 * pi11 + pi01 * pi10 ;

	Matrix const& G = m_genotype_table ;
	Matrix V( 3, 3 ) ;
	V <<
		2.0 * lpi00,		lpi00 + lpi01,				2.0 * lpi01,
		lpi00 + lpi10,		log( het_probability ),		lpi01 + lpi11,
		2.0 * lpi10,		lpi10 + lpi11,				2.0 * lpi11 ;

	m_ll = 0.0 ;
	for( int i = 0; i < 3; ++i ) {
		for( int j = 0; j < 3; ++j ) {
			if( G( i, j ) != 0.0 ) {
				m_ll += G( i, j ) * V( i, j ) ;
			}
		}
	}

/*	m_D_ll
		+ G( 1, 0 ) * ( ( m_dpi[1][0] / pi10 ) )
		+ G( 1, 1 ) * ( ( m_dpi[0][0] * pi11 + m_dpi[1][1] * pi00 + m_dpi[0][1] * pi10 + m_dpi[1][0] * pi01 ) / ( pi00 * pi11 + pi01 * pi10 ) )
		+ G( 1, 2 ) * (  ( m_dpi[1][1] / pi11 ) )
		+ G( 2, 0 ) * 2.0 * ( m_dpi[1][0] / pi10 )
		+ G( 2, 1 ) * ( ( m_dpi[1][0] / pi10 ) + ( m_dpi[1][1] / pi11 ) )
		+ G( 2, 2 ) * 2.0 * ( m_dpi[1][1] / pi11 ) ;
*/
	Matrix c( 2, 2 ) ;
	c <<
		( 2.0 * G( 0, 0 ) + G( 0, 1 ) + G( 1, 0 ) ), 		( 2.0 * G( 0, 2 ) + G( 0, 1 ) + G( 1, 2 ) ),
		( 2.0 * G( 2, 0 ) + G( 1, 0 ) + G( 2, 1 ) ),		( 2.0 * G( 2, 2 ) + G( 1, 2 ) + G( 2, 1 ) ) ;
	RowVector Dhet_probability = ( m_dpi[0][0] * pi11 + m_dpi[1][1] * pi00 + m_dpi[0][1] * pi10 + m_dpi[1][0] * pi01 ) ;
	m_D_ll.setZero() ; // derivative of loglikelihood.
	m_DDt_ll.setZero() ; // derivative of (transpose of) derivative
	if( c( 0, 0 ) != 0.0 ) {
		m_D_ll += c( 0, 0 ) * ( m_dpi[0][0] / pi00 ) ;
		m_DDt_ll += c( 0, 0 ) * -m_dpi[0][0].transpose() * m_dpi[0][0] / ( pi00 * pi00 ) ;
	}
	if( c( 0, 1 ) != 0.0 ) {
		m_D_ll += c( 0, 1 ) * ( m_dpi[0][1] / pi01 ) ;
		m_DDt_ll += c( 0, 1 ) * -m_dpi[0][1].transpose() * m_dpi[0][1] / ( pi01 * pi01 ) ;
	}
	if( c( 1, 0 ) != 0.0 ) {
		m_D_ll += c( 1, 0 ) * ( m_dpi[1][0] / pi10 ) ;
		m_DDt_ll += c( 1, 0 ) * -m_dpi[1][0].transpose() * m_dpi[1][0] / ( pi10 * pi10 ) ;
	}
	if( c( 1, 1 ) != 0.0 ) {
		m_D_ll += c( 1, 1 ) * ( m_dpi[1][1] / pi11 ) ;
		m_DDt_ll += c( 1, 1 ) * -m_dpi[1][1].transpose() * m_dpi[1][1] / ( pi11 * pi11 ) ;
	}
	if( G( 1, 1 ) != 0.0 ) {
		m_D_ll += G( 1, 1 ) * ( Dhet_probability / het_probability ) ;
		m_DDt_ll += G( 1, 1 ) * (
			(
				- m_dpi[1][1].transpose() * m_dpi[0][0] / ( pi11 * pi11 )
				- m_dpi[0][0].transpose() * m_dpi[1][1] / ( pi00 * pi00 )
				- m_dpi[1][0].transpose() * m_dpi[0][1] / ( pi10 * pi10 )
				- m_dpi[0][1].transpose() * m_dpi[1][0] / ( pi01 * pi01 )
			) / het_probability
			- (
				Dhet_probability.transpose() * Dhet_probability ) / ( het_probability * het_probability )
		) ;
	}
}
double HaplotypeFrequencyLogLikelihood::get_value_of_function() const {
	return m_ll ;
}

HaplotypeFrequencyLogLikelihood::Vector HaplotypeFrequencyLogLikelihood::get_value_of_first_derivative() {
	return m_D_ll.transpose() ;
}

HaplotypeFrequencyLogLikelihood::Matrix HaplotypeFrequencyLogLikelihood::get_value_of_second_derivative() {
	return m_DDt_ll ;
}

