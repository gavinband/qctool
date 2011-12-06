#include <string>
#include <boost/function.hpp>
#include <Eigen/Core>
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
	m_ll
		= G( 0, 0 ) * 2.0 * lpi00 
		+ G( 0, 1 ) * ( lpi00 + lpi01 )
		+ G( 0, 2 ) * 2.0 * ( lpi01 )
		+ G( 1, 0 ) * ( lpi00 + lpi10 )
		+ G( 1, 1 ) * log( het_probability )
		+ G( 1, 2 ) * ( lpi01 + lpi11 )
		+ G( 2, 0 ) * 2.0 * lpi10
		+ G( 2, 1 ) * ( lpi10 + lpi11 )
		+ G( 2, 2 ) * 2.0 * lpi11 ;

/*	m_D_ll
		+ G( 1, 0 ) * ( ( m_dpi[1][0] / pi10 ) )
		+ G( 1, 1 ) * ( ( m_dpi[0][0] * pi11 + m_dpi[1][1] * pi00 + m_dpi[0][1] * pi10 + m_dpi[1][0] * pi01 ) / ( pi00 * pi11 + pi01 * pi10 ) )
		+ G( 1, 2 ) * (  ( m_dpi[1][1] / pi11 ) )
		+ G( 2, 0 ) * 2.0 * ( m_dpi[1][0] / pi10 )
		+ G( 2, 1 ) * ( ( m_dpi[1][0] / pi10 ) + ( m_dpi[1][1] / pi11 ) )
		+ G( 2, 2 ) * 2.0 * ( m_dpi[1][1] / pi11 ) ;
*/
	double c00 = ( 2.0 * G( 0, 0 ) + G( 0, 1 ) + G( 1, 0 ) ) ;
	double c01 = ( 2.0 * G( 0, 2 ) + G( 0, 1 ) + G( 1, 2 ) ) ;
	double c10 = ( 2.0 * G( 2, 0 ) + G( 1, 0 ) + G( 2, 1 ) ) ;
	double c11 = ( 2.0 * G( 2, 2 ) + G( 1, 2 ) + G( 2, 1 ) ) ;
	RowVector Dhet_probability = ( m_dpi[0][0] * pi11 + m_dpi[1][1] * pi00 + m_dpi[0][1] * pi10 + m_dpi[1][0] * pi01 ) ;
	m_D_ll // derivative of loglikelihood.
		= c00 * ( m_dpi[0][0] / pi00 )
		+ c01 * ( m_dpi[0][1] / pi01 )
		+ c10 * ( m_dpi[1][0] / pi10 )
		+ c11 * ( m_dpi[1][1] / pi11 )
		+ G( 1, 1 ) * ( Dhet_probability / het_probability )
	;

	m_DDt_ll
		= c00 * -m_dpi[0][0].transpose() * m_dpi[0][0] / ( pi00 * pi00 )
		+ c01 * -m_dpi[0][1].transpose() * m_dpi[0][1] / ( pi01 * pi01 )
		+ c10 * -m_dpi[1][0].transpose() * m_dpi[1][0] / ( pi10 * pi10 )
		+ c11 * -m_dpi[1][1].transpose() * m_dpi[1][1] / ( pi11 * pi11 )
		+ G( 1, 1 ) * (
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
double HaplotypeFrequencyLogLikelihood::get_value_of_function() const {
	return m_ll ;
}

HaplotypeFrequencyLogLikelihood::Vector HaplotypeFrequencyLogLikelihood::get_value_of_first_derivative() {
	return m_D_ll.transpose() ;
}

HaplotypeFrequencyLogLikelihood::Matrix HaplotypeFrequencyLogLikelihood::get_value_of_second_derivative() {
	return m_DDt_ll ;
}

