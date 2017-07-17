
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <iomanip>
#include <boost/math/distributions/binomial.hpp>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/Error.hpp"
#include "components/HaplotypeFrequencyComponent/HaplotypeFrequencyComponent.hpp"

// #define DEBUG_HAPLOTYPE_FREQUENCY_LL 1
HaplotypeFrequencyLogLikelihood::HaplotypeFrequencyLogLikelihood(
	Matrix const& genotype_table,
	Matrix const& haplotype_table
):
	m_genotype_table( genotype_table ),
	m_haplotype_table( haplotype_table ),
	m_D_ll( 3 ),
	m_DDt_ll( 3, 3 )
{
	if( m_genotype_table.array().abs().maxCoeff() == 0 ) {
		throw genfile::BadArgumentError( "HaplotypeFrequencyLogLikelihood::HaplotypeFrequencyLogLikelihood()", "genotyped_table = 0" ) ;
	}
	// store derivatives of elements of the parameters pi with respect to pi.
	m_dpi.push_back( std::vector< RowVector >( 2, RowVector::Zero( 3 ) )) ;
	m_dpi.push_back( std::vector< RowVector >( 2, RowVector::Zero( 3 ) )) ;
	m_dpi[0][0] = RowVector::Constant( 3, -1 ) ;
	m_dpi[0][1]( 0 ) = 1 ;
	m_dpi[1][0]( 1 ) = 1 ;
	m_dpi[1][1]( 2 ) = 1 ;
}

HaplotypeFrequencyLogLikelihood::Vector HaplotypeFrequencyLogLikelihood::get_MLE_by_EM() const {
	//
	// This method treats the number of individuals X in the G(1,1) cell
	// that have AB/ab genotypes as binomially distributed with some unknown
	// parameter p that must be estimated.  We pick a value of p and sum over
	// possible values of X to get a new estimate of the four parameters
	// pi00, pi01, pi10, pi11.  Then p is re-computed from the parameter values.
	//
	// We stop when the parameter estimate does not change, up to some tolerance.
	//
	Matrix const& G = m_genotype_table ;
	Vector pi = estimate_parameters( 0.5 ) ;
#if DEBUG_HAPLOTYPE_FREQUENCY_LL 
	std::cerr << std::resetiosflags( std::ios::floatfield ) << "Maximising likelihood for genotypes:\n"
		<< m_genotype_table
		<< "\nand haplotypes:\n"
		<< m_haplotype_table
		<< "\n" ;
	std::cerr << "pi = " << pi(0) << " " << pi(1) << " " << pi(2) << " " << pi(3) << " " << ", p = 0.5.\n" ;
#endif
	if( G(1,1) != 0.0 ) {
		Vector old_pi ;
		std::size_t count = 0 ;
		std::size_t const max_count = 100000 ;
		double const tolerance = 0.0000001 ;
		do {
			old_pi = pi ;

			double p = ( pi(0) * pi( 3 ) ) ;
			if( p != 0 ) {
				p = p / ( p + ( pi( 1 ) * pi( 2 )) ) ;
			}

			pi = estimate_parameters( p ) ;
#if DEBUG_HAPLOTYPE_FREQUENCY_LL 
			std::cerr << "pi = " << pi(0) << " " << pi(1) << " " << pi(2) << " " << pi(3) << " " << ", p = " << p << ".\n" ;
#endif
		}
		while( ( pi - old_pi ).array().abs().maxCoeff() > tolerance && ++count < max_count ) ;
		if( count == max_count ) {
			throw genfile::OperationFailedError(
				"HaplotypeFrequencyLogLikelihood::maximise_by_EM()",
				"object of type HaplotypeFrequencyLogLikelihood",
				"convergence"
			) ;
		}
	}
	return pi.tail( 3 ) ;
}

HaplotypeFrequencyLogLikelihood::Vector HaplotypeFrequencyLogLikelihood::estimate_parameters(
	double const p
) const {
	Matrix const& G = m_genotype_table ;
	Matrix const& H = m_haplotype_table ;
	
	// p is proportion of het / het genotypes
	// that are due to 00+11 haplotype combinations.
	double expected_00_11 = G( 1, 1 ) * p ;
	double expected_01_10 = G( 1, 1 ) - expected_00_11;
	
	// params are pi00 pi01 pi10 pi11
	Vector result( 4 ) ;
	result <<
		H(0,0) + G( 0, 1 ) + 2 * G( 0, 0 ) + G( 1, 0 ) + expected_00_11,
		H(0,1) + G( 0, 1 ) + 2 * G( 0, 2 ) + G( 1, 2 ) + expected_01_10,
		H(1,0) + G( 2, 1 ) + 2 * G( 2, 0 ) + G( 1, 0 ) + expected_01_10,
		H(1,1) + G( 1, 2 ) + 2 * G( 2, 2 ) + G( 2, 1 ) + expected_00_11
	;
	result /= result.sum() ;
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
	Matrix const& H = m_haplotype_table ;
	Matrix VG( 3, 3 ) ;
	VG <<
		2.0 * lpi00,		lpi00 + lpi01,				2.0 * lpi01,
		lpi00 + lpi10,		log( het_probability ),		lpi01 + lpi11,
		2.0 * lpi10,		lpi10 + lpi11,				2.0 * lpi11
	;

	Matrix VH( 2, 2 ) ;
	VG <<
		lpi00,		lpi01,
		lpi10,		lpi11
	;

	// cells with count 0 contribute 0 to the loglikelihood/
	// Since log(0) = -inf we have to handle this specially.
	// We do this by replacing these parameter entries with 0.
	VG = VG.array() * ( G.array() > 0.0 ).cast< double >() ;
	VH = VH.array() * ( H.array() > 0.0 ).cast< double >() ;

	m_ll = (
		(G.array() * VG.array())
		+ (H.array() * VH.array())
	).sum() ;
}
double HaplotypeFrequencyLogLikelihood::get_value_of_function() const {
	return m_ll ;
}

HaplotypeFrequencyLogLikelihood::Vector HaplotypeFrequencyLogLikelihood::get_value_of_first_derivative() {
	assert(0) ;
	return m_D_ll ;
}

HaplotypeFrequencyLogLikelihood::Matrix HaplotypeFrequencyLogLikelihood::get_value_of_second_derivative() {
	assert(0) ;
	return m_DDt_ll ;
}

