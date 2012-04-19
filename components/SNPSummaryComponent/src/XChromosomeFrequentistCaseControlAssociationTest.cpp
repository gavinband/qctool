
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/math/distributions/chi_squared.hpp>
#include <Eigen/Core>
#include "genfile/SNPIdentifyingData.hpp"
#include "snptest/case_control/LogLikelihood.hpp"
#include "snptest/case_control/NullModelLogLikelihood.hpp"
#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"
#include "integration/maximisation.hpp"
#include "components/SNPSummaryComponent/AssociationTest.hpp"
#include "components/SNPSummaryComponent/XChromosomeFrequentistCaseControlAssociationTest.hpp"

XChromosomeFrequentistCaseControlAssociationTest::XChromosomeFrequentistCaseControlAssociationTest(
	genfile::CohortIndividualSource const& samples,
	Vector const& phenotypes,
	Matrix const& covariates,
	bool const with_X_inactivation
):
	FrequentistCaseControlAssociationTest( phenotypes, covariates ),
	m_samples( samples ),
	m_sexes( get_sexes( m_samples ) ),
	m_with_X_inactivation( with_X_inactivation )
{
	m_samples_by_sex[ 'm' ] ;
	m_samples_by_sex[ 'f' ] ;
	m_samples_by_sex[ '.' ] ;

	for( std::size_t i = 0; i < m_sexes.size(); ++i ) {
		m_samples_by_sex[ m_sexes[i] ].push_back( i ) ;
	}
}

std::vector< char > XChromosomeFrequentistCaseControlAssociationTest::get_sexes( genfile::CohortIndividualSource const& samples ) const {
	std::vector< char > result( samples.get_number_of_individuals(), '.' ) ;
	genfile::CohortIndividualSource::ColumnSpec column_spec = samples.get_column_spec() ;
	if( column_spec.check_for_column( "sex" ) ) {
		if( column_spec[ "sex" ].is_discrete() ) {
			for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
				genfile::CohortIndividualSource::Entry const entry = samples.get_entry( i, "sex" ) ;
				if( !entry.is_missing() ) {
					std::string const sex = genfile::string_utils::to_lower( entry.as< std::string >() ) ;
					if( sex == "m" || sex == "male" ) {
						result[i] = 'm' ;
					} else if( sex == "f" || sex == "female" ) {
						result[i] = 'f' ;
					}
					else {
						throw genfile::MalformedInputError( samples.get_source_spec(), i+2, column_spec.find_column( "sex" )) ;
					}
				}
			}
		}
		else {
			std::cerr << "!! (XChromosomeFrequentistCaseControlAssociationTest::get_sexes): sex column found but it has the wrong type!\n" ;
			throw genfile::MalformedInputError( samples.get_source_spec(), 1, column_spec.find_column( "sex" )) ;
		}
	}
	return result ;
}

void XChromosomeFrequentistCaseControlAssociationTest::operator()(
	SNPIdentifyingData const& snp,
	Matrix const& genotypes,
	genfile::VariantDataReader&,
	ResultCallback callback
) {
	if( snp.get_position().chromosome() != genfile::Chromosome( "0X" ) ) {
		return ;
	}
	
	Matrix predictor_probs = genotypes ;
	
	// 1. figure out if males are coded like heterozygote or homozygote females.
	// Then make them like heterozygote females.
	{
		std::vector< int > const& males = m_samples_by_sex[ 'm' ] ;

		if( males.size() > 0 ) {
			if( determine_male_coding_column( snp, genotypes ) == 2 ) {
				for( std::size_t i = 0; i < males.size(); ++i ) {
					predictor_probs( males[i], 1 ) = predictor_probs( males[i], 2 ) ;
					predictor_probs( males[i], 2 ) = 0 ;
				}
			}
		}
	}
	
	// 2. If we are doing X chromosome inactivation, replace female probs here with
	// c_0 + 1/2 c_1 and 1/2 c_1 + c_2
	if( m_with_X_inactivation ) {
		std::vector< int > const& females = m_samples_by_sex['f'] ;
		for( std::size_t i = 0; i < females.size(); ++i ) {
			predictor_probs( females[i], 1 ) *= 0.5 ;
			predictor_probs( females[i], 0 ) += predictor_probs( females[i], 1 ) ;
			predictor_probs( females[i], 1 ) += predictor_probs( females[i], 2 ) ;
			predictor_probs( females[i], 2 ) = 0 ;
		}
	}
	
	// 3. get rid of anyone without sex 
	{
		std::vector< int > const& unknowns = m_samples_by_sex['f'] ;
		for( std::size_t i = 0; i < unknowns.size(); ++i ) {
			predictor_probs.row( unknowns[i] ).setZero() ;
		}
	}

	test( snp, predictor_probs, callback ) ;
}

int XChromosomeFrequentistCaseControlAssociationTest::determine_male_coding_column(
	SNPIdentifyingData const& snp,
	Matrix const& genotypes
) const {
	int column = -1 ;
	std::size_t column_determining_sample ;

	std::vector< int > const& males = m_samples_by_sex.find( 'm' )->second ;
	for( std::size_t i = 0; i < males.size(); ++i ) {
		if( genotypes( males[i], 1 ) != 0 ) {
			if( genotypes( males[i], 2 ) != 0 ) {
				std::cerr << "!! (XChromosomeFrequentistCaseControlAssociationTest::operator()): at X chromosome SNP "
					<< snp
					<< ", sample #"
					<< (males[i]+1)
					<< " (" << m_samples.get_entry( males[i], "ID_1" ) << ") "
					<< " has nonzero heterozygote and homozygote call probabilities!\n" ;
				throw genfile::BadArgumentError( "XChromosomeFrequentistCaseControlAssociationTest::operator()()", "genotypes" ) ;
			}
			else {
				if( column == -1 ) {
					column = 1 ;
					column_determining_sample = males[i] ;
				}
				else if( column != 1 ) {
					std::cerr << "!! (XChromosomeFrequentistCaseControlAssociationTest::operator()): at X chromosome SNP "
						<< snp
						<< ", samples "
						<< (column_determining_sample+1)
						<< " (" << m_samples.get_entry( column_determining_sample, "ID_1" ) << ") and "
						<< (males[i]+1)
						<< " (" << m_samples.get_entry( males[i], "ID_1" ) << ") "
						<< "are coded differently (one heterozygote, one homozygote.)\n" ;
					throw genfile::BadArgumentError( "XChromosomeFrequentistCaseControlAssociationTest::operator()()", "genotypes" ) ;
				}
			}
		}
	}
	return column ;
}
