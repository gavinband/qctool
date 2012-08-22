
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "genfile/Error.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputationManager.hpp"
#include "components/SNPSummaryComponent/StratifyingSNPSummaryComputation.hpp"

SNPSummaryComputationManager::SNPSummaryComputationManager( genfile::CohortIndividualSource const& samples ):
	m_samples( samples ),
	m_sexes( get_sexes( samples )),
	m_samples_by_sex( get_samples_by_sex( m_sexes ) )
{}

std::vector< char > SNPSummaryComputationManager::get_sexes( genfile::CohortIndividualSource const& samples ) const {
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
			std::cerr << "!! (SNPSummaryComputationManager::get_sexes): sex column found but it has the wrong type!\n" ;
			throw genfile::MalformedInputError( samples.get_source_spec(), 1, column_spec.find_column( "sex" )) ;
		}
	}
	return result ;
}

std::map< char, std::vector< int > > SNPSummaryComputationManager::get_samples_by_sex( std::vector< char > const& sexes ) const {
	std::map< char, std::vector< int > > result ;
	result[ 'm' ] ;
	result[ 'f' ] ;
	result[ '.' ] ;

	for( std::size_t i = 0; i < sexes.size(); ++i ) {
		result[ sexes[i] ].push_back( i ) ;
	}

	return result ;
}

void SNPSummaryComputationManager::add_computation( std::string const& name, SNPSummaryComputation::UniquePtr computation ) {
	m_computations.insert( name, computation ) ;
}

void SNPSummaryComputationManager::add_result_callback( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}

void SNPSummaryComputationManager::begin_processing_snps( std::size_t number_of_samples ) {
	m_snp_index = 0 ;
	m_genotypes.resize( number_of_samples, 3 ) ;
}

void SNPSummaryComputationManager::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	{
		genfile::vcf::GenotypeSetter< Eigen::MatrixBase< SNPSummaryComputation::Genotypes > > setter( m_genotypes ) ;
		data_reader.get( "genotypes", setter ) ;
	}

	if( !snp.get_position().chromosome().is_missing() && snp.get_position().chromosome().is_sex_determining() ) {
		try {
			fix_sex_chromosome_genotypes( snp, m_genotypes ) ;
		}
		catch( genfile::BadArgumentError const& e ) {
			m_result_signal( snp, "comment", "Unable to determine genotype coding for males/females.  Calling may be wrong so I will treat the genotypes as missing." ) ;
			m_genotypes.setZero() ;
		}
	}

	Computations::iterator i = m_computations.begin(), end_i = m_computations.end() ;
	for( ; i != end_i; ++i ) {
		i->second->operator()(
			snp,
			m_genotypes,
			m_sexes,
			data_reader,
			boost::bind(
				boost::ref( m_result_signal ),
				snp,
				_1,
				_2
			)
		) ;
	}
	++m_snp_index ;
}

// On X and Y chromosome, recode males if necessary so the genotypes are 0/1.
// On Y chromosome check all female calls are 0.
void SNPSummaryComputationManager::fix_sex_chromosome_genotypes( genfile::SNPIdentifyingData const& snp, SNPSummaryComputation::Genotypes& genotypes ) const {
	genfile::Chromosome const& chromosome = snp.get_position().chromosome() ;
	if( chromosome != genfile::Chromosome( "0X" ) && chromosome != genfile::Chromosome( "0Y" ) ) {
		throw genfile::BadArgumentError( "SNPSummaryComputationManager::fix_sex_chromosome_genotypes()", "snp=\"" + genfile::string_utils::to_string( snp ) = "\"" ) ;
	}

	{
		std::vector< int > const& males = m_samples_by_sex.find( 'm' )->second ;
		std::cerr << "SNPSummaryComputationManager::fix_sex_chromosome_genotypes(): examining genotypes for " << males.size() << " males...\n" ;
		if( males.size() > 0 ) {
			std::cerr << "SNPSummaryComputationManager::fix_sex_chromosome_genotypes(): male coding column is " << determine_male_coding_column( snp, genotypes, males ) << "!\n" ;
			if( determine_male_coding_column( snp, genotypes, males ) == 2 ) {
				std::cerr << "SNPSummaryComputationManager::fix_sex_chromosome_genotypes(): male coding is 2...\n" ;
				for( std::size_t i = 0; i < males.size(); ++i ) {
					genotypes( males[i], 1 ) = genotypes( males[i], 2 ) ;
					genotypes( males[i], 2 ) = 0 ;
				}
			}
		}
	}
	
	if( chromosome == genfile::Chromosome( "0Y" )) {
		std::vector< int > const& females = m_samples_by_sex.find( 'f' )->second ;
		for( std::size_t i = 0; i < females.size(); ++i ) {
			if( genotypes.row( females[i] ).array().abs().maxCoeff() != 0 ) {
				std::cerr << "!! (SNPSummaryComputationManager::fix_sex_chromosome_genotypes()): at Y chromosome SNP "
					<< snp
					<< ", sample #"
					<< (females[i]+1)
					<< " (" << m_samples.get_entry( females[i], "ID_1" ) << ") "
					<< " has nonzero genotype call!\n" ;
				throw genfile::BadArgumentError( "SNPSummaryComputationManager::fix_sex_chromosome_genotypes()", "genotypes" ) ;
			}
		}
	}
}

// Figure out if males are coded like heterozygote or homozygote females.
int SNPSummaryComputationManager::determine_male_coding_column(
	genfile::SNPIdentifyingData const& snp,
	SNPSummaryComputation::Genotypes const& genotypes,
	std::vector< int > const& males
) const {
	int column = -1 ;
	std::size_t column_determining_sample ;
	for( std::size_t i = 0; i < males.size(); ++i ) {
		if( genotypes( males[i], 1 ) != 0 && genotypes( males[i], 2 ) != 0 ) {
				std::cerr << "!! (SNPSummaryComputationManager::determine_male_coding_column()): at X chromosome SNP "
					<< snp
					<< ", sample #"
					<< (males[i]+1)
					<< " (" << m_samples.get_entry( males[i], "ID_1" ) << ") "
					<< " has nonzero heterozygote and homozygote call probabilities!\n" ;
				throw genfile::BadArgumentError( "SNPSummaryComputationManager::determine_male_coding_column()", "genotypes" ) ;
		}
		
		for( int g = 1; g < 3; ++g ) {
			if( genotypes( males[i], g ) != 0 ) {
				if( column == -1 ) {
					column = g ;
					column_determining_sample = males[i] ;
				}
				else if( column != g ) {
					std::cerr << "!! (SNPSummaryComputationManager::determine_male_coding_column()): at X chromosome SNP "
						<< snp
						<< ":\n"
						<< "!! (SNPSummaryComputationManager::determine_male_coding_column()): male sample "
						<< m_samples.get_entry( column_determining_sample, "ID_1" )
						<< " (#" << (column_determining_sample+1) << ")"
						<< " is coded as a " << ( ( g == 1 ) ? "homozygote" : "heterozygote" ) << ",\n"
						<< "!! (SNPSummaryComputationManager::determine_male_coding_column()): but male sample "
						<< m_samples.get_entry( males[i], "ID_1" )
						<< " (#" << (males[i]+1) << ")"
						<< " is coded as a " << ( ( g == 1 ) ? "heterozygote" : "homozygote" )
						<< ".\n" ;
					throw genfile::BadArgumentError( "SNPSummaryComputationManager::determine_male_coding_column()", "genotypes" ) ;
				}
				break ; // no need to do both genotypes due to the check above.
			}
		}
	}
	return column ;
}

void SNPSummaryComputationManager::end_processing_snps() {}

void SNPSummaryComputationManager::stratify_by( StrataMembers const& strata, std::string const& stratification_name ) {
	Computations::iterator computation_i = m_computations.begin(), end_computations = m_computations.end() ;
	for( ; computation_i != end_computations; computation_i ) {
		Computations::iterator this_computation_i = computation_i++ ;
		std::string const this_name = this_computation_i->first ;
		SNPSummaryComputation::UniquePtr computation( m_computations.release( this_computation_i ).release() ) ;
		m_computations.insert(
			this_name,
			SNPSummaryComputation::UniquePtr( new StratifyingSNPSummaryComputation( computation, stratification_name, strata ))
		) ;
	}
}

std::string SNPSummaryComputationManager::get_summary( std::string const& prefix, std::size_t column_width ) {
	std::string result ;
	Computations::const_iterator i = m_computations.begin(), end_i = m_computations.end() ;
	for( ; i != end_i; ++i ) {
		result += i->second->get_summary( prefix, column_width ) + "\n";
	}
	return result ;
}
