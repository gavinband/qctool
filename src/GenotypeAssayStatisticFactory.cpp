#include "GenotypeAssayStatisticFactory.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "GenotypeAssayStatisticArithmetic.hpp"
#include "SimpleGenotypeAssayStatistics.hpp"
#include "GenRowStatistics.hpp"
#include "SampleRowStatistics.hpp"
#include "GenotypeAssayStatisticArithmetic.hpp"
#include "SNPHWE.hpp"
#include "HardyWeinbergExactTestStatistic.hpp"
#include "LikelihoodRatioTestStatistic.hpp"	
#include "InformationStatistic.hpp"
#include "string_utils/string_utils.hpp"


void GenotypeAssayStatisticFactory::add_statistics( std::vector< std::string > statistic_specs, GenotypeAssayStatistics& statistics ) {
	for( std::size_t i = 0; i < statistic_specs.size(); ++i ) {
		statistics.add_statistic( statistic_specs[i], GenotypeAssayStatisticFactory::create_statistic( statistic_specs[i] )) ;
	}
}

std::auto_ptr< GenotypeAssayStatistic > GenotypeAssayStatisticFactory::create_statistic( std::string statistic_spec ) {
	statistic_spec = string_utils::strip( statistic_spec ) ;

	// Handle statistic arithmetic
	std::vector< std::string > bits = string_utils::split( statistic_spec, "/" ) ;
	if (bits.size() == 2) {
		try {
			std::auto_ptr< GenotypeAssayStatistic > first_stat_ptr, second_stat_ptr ;
			first_stat_ptr = GenotypeAssayStatisticFactory::create_statistic( bits[0] ) ;
			second_stat_ptr = GenotypeAssayStatisticFactory::create_statistic( bits[1] ) ;
			return std::auto_ptr< GenotypeAssayStatistic >( new StatisticRatio( first_stat_ptr, second_stat_ptr )) ;
		}
		catch( StatisticNotFoundException const& e ) {
		}
	}	

	bits = string_utils::split( statistic_spec, "*" ) ;
	if (bits.size() == 2) {
		try {
			std::auto_ptr< GenotypeAssayStatistic > first_stat_ptr, second_stat_ptr ;
			first_stat_ptr = GenotypeAssayStatisticFactory::create_statistic( bits[0] ) ;
			second_stat_ptr = GenotypeAssayStatisticFactory::create_statistic( bits[1] ) ;
			return std::auto_ptr< GenotypeAssayStatistic >( new StatisticProduct( first_stat_ptr, second_stat_ptr )) ;
		}
		catch( StatisticNotFoundException const& e ) {
		}
	}	

	bits = string_utils::split( statistic_spec, "+" ) ;
	if (bits.size() == 2) {
		std::auto_ptr< GenotypeAssayStatistic > first_stat_ptr, second_stat_ptr ;
		try {
			first_stat_ptr = GenotypeAssayStatisticFactory::create_statistic( bits[0] ) ;
			second_stat_ptr = GenotypeAssayStatisticFactory::create_statistic( bits[1] ) ;
			return std::auto_ptr< GenotypeAssayStatistic >( new StatisticSum( first_stat_ptr, second_stat_ptr )) ;
		}
		catch( StatisticNotFoundException const& e ) {
		}
	}	

	bits = string_utils::split( statistic_spec, "-" ) ;
	if (bits.size() == 2) {
		try {
			std::auto_ptr< GenotypeAssayStatistic > first_stat_ptr, second_stat_ptr ;
			first_stat_ptr = GenotypeAssayStatisticFactory::create_statistic( bits[0] ) ;
			second_stat_ptr = GenotypeAssayStatisticFactory::create_statistic( bits[1] ) ;
			return std::auto_ptr< GenotypeAssayStatistic >( new StatisticDifference( first_stat_ptr, second_stat_ptr )) ;
		}
		catch( StatisticNotFoundException const& e ) {
		}
	}	
	
	// See if it's a special type.
	// TODO: rework this dispatching mechanism.
	std::auto_ptr< GenotypeAssayStatistic > result = GenRowStatisticFactory::create_statistic( statistic_spec ) ;
	if( !result.get() ) {
		result = SampleRowStatisticFactory::create_statistic( statistic_spec ) ;
	}
	if( result.get() ) {
		return result ;
	}

	// if none of the above, just create a single stat.
	
	// Stats from whiteboard 26/06/09:
	if( statistic_spec == "MAF" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new MinorAlleleProportionStatistic ) ;
	}
	else if( statistic_spec == "HWE" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new MinusLog10SNPHWEStatistic ) ;
	}
	else if( statistic_spec == "information" ) {
		return std::auto_ptr< GenotypeAssayStatistic > ( new PlainInformationStatistic() ) ;
	}
	else if( statistic_spec == "filled_information" ) {
		return std::auto_ptr< GenotypeAssayStatistic > ( new FillingInformationStatistic()) ;
	}
	else if( statistic_spec == "scaled_information" ) {
		return std::auto_ptr< GenotypeAssayStatistic > ( new ScalingInformationStatistic()) ;
	}
	else if( statistic_spec == "missing" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new MissingDataProportionStatistic ) ; 
	}
	else if( statistic_spec == "heterozygosity" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new HeterozygosityStatistic ) ; 
	}

	// Other stats, these might not be used.
	if( statistic_spec == "AA" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new AAStatistic ) ;
	}
	if( statistic_spec == "AB" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new ABStatistic ) ;
	}
	if( statistic_spec == "BB" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new BBStatistic ) ;
	}
	else if( statistic_spec == "HWE(slow)" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new HardyWeinbergExactTestStatistic ) ;
	}
	else if( statistic_spec == "MLIG" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new LogMaximumLikelihoodForIndependentGenotypes ) ;
	}
	else if( statistic_spec == "MLIGHW" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new LogMaximumLikelihoodForIndependentGenotypesInHardyWeinberg ) ;
	}
	else if( statistic_spec == "HWLR" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new HardyWeinbergLikelihoodRatioTestStatistic ) ;
	}
	else if( statistic_spec == "NULL" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new NullStatistic ) ;
	}
	else {
		throw StatisticNotFoundException( "Unable to construct statistic \"" + statistic_spec
			+ "\" -- possible values are \"MAF\", \"HWE\", \"missing-rate\", \"heterozygosity\", \"SNPID\", \"RSID\", \"alleles\", \"ID1\", \"ID2\"." ) ;
	}
}

std::auto_ptr< GenotypeAssayStatistic > GenRowStatisticFactory::create_statistic( std::string statistic_spec ) {
	if( statistic_spec == "position" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new GenRowSNPPosition ) ;
	}
	else if( statistic_spec == "SNPID" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new GenRowSNPID ) ;
	}
	else if( statistic_spec == "RSID" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new GenRowRSID ) ;
	}
	else if( statistic_spec == "minor_allele" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new GenRowAllele( GenRowAllele::minor ) ) ;
	}
	else if( statistic_spec == "major_allele" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new GenRowAllele( GenRowAllele::major ) ) ;
	}
	else {
		return std::auto_ptr< GenotypeAssayStatistic >() ;
	}
}

std::auto_ptr< GenotypeAssayStatistic > SampleRowStatisticFactory::create_statistic( std::string statistic_spec ) {
	statistic_spec = string_utils::strip( statistic_spec ) ;
	if( statistic_spec == "ID1" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new SampleRowID1 ) ;
	}
	else if( statistic_spec == "ID2" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new SampleRowID2 ) ;
	}
	else {
		return std::auto_ptr< GenotypeAssayStatistic >() ;
	}	
}

void SampleRowStatisticFactory::add_statistics( std::vector< std::string > statistic_specs, GenotypeAssayStatistics& statistics ) {
	GenotypeAssayStatisticFactory::add_statistics( statistic_specs, statistics ) ;
}

void GenRowStatisticFactory::add_statistics( std::vector< std::string > statistic_specs, GenotypeAssayStatistics& statistics ) {
	GenotypeAssayStatisticFactory::add_statistics( statistic_specs, statistics ) ;
}

