#include "GenotypeAssayStatisticFactory.hpp"
#include "GenotypeAssayStatistics.hpp"
#include "GenRowStatistics.hpp"
#include "SNPHWE.hpp"
#include "HardyWeinbergExactTestStatistic.hpp"
#include "LikelihoodRatioTestStatistic.hpp"	
#include "string_utils.hpp"

std::auto_ptr< GenotypeAssayStatistic > GenotypeAssayStatisticFactory::create_statistic( std::string statistic_spec ) {
	statistic_spec = strip( statistic_spec ) ;
	if( statistic_spec == "SNP_position" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new GenRowSNPPosition ) ;
	}
	if( statistic_spec == "SNPID" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new GenRowSNPID ) ;
	}
	if( statistic_spec == "RSID" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new GenRowRSID ) ;
	}
	else if( statistic_spec == "#Samples" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new NumberOfSamplesStatistic ) ;
	}
	else if( statistic_spec == "AA" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new AAStatistic ) ;
	}
	else if( statistic_spec == "AB" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new ABStatistic ) ;
	}
	else if( statistic_spec == "BB" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new BBStatistic ) ;
	}
	else if( statistic_spec == "Missing" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new MissingDataStatistic ) ;
	}
	else if( statistic_spec == "%Missing" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new MissingDataProportionStatistic ) ;
	}
	else if( statistic_spec == "HWE" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new HardyWeinbergExactTestStatistic ) ;
	}
	else if( statistic_spec == "Wigginton" ) {
		return std::auto_ptr< GenotypeAssayStatistic >( new SNPHWEStatistic ) ;
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
			+ "\" -- possible values are \"SNPID\", \"RSID\", \"HWE\", \"Wigginton\", \"MLIG\", \"MLIGHW\", \"HWLR\", \"HWLR test stat\"." ) ;
	}
}

