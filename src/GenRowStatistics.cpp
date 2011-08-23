#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <numeric>

#include "GenRowStatistics.hpp"
#include "floating_point_utils.hpp"


GenRowStatistics::GenRowStatistics()
	: m_row(0)
{
}

void GenRowStatistics::process( GenRow const& row ) {
	m_row = &row ;
	base_t::process( row.begin_genotype_proportions(), row.end_genotype_proportions() ) ;
}

double GenRowSpecificStatistic::calculate_value( GenotypeAssayStatistics const& statistics ) const {
	return calculate_value( get_row_statistics( statistics )) ;
}

std::string GenRowSpecificStatistic::calculate_string_value( GenotypeAssayStatistics const& statistics ) const {
	return calculate_string_value( get_row_statistics( statistics )) ;
}

std::string GenRowSpecificStatistic::calculate_string_value( GenRowStatistics const& row_statistics ) const {
	return base_t::calculate_string_value( row_statistics ) ;
}

GenRowStatistics const& GenRowSpecificStatistic::get_row_statistics( GenotypeAssayStatistics const& statistics ) const {
	try {
		return dynamic_cast< GenRowStatistics const& >( statistics ) ;
	}
	catch (std::bad_cast const& e ) {
		throw GenotypeAssayStatisticException( "Unable to cast statistics to type GenRowStatistic.  This statistic only applies to gen rows" ) ;
	}
}

double GenRowChromosome::calculate_value( GenRowStatistics const& row_statistics ) const {
	throw GenotypeAssayStatisticException( "GenRowChromosome does not support double values." ) ;
}

std::string GenRowChromosome::calculate_string_value( GenRowStatistics const& row_statistics ) const {
	std::ostringstream oStream ;
	oStream << row_statistics.row().chromosome() ;
	return oStream.str() ;
}

double GenRowSNPPosition::calculate_value( GenRowStatistics const& row_statistics ) const {
	return row_statistics.row().SNP_position() ;
}

std::string GenRowSNPPosition::calculate_string_value( GenRowStatistics const& row_statistics ) const {
	std::ostringstream oStream ;
	oStream << std::fixed << std::noshowpoint << row_statistics.row().SNP_position() ;
	return oStream.str() ;
}

double GenRowSNPID::calculate_value( GenRowStatistics const& row_statistics ) const {
	throw GenotypeAssayStatisticException( "GenRowSNPID does not support double values." ) ;
}

std::string GenRowSNPID::calculate_string_value( GenRowStatistics const& row_statistics ) const {
	return row_statistics.row().SNPID() ;
}

double GenRowRSID::calculate_value( GenRowStatistics const& row_statistics ) const {
	throw GenotypeAssayStatisticException( "GenRowRSID does not support double values." ) ;
}

std::string GenRowRSID::calculate_string_value( GenRowStatistics const& row_statistics ) const {
	return row_statistics.row().RSID() ;
}

GenRowAllele::GenRowAllele( AlleleSelector const& selector ) : m_selector( selector ) {}

double GenRowAllele::calculate_value( GenRowStatistics const& row_statistics ) const {
	throw GenotypeAssayStatisticException( "GenRowAllele does not support double values." ) ;	
}

std::string GenRowAllele::calculate_string_value( GenRowStatistics const& row_statistics ) const {
	std::string result ;
	switch( m_selector ) {
		case first_allele:
			result = std::string( std::size_t(1), row_statistics.row().first_allele() ) ;
			break ;
		case second_allele:
			result = std::string( std::size_t(1), row_statistics.row().second_allele() ) ;
			break ;
		case minor_allele:
		case major_allele:
			AlleleProportions allele_proportions = row_statistics.get_mean_allele_proportions() ;
			bool choose_first_allele
				= ((m_selector == minor_allele) && ( allele_proportions.A() <= allele_proportions.B()))
				| ((m_selector == major_allele) && ( allele_proportions.A() >= allele_proportions.B())) ;
			if( choose_first_allele ) {
				result = std::string( std::size_t(1), row_statistics.row().first_allele() ) ;
			}
			else {
				result = std::string( std::size_t(1), row_statistics.row().second_allele() ) ;
			}
			break ;
	}
	return result ;
}


