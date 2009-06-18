#include <cassert>
#include <iostream>
#include <iomanip>
#include <numeric>

#include "GenotypeAssayBasicStatistics.hpp"
#include "floating_point_utils.hpp"

void GenotypeAssayBasicStatistics::floor_genotype_amounts() {
	m_genotype_amounts.floor() ;
}

void GenotypeAssayBasicStatistics::round_genotype_amounts() {
	m_genotype_amounts.round() ;
}

std::ostream& GenotypeAssayBasicStatistics::format_column_headers( std::ostream& aStream ) const {
	return aStream << "Total    AA       AB       BB       ??        " ;
}

std::ostream& GenotypeAssayBasicStatistics::format_statistic_values( std::ostream& aStream ) const {
	aStream 
		<< std::setw(7) << std::left << number_of_samples()
		<< "  "
		<< std::setw(7) << std::left << get_genotype_amounts().AA()
		<< "  "
		<< std::setw(7) << std::left << get_genotype_amounts().AB()
		<< "  "
		<< std::setw(7) << std::left << get_genotype_amounts().BB()
		<< "  "
		<< std::setw(8) << std::left << std::noshowpoint << std::fixed << std::setprecision( 2 ) << amount_of_missing_data()
		<< "  " ;

		return aStream ;
}

std::ostream& operator<<( std::ostream& aStream, GenotypeAssayBasicStatistics const& statistics ) {
	return statistics.format_statistic_values( aStream ) ;
}
