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
	return aStream ;
}

std::ostream& GenotypeAssayBasicStatistics::format_statistic_values( std::ostream& aStream ) const {
		return aStream ;
}

std::ostream& operator<<( std::ostream& aStream, GenotypeAssayBasicStatistics const& statistics ) {
	return statistics.format_statistic_values( aStream ) ;
}
