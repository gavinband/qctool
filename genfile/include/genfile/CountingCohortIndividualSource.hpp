#ifndef GENFILE_COUNTING_COHORT_INDIVIDUAL_SOURCE_HPP
#define GENFILE_COUNTING_COHORT_INDIVIDUAL_SOURCE_HPP

#include "genfile/CohortIndividualSource.hpp"

namespace genfile {
	struct CountingCohortIndividualSource: public CohortIndividualSource
	{
		CountingCohortIndividualSource( std::size_t number_of_samples, std::string const& sample_name_template ) ;
		CountingCohortIndividualSource( std::vector< std::string > const& sample_ids ) ;
		std::size_t get_number_of_individuals() const ;
		ColumnSpec get_column_spec() const ;
		bool check_for_column( std::string const& column_name ) const ;
		Entry get_entry( std::size_t sample_i, std::string const& column_name ) const ;
		std::string get_source_spec() const ;
	private:
		std::vector< std::string > m_sample_ids ;
	} ;
}
#endif

