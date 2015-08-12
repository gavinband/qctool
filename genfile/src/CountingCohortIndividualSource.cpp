#include <string>
#include <cassert>
#include <boost/format.hpp>
#include "genfile/string_utils/string_utils.hpp"
#include "genfile/CountingCohortIndividualSource.hpp"

namespace genfile {
	CountingCohortIndividualSource::CountingCohortIndividualSource(
		std::size_t number_of_samples,
		std::string const& sample_name_template
	):
	 	m_sample_ids( number_of_samples )
	{
		assert( sample_name_template.find( "%d" ) != std::string::npos ) ;
		for( std::size_t i = 0; i < m_sample_ids.size(); ++i ) {
			m_sample_ids[i] = ( boost::format( sample_name_template ) % i ).str() ;
		}
	}

	CountingCohortIndividualSource::CountingCohortIndividualSource( std::vector< std::string > const& sample_ids ):
		m_sample_ids( sample_ids )
	{}

	std::size_t CountingCohortIndividualSource::get_number_of_individuals() const {
		return m_sample_ids.size() ;
	}

	CohortIndividualSource::ColumnSpec CountingCohortIndividualSource::get_column_spec() const {
		ColumnSpec result ;
		result.add_column( "ID_1", e_ID_COLUMN ) ;
		result.add_column( "ID_2", e_ID_COLUMN ) ;
		result.add_column( "missing", e_MISSINGNESS_COLUMN ) ;
		return result ;
	}

	bool CountingCohortIndividualSource::check_for_column( std::string const& column_name ) const {
		std::string const lower_column_name = string_utils::to_lower( column_name ) ;
		return column_name == "id_1" || column_name == "id_2" || column_name == "missing" ;
	}

	CohortIndividualSource::Entry CountingCohortIndividualSource::get_entry( std::size_t sample_i, std::string const& column_name ) const {
		Entry result ;
		std::string const lower_column_name = string_utils::to_lower( column_name ) ;
		if( lower_column_name == "id_1" || lower_column_name == "id_2" ) {
			result = m_sample_ids[sample_i] ;
		} else if( lower_column_name == "missing" ) {
			result = 0.0 ;
		} else {
			assert(0) ;
		}
		return result ;
	}

	std::string CountingCohortIndividualSource::get_source_spec() const {
		return "CountingCohortIndividualSource(" + string_utils::to_string( m_sample_ids.size() ) + ")" ;
	}
}

