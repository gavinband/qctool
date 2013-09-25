
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <string>
#include <iosfwd>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/DataLinkingCohortIndividualSource.hpp"

namespace genfile {
	// class DataLinkingCohortIndividualSource
	// This class implements a source which links data using the ID_1 columns
	// between several sources passed to it.
	// The first source, which must be passed to the constructor, is regarded as the 'master' source,
	// and contains the sample IDs which the source representes.
	// Subsequent sources may be passed to the add_source() method and can contain any sample IDs at all
	// (but sample IDs not in the first source will be ignored.)
	// It is legal to have two columns of the same name, but only the first will ever be matched by name.
	DataLinkingCohortIndividualSource::DataLinkingCohortIndividualSource( CohortIndividualSource::UniquePtr master )
	{
		m_sources.push_back( master ) ;
		for( std::size_t i = 0; i < m_sources[0].get_number_of_individuals(); ++i ) {
			m_sample_indices[ i ].push_back( i ) ;
		}
		m_column_spec = m_sources[0].get_column_spec() ;
	}

	DataLinkingCohortIndividualSource::~DataLinkingCohortIndividualSource()
	{}

	std::size_t DataLinkingCohortIndividualSource::get_number_of_individuals() const {
		return m_sources[0].get_number_of_individuals() ;
	}

	void DataLinkingCohortIndividualSource::add_source( CohortIndividualSource::UniquePtr source ) {
		m_sources.push_back( source ) ;
		// link samples
		for( std::size_t i = 0; i < m_sources[0].get_number_of_individuals(); ++i ) {
			std::string const ID = m_sources[0].get_entry( i, "ID_1" ).as< std::string >() ;
			std::vector< std::size_t > samples = m_sources.back().find_samples_by_value( "ID_1", ID ) ;
			assert( samples.size() <= 1 ) ;
			boost::optional< std::size_t > index ;
			if( samples.size() == 1 ) {
				m_sample_indices[i].push_back( samples[0] ) ;
			} else {
				m_sample_indices[i].push_back( boost::optional< std::size_t >() ) ;
			}
		}
		
		// link columns
		CohortIndividualSource::ColumnSpec spec = m_column_spec ;
		CohortIndividualSource::ColumnSpec source_spec = m_sources.back().get_column_spec() ;
		for( std::size_t i = 0; i < source_spec.size(); ++i ) {
			std::string const name = source_spec[i].name() ;
			if( !spec.check_for_column( name ) ) {
				spec.add_column( name, source_spec[i].type() ) ;
			} else if( source_spec[i].type() != spec[ name ].type() ) {
				throw DuplicateEntryError( "DataLinkingCohortIndividualSource::add_source()", get_source_spec(), "column types", name, "Conflicting type for column \"" + name + "\" in source \"" + m_sources.back().get_source_spec() + "\"" ) ; 
			}
		}
		m_column_spec = spec ;
	}

	CohortIndividualSource::ColumnSpec DataLinkingCohortIndividualSource::get_column_spec() const {
		return m_column_spec ;
	}

	CohortIndividualSource::Entry DataLinkingCohortIndividualSource::get_entry( std::size_t sample_i, std::string const& column_name ) const {
		using genfile::string_utils::to_string ;
		
		Entry result ;
		bool column_found = false ;
		bool entry_found = false ;

		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			ColumnSpec const spec = m_sources[i].get_column_spec() ;
			if( spec.check_for_column( column_name )) {
				boost::optional< std::size_t > sample_index = sample_i ;
				if( i > 0 ) {
					IndexMap::const_iterator where = m_sample_indices.find( sample_i ) ;
					assert( where != m_sample_indices.end() ) ;
					sample_index = where->second[i] ;
				}
				
				if( sample_index ) {
					Entry const source_entry = m_sources[i].get_entry( *sample_index, column_name ) ;
					if( entry_found ) {
						if( source_entry != result ) {
							std::string const ID = m_sources[0].get_entry( sample_i, "ID_1" ).as< std::string >() ;
							
							throw genfile::DuplicateEntryError(
								"DataLinkingCohortIndividualSource::get_entry()",
								get_source_spec(),
								ID,
								column_name,
								"Value " + to_string( source_entry ) + " for sample \"" + ID + "\", variable \"" + column_name + "\" conflicts with previous value " + to_string( result )
							) ;
						}
					} else {
						result = source_entry ;
						entry_found = true ;
					}
				}

				column_found = true ;
			}
		}
		assert( column_found ) ;
		return result ;
	}

	CohortIndividualSource const& DataLinkingCohortIndividualSource::get_base_source() const {
		return m_sources[0].get_base_source() ;
	}
	CohortIndividualSource const& DataLinkingCohortIndividualSource::get_parent_source() const {
		return m_sources[0] ;
	}

	std::string DataLinkingCohortIndividualSource::get_source_spec() const {
		std::ostringstream str ;
		str << "DataLinkingCohortIndividualSource( " ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			if( i > 0 ) {
				str << ", " ;
			}
			str << m_sources[i].get_source_spec() ;
		}
		str << " )" ;
		return str.str() ;
	}
}

