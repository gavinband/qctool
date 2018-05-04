
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include "metro/SampleRange.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/string_utils/string_utils.hpp"
#include "metro/SampleDataInterpretation.hpp"

namespace metro {
	SampleDataInterpretation::SampleDataInterpretation( genfile::CohortIndividualSource const& samples, std::vector< std::string > column_names ):
		m_samples( samples ),
		m_column_spec( m_samples.get_column_spec() ),
		m_columns( column_names )
	{
		for( std::size_t i = 0; i < m_columns.size(); ++i ) {
			m_mappings.push_back( 
				genfile::CrossCohortCovariateValueMapping::create( m_column_spec[i] )
			) ;
			m_mappings.back().add_source( samples ) ;
		}
	}
	
	namespace {
		/*
		* Compute the layout of columns in the result matrix
		* layout means a vector v of integers of length d+1 (if there are d columns to lay out), such that
		* the matrix columns [ v[i], v(i+1) ] contain the expanded data.
		*/
		void get_layout(
			genfile::CohortIndividualSource::ColumnSpec const& column_spec,
			std::vector< std::string > const& column_names,
			boost::ptr_vector< genfile::CrossCohortCovariateValueMapping > const& mappings,
			std::vector< int >* layout,
			std::vector< std::string >* result_column_names
		) {
			assert( layout != 0 ) ;
			assert( result_column_names != 0 ) ;
			assert( layout->empty() ) ;
			assert( result_column_names->empty() ) ;
			layout->push_back( 0 ) ;
	
			for( std::size_t i = 0; i < column_names.size(); ++i ) {
				std::string const& column_name = column_names[i] ;
				genfile::CrossCohortCovariateValueMapping const& mapping = mappings[i] ;
				if( column_spec[ i ].is_discrete() ) {
					int const number_of_levels = mapping.get_number_of_distinct_mapped_values() ;
					if( number_of_levels > 1 ) {
						layout->push_back( layout->back() + number_of_levels - 1 ) ;
						for(
							int i = 1; // we need only N-1 levels as covariates, as we also have a baseline term.
							i < mapping.get_number_of_distinct_mapped_values() ;
							++i
						) {
							result_column_names->push_back(
								column_names[i] + "=" + genfile::string_utils::to_string( mapping.get_unmapped_value( i + 1 ))
							) ;
						}
					}
				} else {
					layout->push_back( layout->back() + 1 ) ;
					result_column_names->push_back( column_names[i] ) ;
				}
			}
		}
	}

	SampleDataInterpretation::Matrix SampleDataInterpretation::compute_result() const {
		
	}
	std::vector< std::string > const& SampleDataInterpretation::column_names() const {
		
	}
	SampleDataInterpretation::SampleRanges SampleDataInterpretation::compute_nonmissing_samples() const {
		
	}
}

