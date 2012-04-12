
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <sstream>
#include <iostream>
#include <iomanip>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/CrossCohortCovariateValueMapping.hpp"
#include "genfile/ValueMappingCohortIndividualSource.hpp"

namespace genfile {
	ValueMappingCohortIndividualSource::ValueMappingCohortIndividualSource(
		CohortIndividualSource::UniquePtr source
	):
		m_source( source )
	{
		assert( m_source.get() ) ;
	}
	
	ValueMappingCohortIndividualSource::~ValueMappingCohortIndividualSource() {
		Mappings::const_iterator
			i = m_mappings.begin(),
			end_i = m_mappings.end() ;
		for( ; i != end_i; ++i ) {
			delete i->second.second ;
		}
	}
		
	void ValueMappingCohortIndividualSource::add_mapping(
		std::string const& source_column_name,
		std::string const& target_column_name,
		CrossCohortCovariateValueMapping::UniquePtr mapping
	) {
		if( !m_source->check_for_column( source_column_name )) {
			throw BadArgumentError( "genfile::ValueMappingCohortIndividualSource::add_mapping()", "source_column_name = \"" + source_column_name + "\"" ) ;
		}
		Mappings::const_iterator i = m_mappings.find( target_column_name ) ;
		if( i != m_mappings.end() ){
			throw DuplicateKeyError( "Mapped columns of " + get_source_spec(), target_column_name ) ;
		}
		assert( mapping.get() ) ;
		m_mappings[ target_column_name ] = std::make_pair( source_column_name, mapping.release() ) ;
	}
		
	CohortIndividualSource::Entry ValueMappingCohortIndividualSource::get_entry( std::size_t sample_i, std::string const& column_name ) const {
		Entry result ;
		Mappings::const_iterator i = m_mappings.find( column_name ) ;
		if( i != m_mappings.end() ) {
			std::string source_column_name = i->second.first ;
			result = m_source->get_entry( sample_i, source_column_name ) ;
			if( !result.is_missing() && i != m_mappings.end() ) {
				result = i->second.second->get_mapped_value( result ) ;
			}
		}
		else {
			result = m_source->get_entry( sample_i, column_name ) ;
		}
		return result ;
	}

	std::string ValueMappingCohortIndividualSource::get_source_spec() const {
		return "mapped:" + m_source->get_source_spec() ;
	}
	
	std::string ValueMappingCohortIndividualSource::get_summary( std::string const& prefix ) const {
		std::ostringstream ostr ;
		ostr << "Mapped values:\n" ;
		Mappings::const_iterator
			i = m_mappings.begin(),
			end_i = m_mappings.end() ;
		std::size_t max_column_name_length = 0 ;
		for( ; i != end_i; ++i ) {
			max_column_name_length = std::max( max_column_name_length, i->first.size() ) ;
		}
		for( i = m_mappings.begin(); i != end_i; ++i ) {
			ostr << prefix << std::setw( max_column_name_length ) << i->first << ": "
				<< i->second.second->get_summary( prefix + std::string( max_column_name_length + 2, ' ' ))
				<< "\n\n" ;
		}
		return ostr.str() ;
	}
}
