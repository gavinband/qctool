
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <memory>
#include <boost/bimap.hpp>
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "qcdb/Storage.hpp"
#include "qcdb/FlatFileMultiVariantOutputter.hpp"

namespace qcdb {
	namespace {
		void append_to_string( std::string* target, std::string const& value ) {
			(*target) += ( target->size() > 0 ? "," : "" ) + value ;
		}

		std::vector< std::string > generate_key_entry_names( std::size_t n ) {
			std::vector< std::string > result(n) ;
			for( std::size_t i = 0; i < n; ++i ) {
				result[i] = (boost::format( "variant%d" ) % (i+1)).str() ;
			}
			return result ;
		}
	}

	FlatFileMultiVariantOutputter::UniquePtr FlatFileMultiVariantOutputter::create(
		std::string const& filename,
		std::size_t const number_of_key_entries,
		std::string const& analysis_name,
		Metadata const& metadata
	) {
		return UniquePtr( new FlatFileMultiVariantOutputter( filename, number_of_key_entries, analysis_name, metadata ) ) ;
	}

	FlatFileMultiVariantOutputter::FlatFileMultiVariantOutputter(
		std::string const& filename,
		std::size_t const number_of_key_entries,
		std::string const& analysis_name,
		Metadata const& metadata
	):
		m_filename( filename ),
		m_key_entry_names( generate_key_entry_names(number_of_key_entries) ),
		m_analysis_name( analysis_name ),
		m_metadata( metadata ),
		m_max_keys_per_block( 1000 )
	{
		m_keys.reserve( 1000 ) ;
	}
	
	FlatFileMultiVariantOutputter::~FlatFileMultiVariantOutputter() {
	}
	
	void FlatFileMultiVariantOutputter::finalise( long ) {
		store_block() ;
		m_keys.clear() ;
		m_values.clear() ;
		m_sink->write_comment( "Completed successfully at " + appcontext::get_current_time_as_string() ) ;
	}

	FlatFileMultiVariantOutputter::AnalysisId FlatFileMultiVariantOutputter::analysis_id() const {
		// A flat file only ever has one analysis.
		return 0 ;
	}

	void FlatFileMultiVariantOutputter::set_variant_names( std::vector< std::string > const& names ) {
		// Ensure we don't change names if we have already written data
		if( m_sink.get() && names != m_key_entry_names ) {
			throw genfile::BadArgumentError( "qcdb::FlatFileMultiVariantOutputter::set_variant_names()", "names" ) ;
		}
		m_key_entry_names = names ;
	}

	void FlatFileMultiVariantOutputter::add_variable(
		std::string const& variable
	) {
		VariableMap::left_const_iterator where = m_variables.left.find( variable ) ;
		if( where == m_variables.left.end() ) {
			if( m_sink.get() ) {
				// Uh-oh, have already written a header.
				throw genfile::BadArgumentError( "qcdb::FlatFileMultiVariantOutputter::add_variable()", "variable=\"" + variable + "\"" ) ;
			}
			else {
				// Still have time to add the variable to our list of variables, retaining the order of addition.
				where = m_variables.left.insert( VariableMap::left_value_type( variable, m_variables.size() ) ).first ;
			}
		}
	}

	void FlatFileMultiVariantOutputter::create_new_key( MultiVariantStorage::Key const& key ) {
		if( m_keys.size() == m_max_keys_per_block ) {
			store_block() ;
			m_keys.clear() ;
			m_values.clear() ;
		}
		m_keys.push_back( key ) ;
	}

	void FlatFileMultiVariantOutputter::store_data_for_key(
		MultiVariantStorage::Key const& key,
		std::string const& variable,
		genfile::VariantEntry const& value
	) {
		assert( key.size() <= m_key_entry_names.size() ) ;
		bool const new_key = m_keys.empty() || key != m_keys.back() ;
		if( new_key ) {
			// If we have a whole block's worth of data, store it now.
			if( m_keys.size() == m_max_keys_per_block ) {
				store_block() ;
				m_keys.clear() ;
				m_values.clear() ;
			}
			m_keys.push_back( key ) ;
		}

		VariableMap::left_const_iterator where = m_variables.left.find( variable ) ;
		if( where == m_variables.left.end() ) {
			if( m_sink.get() ) {
				// Uh-oh, have already written a header.
				throw genfile::BadArgumentError( "qcdb::FlatFileMultiVariantOutputter::store_per_variant_data()", "variable=\"" + variable + "\"" ) ;
			}
			else {
				// Still have time to add the variable to our list of variables, retaining the order of addition.
				where = m_variables.left.insert( VariableMap::left_value_type( variable, m_variables.size() ) ).first ;
			}
		}

		// Store the value of this variable.
		m_values[ std::make_pair( m_keys.size() - 1, where->second ) ] = value ;
	}
	
	void FlatFileMultiVariantOutputter::store_block() {
		if( !m_sink.get() ) {
			m_sink = statfile::BuiltInTypeStatSink::open( m_filename ) ;
			m_sink->write_metadata( format_metadata() ) ;
			create_columns() ;
		}
		for( std::size_t key_i = 0; key_i < m_keys.size(); ++key_i ) {
			// output key
			Key const& key = m_keys[ key_i ] ;
			for( std::size_t i = 0; i < key.size(); ++i ) {
				genfile::VariantIdentifyingData const& variant = key[i] ;
				std::string SNPID ;
				if( variant.number_of_identifiers() == 1 ) {
					SNPID = "NA" ;
				} else {
					variant.get_identifiers( boost::bind( &append_to_string, &SNPID, _1 ), 1 ) ;
				}
				(*m_sink) << SNPID << variant.get_primary_id() << variant.get_position().chromosome()
				<< variant.get_position().position()
				<< variant.get_allele(0)
				<< (( variant.number_of_alleles() < 2 ) ? "." : variant.get_alleles_as_string( ",", 1, variant.number_of_alleles() )) ;
			}
			genfile::MissingValue const NA ;
			for( std::size_t i = key.size(); i < m_key_entry_names.size(); ++i ) {
				(*m_sink)
					<< NA << NA << NA
					<< NA << NA << NA ;
			}

			// output variables
			VariableMap::right_const_iterator
				var_i = m_variables.right.begin(),
				end_var_i = m_variables.right.end() ;
			for( ; var_i != end_var_i; ++var_i ) {
				ValueMap::const_iterator where = m_values.find( std::make_pair( key_i, var_i->first )) ;
				if( where == m_values.end() ) {
					(*m_sink) << genfile::MissingValue() ;
				}
				else {
					(*m_sink) << where->second ;
				}
			}
			(*m_sink) << statfile::end_row() ;
		}
	}
	
	std::string FlatFileMultiVariantOutputter::format_metadata() const {
		std::ostringstream str ;
		str << "Analysis: \"" << m_analysis_name << "\"\n"
			<< " started: " << appcontext::get_current_time_as_string() << "\n" ;
		str << "\nAnalysis properties:\n" ;
		for( Metadata::const_iterator i = m_metadata.begin(); i != m_metadata.end(); ++i ) {
			str << "  "
				<< i->first
				<< " "
				<< genfile::string_utils::join( i->second.first, " " )
				<< " (" + i->second.second + ")\n" ;
		}
		return str.str() ;
	}
	
	void FlatFileMultiVariantOutputter::create_columns() {
		using genfile::string_utils::to_string ;
		for( std::size_t i = 0; i < m_key_entry_names.size(); ++i ) {
			std::string const& stub = m_key_entry_names[i] ;
			(*m_sink) | (stub+":alternate_ids") | (stub+":rsid") | (stub+":chromosome") | (stub+":position") | (stub+":alleleA") | (stub+":alleleB") ;
		}
		VariableMap::right_const_iterator
			i = m_variables.right.begin(),
			end_i = m_variables.right.end() ;
		for( ; i != end_i; ++i ) {
			(*m_sink).add_column( i->second ) ;
		}
		(*m_sink) << statfile::begin_data() ;
	}
}
