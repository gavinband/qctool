
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread/thread.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "components/HaplotypeFrequencyComponent/FlatTableDBOutputter.hpp"

// #define DEBUG_FLATTABLEDBOUTPUTTER 1

namespace haplotype_frequency_component {
	FlatTableDBOutputter::UniquePtr FlatTableDBOutputter::create(
		std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata,
		std::string const& snp_match_fields
	) {
		return UniquePtr( new FlatTableDBOutputter( filename, analysis_name, analysis_description, metadata, snp_match_fields ) ) ;
	}

	FlatTableDBOutputter::FlatTableDBOutputter(
		std::string const& filename,
		std::string const& analysis_name,
		std::string const& analysis_description,
		Metadata const& metadata,
		std::string const& snp_match_fields
	):
		m_outputter( filename, analysis_name, analysis_description, metadata, boost::optional< db::Connection::RowId >(), snp_match_fields ),
		m_table_name( "ld_analysis" + genfile::string_utils::to_string( m_outputter.analysis_id() ) ),
		m_max_variants_per_block( 1000 )
	{}

	FlatTableDBOutputter::~FlatTableDBOutputter() {
	}
	
	void FlatTableDBOutputter::set_table_name( std::string const& table_name ) {
		if( table_name.find_first_not_of( "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_" ) != std::string::npos ) {
			throw genfile::BadArgumentError( "qcdb::FlatTableDBOutputter::set_table_name()", "table_name=\"" + table_name + "\"" ) ;
		}
		m_table_name = table_name ;
	}

	FlatTableDBOutputter::AnalysisId FlatTableDBOutputter::analysis_id() const {
		return m_outputter.analysis_id() ;
	}

	void FlatTableDBOutputter::finalise( long options ) {
		store_block() ;
		m_variants.clear() ;
		m_values.clear() ;

		if( options & qcdb::eCreateIndices ) {
			db::Connection::ScopedTransactionPtr transaction = m_outputter.connection().open_transaction( 7200 ) ;
			m_outputter.connection().run_statement( "CREATE INDEX IF NOT EXISTS " + m_table_name + "_index ON " + m_table_name + "( variant1_id, variant2_id )" ) ;
		}

		m_outputter.finalise( options ) ;
	}

	void FlatTableDBOutputter::add_variable( std::string const& variable ) {
		VariableMap::left_const_iterator where = m_variables.left.find( variable ) ;
		if( where == m_variables.left.end() ) {
			if( m_insert_data_sql.get() ) {
				// Uh-oh, table columns are already fixed.
				// TODO: alter to put this data in SummaryData table?
				throw genfile::BadArgumentError( "qcdb::FlatTableDBOutputter::add_variable()", "variable=\"" + variable + "\"" ) ;
			}
			else {
				// Still have time to add the variable to our list of variables, retaining the order of addition.
				m_variables.left.insert( VariableMap::left_value_type( variable, m_variables.size() ) ) ;
			}
		}
	}

	void FlatTableDBOutputter::create_new_variant_pair(
		genfile::VariantIdentifyingData const& variant1,
		genfile::VariantIdentifyingData const& variant2
	 ) {
		if( m_variants.size() == m_max_variants_per_block ) {
			store_block() ;
			m_variants.clear() ;
			m_values.clear() ;
		}
		m_variants.push_back( std::make_pair( variant1, variant2 ) ) ;
	}

	void FlatTableDBOutputter::store_per_variant_pair_data(
		genfile::VariantIdentifyingData const& variant1,
		genfile::VariantIdentifyingData const& variant2,
		std::string const& variable,
		genfile::VariantEntry const& value
	) {
		std::pair< genfile::VariantIdentifyingData, genfile::VariantIdentifyingData > variant_pair = std::make_pair( variant1, variant2 ) ;
		bool const new_variant = m_variants.empty() || variant_pair != m_variants.back() ;
		if( new_variant ) {
			create_new_variant_pair( variant1, variant2 ) ;
		}

		VariableMap::left_const_iterator where = m_variables.left.find( variable ) ;
		if( where == m_variables.left.end() ) {
			if( m_insert_data_sql.get() ) {
				// Uh-oh, table columns are already fixed.
				// TODO: alter to put this data in SummaryData table?
				throw genfile::BadArgumentError( "qcdb::FlatTableDBOutputter::store_per_variant_pair_data()", "variable=\"" + variable + "\"" ) ;
			}
			else {
				// Still have time to add the variable to our list of variables, retaining the order of addition.
				where = m_variables.left.insert( VariableMap::left_value_type( variable, m_variables.size() ) ).first ;
			}
		}

		// Store the value of this variable
		m_values[ std::make_pair( m_variants.size() - 1, where->second ) ] = value ;
	}

	void FlatTableDBOutputter::store_block() {
		db::Connection::ScopedTransactionPtr transaction = m_outputter.connection().open_transaction( 7200 ) ;

		if( !m_insert_data_sql.get() ) {
			create_schema() ;
			// create_variables() ;
		}
		for( std::size_t i = 0; i < m_variants.size(); ++i ) {
			db::Connection::RowId const variant1_id = m_outputter.get_or_create_variant( m_variants[i].first ) ;
			db::Connection::RowId const variant2_id = m_outputter.get_or_create_variant( m_variants[i].second ) ;
			store_data_for_variants( i, m_outputter.analysis_id(), variant1_id, variant2_id ) ;
		}
	}

	std::string FlatTableDBOutputter::get_table_name() const {
		return m_table_name ;
	}

	void FlatTableDBOutputter::create_schema() {
		using genfile::string_utils::to_string ;
		std::ostringstream schema_sql ;
		std::ostringstream index_sql ;
		std::ostringstream insert_data_sql ;
		std::ostringstream insert_data_sql_columns ;
		std::ostringstream insert_data_sql_values ;
		std::string const& table_name = m_table_name ;
		schema_sql << "CREATE TABLE IF NOT EXISTS "
			<< table_name
			<< " ( "
			"analysis_id INT NOT NULL REFERENCES Entity( id ), "
			"variant1_id INT NOT NULL REFERENCES Variant( id ), "
			"variant2_id INT NOT NULL REFERENCES Variant( id )"
		;
		
		insert_data_sql << "INSERT INTO "
			<< table_name ;
		insert_data_sql_columns << "( analysis_id, variant1_id, variant2_id" ;
		insert_data_sql_values << "VALUES( ?1, ?2, ?3" ;
		
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;

		for( std::size_t bind_i = 4; var_i != end_var_i; ++var_i, ++bind_i ) {
			schema_sql << ", " ;
			insert_data_sql_columns << ", " ;
			insert_data_sql_values << ", " ;
			schema_sql
				<< '"'
				<< var_i->second
				<< '"'
				<< " NULL" ;
				
			insert_data_sql_columns << '"' << var_i->second << '"' ;
			insert_data_sql_values << "?" << to_string( bind_i ) ;
		}

		schema_sql << " ) ; " ;

		insert_data_sql_columns << ") " ;
		insert_data_sql_values << ") ; " ;

#if DEBUG_FLATTABLEDBOUTPUTTER
		std::cerr << "Creating table " << table_name << " using this SQL:\n"
			<< schema_sql.str()
			<< "\n" ;
#endif		
		m_outputter.connection().run_statement( schema_sql.str() ) ;

		std::ostringstream view_sql ;	
		view_sql
			<< "CREATE VIEW IF NOT EXISTS `" << table_name << "View` AS "
			<< "SELECT V1.rsid AS variant1_rsid, V1.chromosome AS variant1_chromosome, V1.position AS variant1_position, "
			<< "V2.rsid AS variant2_rsid, V2.chromosome AS variant2_chromosome, V2.position AS variant2_position, A.name AS analysis, "
			<< "T.* FROM `"
			<< table_name << "` T "
			<< "INNER JOIN Variant V1 ON V1.id = T.variant1_id "
			<< "INNER JOIN Variant V2 ON V2.id = T.variant2_id "
			<< "INNER JOIN Analysis A ON A.id = T.analysis_id"
		;

#if DEBUG_FLATTABLEDBOUTPUTTER
		std::cerr << "Creating view using SQL:\n" << view_sql.str() << "\n" ;
#endif
		m_outputter.connection().run_statement( view_sql.str() ) ;

		insert_data_sql << " " << insert_data_sql_columns.str() << " " << insert_data_sql_values.str() ;
#if DEBUG_FLATTABLEDBOUTPUTTER
		std::cerr << "Inserts will use this SQL:\n"
			<< insert_data_sql.str()
			<< "\n" ;
#endif
		
		m_insert_data_sql = m_outputter.connection().get_statement(
			insert_data_sql.str()
		) ;
	}

	/*
	void FlatTableDBOutputter::create_variables() {
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;
		for( ; var_i != end_var_i; ++var_i ) {
			m_outputter.get_or_create_column( m_analysis_id,, var_i->second, variable_class_id ) ;
		}
	}
	*/
	
	void FlatTableDBOutputter::store_data_for_variants(
		std::size_t const variant_i,
		db::Connection::RowId const analysis_id,
		db::Connection::RowId const variant1_id,
		db::Connection::RowId const variant2_id
	) {
		m_insert_data_sql->bind( 1, analysis_id ) ;
		m_insert_data_sql->bind( 2, variant1_id ) ;
		m_insert_data_sql->bind( 3, variant2_id ) ;
		
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;
		
		for( std::size_t bind_i = 4; var_i != end_var_i; ++var_i, ++bind_i ) {
			ValueMap::const_iterator where = m_values.find( std::make_pair( variant_i, var_i->first )) ;
			if( where == m_values.end() ) {
				m_insert_data_sql->bind( bind_i, genfile::MissingValue() ) ;
			} else {
				m_insert_data_sql->bind( bind_i, where->second ) ;
			}
		}
		
		m_insert_data_sql->step() ;
		m_insert_data_sql->reset() ;
	}
}
