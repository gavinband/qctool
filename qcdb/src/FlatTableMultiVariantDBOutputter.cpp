
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread/thread.hpp>
#include <boost/format.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/db/Connection.hpp"
#include "genfile/db/SQLStatement.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "qcdb/FlatTableMultiVariantDBOutputter.hpp"

// #define DEBUG_FLATTABLEDBOUTPUTTER 1

namespace qcdb {
	namespace {
		std::vector< std::string > generate_key_entry_names( std::size_t n ) {
			std::vector< std::string > result(n) ;
			for( std::size_t i = 0; i < n; ++i ) {
				result[i] = (boost::format( "variant%d" ) % (i+1)).str() ;
			}
			return result ;
		}
	}

	FlatTableMultiVariantDBOutputter::UniquePtr FlatTableMultiVariantDBOutputter::create(
		std::string const& filename,
		std::size_t number_of_key_fields,
		std::string const& analysis_name,
		std::string const& analysis_description,
		Metadata const& metadata,
		std::string const& snp_match_fields,
		boost::optional< genfile::db::Connection::RowId > analysis_id
	) {
		return UniquePtr(
			new FlatTableMultiVariantDBOutputter( filename, number_of_key_fields, analysis_name, analysis_description, metadata, snp_match_fields, analysis_id )
		) ;
	}

	FlatTableMultiVariantDBOutputter::FlatTableMultiVariantDBOutputter(
		std::string const& filename,
		std::size_t number_of_key_fields,
		std::string const& analysis_name,
		std::string const& analysis_description,
		Metadata const& metadata,
		std::string const& snp_match_fields,
		boost::optional< genfile::db::Connection::RowId > analysis_id
	):
		m_outputter( filename, analysis_name, analysis_description, metadata, analysis_id, snp_match_fields ),
		m_key_entry_names( generate_key_entry_names( number_of_key_fields )),
		m_table_name( "analysis" + genfile::string_utils::to_string( m_outputter.analysis_id() ) ),
		m_without_rowid( false ),
		m_max_variants_per_block( 1000 )
	{}

	FlatTableMultiVariantDBOutputter::~FlatTableMultiVariantDBOutputter() {
	}
	
	void FlatTableMultiVariantDBOutputter::set_table_name( std::string const& table_name ) {
		if( table_name.find_first_not_of( "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_" ) != std::string::npos ) {
			throw genfile::BadArgumentError( "qcdb::FlatTableMultiVariantDBOutputter::set_table_name()", "table_name=\"" + table_name + "\"" ) ;
		}
		m_table_name = table_name ;
	}

	void FlatTableMultiVariantDBOutputter::set_without_rowid() {
		m_without_rowid = true ;
	}

	FlatTableMultiVariantDBOutputter::AnalysisId FlatTableMultiVariantDBOutputter::analysis_id() const {
		return m_outputter.analysis_id() ;
	}

	void FlatTableMultiVariantDBOutputter::finalise( long options ) {
		store_block() ;
		m_keys.clear() ;
		m_values.clear() ;

		if( options & qcdb::eCreateIndices && !m_without_rowid ) {
			genfile::db::Connection::ScopedTransactionPtr transaction = m_outputter.connection().open_transaction( 7200 ) ;
			std::string sql = "CREATE INDEX IF NOT EXISTS " + m_table_name + "_index ON `" + m_table_name + "` (";
			for( std::size_t i = 0; i < m_key_entry_names.size(); ++i ) {
				sql += (i>0 ? ", `" : " `" ) + (m_key_entry_names[i] + "_id`" ) ;
			}
			sql += ")" ;
			m_outputter.connection().run_statement( sql ) ;
		}

		m_outputter.finalise( options ) ;
	}

	void FlatTableMultiVariantDBOutputter::set_variant_names( std::vector< std::string > const& names ) {
		// Ensure we don't change names once we have begun writing data
		if( names != m_key_entry_names && m_insert_data_sql.get() ) {
			throw genfile::BadArgumentError(
				"qcdb::FlatTableMultiVariantDBOutputter::set_variant_names()",
				"names",
				"Can't change table column names because data has already been written."
			) ;
		}
		m_key_entry_names = names ;
	}

	void FlatTableMultiVariantDBOutputter::add_variable( std::string const& variable ) {
		VariableMap::left_const_iterator where = m_variables.left.find( variable ) ;
		if( where == m_variables.left.end() ) {
			if( m_insert_data_sql.get() ) {
				// Uh-oh, table columns are already fixed.
				// TODO: alter to put this data in SummaryData table?
				throw genfile::BadArgumentError( "qcdb::FlatTableMultiVariantDBOutputter::add_variable()", "variable=\"" + variable + "\"" ) ;
			}
			else {
				// Still have time to add the variable to our list of variables, retaining the order of addition.
				m_variables.left.insert( VariableMap::left_value_type( variable, m_variables.size() ) ) ;
			}
		}
	}

	void FlatTableMultiVariantDBOutputter::create_new_key(
		Key const& key
	 ) {
		if( key.size() > m_key_entry_names.size() ) {
			throw genfile::BadArgumentError(
				"qcdb::FlatTableMultiVariantDBOutputter::create_new_key()",
				"key",
				(
					boost::format(
						"number of fields is larger than expected (%d, expected at most %d)"
					) % key.size() % m_key_entry_names.size()
				).str()
			) ;
		}
 		if( m_keys.size() == m_max_variants_per_block ) {
 			store_block() ;
 			m_keys.clear() ;
 			m_values.clear() ;
 		}
		m_keys.push_back( key ) ;
	}

	void FlatTableMultiVariantDBOutputter::store_data_for_key(
		Key const& key,
		std::string const& variable,
		genfile::VariantEntry const& value
	) {
		bool const new_key = m_keys.empty() || key != m_keys.back() ;
		if( new_key ) {
			create_new_key( key ) ;
		}

		VariableMap::left_const_iterator where = m_variables.left.find( variable ) ;
		if( where == m_variables.left.end() ) {
			if( m_insert_data_sql.get() ) {
				// Uh-oh, table columns are already fixed.
				// TODO: alter to put this data in SummaryData table?
				throw genfile::BadArgumentError( "qcdb::FlatTableMultiVariantDBOutputter::store_per_variant_pair_data()", "variable=\"" + variable + "\"" ) ;
			}
			else {
				// Still have time to add the variable to our list of variables, retaining the order of addition.
				where = m_variables.left.insert( VariableMap::left_value_type( variable, m_variables.size() ) ).first ;
			}
		}

		// Store the value of this variable
		m_values[ std::make_pair( m_keys.size() - 1, where->second ) ] = value ;
	}

	void FlatTableMultiVariantDBOutputter::store_block() {
		genfile::db::Connection::ScopedTransactionPtr transaction = m_outputter.connection().open_transaction( 7200 ) ;

		if( !m_insert_data_sql.get() ) {
			create_schema() ;
			// create_variables() ;
		}
		std::vector< genfile::db::Connection::RowId > variant_ids ;
		for( std::size_t key_i = 0; key_i < m_keys.size(); ++key_i ) {
			variant_ids.resize( m_keys[key_i].size() ) ;
			for( std::size_t i = 0; i < m_keys[key_i].size(); ++i ) {
				variant_ids[i] = m_outputter.get_or_create_variant( m_keys[key_i][i] ) ;
			}
			store_data_for_variants( key_i, m_outputter.analysis_id(), variant_ids ) ;
		}
	}

	std::string FlatTableMultiVariantDBOutputter::get_table_name() const {
		return m_table_name ;
	}

	void FlatTableMultiVariantDBOutputter::create_schema() {
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
				;
		for( std::size_t i = 0; i < m_key_entry_names.size(); ++i ) {
			schema_sql << (i>0 ? ", " : " " )
				<< "`" << m_key_entry_names[i] << "_id` INT NOT NULL REFERENCES Variant( id )" ;
		}
		
		insert_data_sql << "INSERT INTO "
			<< table_name ;
		insert_data_sql_columns << "( analysis_id" ;
		insert_data_sql_values << "VALUES( ?1" ;

		for( std::size_t i = 0; i < m_key_entry_names.size(); ++i ) {
			insert_data_sql_columns << ", `" << m_key_entry_names[i] << "_id`" ;
			insert_data_sql_values << boost::format( ", ?%d") % (i+2) ;
		}
		//insert_data_sql_columns << ")" ;
		//insert_data_sql_values << ")" ;

		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;

		for( std::size_t bind_i = m_key_entry_names.size()+2; var_i != end_var_i; ++var_i, ++bind_i ) {
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

		if( m_without_rowid ) {
			schema_sql << ", PRIMARY KEY( `analysis_id`" ;
			for( std::size_t i = 0; i < m_key_entry_names.size(); ++i ) {
				schema_sql  << ", `" << m_key_entry_names[i] << "_id`" ;
			}
			schema_sql << " ) ) WITHOUT ROWID ;" ;
		} else {
			schema_sql << " ) ; " ;
		}

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
			<< "CREATE VIEW IF NOT EXISTS `" << table_name << "View` AS SELECT A.name AS `analysis`, "
		;
		for( std::size_t i = 0; i < m_key_entry_names.size(); ++i ) {
			view_sql << boost::format( " V%d.rsid AS `%s_rsid`, V%d.chromosome AS `%s_chromosome`, V%d.position AS `%s_position`," )
				% (i+1)
				% m_key_entry_names[i]
				% (i+1)
				% m_key_entry_names[i]
				% (i+1)
				% m_key_entry_names[i]
			;
		}
		view_sql  << "T.* FROM `"
			<< table_name << "` T " ;
		for( std::size_t i = 0; i < m_key_entry_names.size(); ++i ) {
			view_sql << (boost::format( "INNER JOIN Variant V%d ON V%d.id = " + m_key_entry_names[i] + "_id " )
				% (i+1) % (i+1) ) ;
		}
		view_sql << "INNER JOIN Analysis A ON A.id = T.analysis_id" ;

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

	void FlatTableMultiVariantDBOutputter::store_data_for_variants(
		std::size_t const key_i,
		genfile::db::Connection::RowId const analysis_id,
		std::vector< genfile::db::Connection::RowId > const& variant_ids
	) {
		m_insert_data_sql->bind( 1, analysis_id ) ;
		std::size_t i = 0 ;
		for( ; i < variant_ids.size(); ++i ) {
			m_insert_data_sql->bind( i+2, variant_ids[i] ) ;
		}
		
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;
		
		for( std::size_t bind_i = i+2; var_i != end_var_i; ++var_i, ++bind_i ) {
			ValueMap::const_iterator where = m_values.find( std::make_pair( key_i, var_i->first )) ;
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
