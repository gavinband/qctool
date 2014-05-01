
//			Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//	  (See accompanying file LICENSE_1_0.txt or copy at
//			http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread/thread.hpp>
#include "genfile/SNPIdentifyingData2.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "qcdb/FlatTableDBOutputter.hpp"

// #define DEBUG_FLATTABLEDBOUTPUTTER 1

namespace qcdb {
	FlatTableDBOutputter::UniquePtr FlatTableDBOutputter::create( std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata ) {
		return UniquePtr( new FlatTableDBOutputter( filename, analysis_name, analysis_description, metadata ) ) ;
	}

	FlatTableDBOutputter::SharedPtr FlatTableDBOutputter::create_shared( std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata ) {
		return SharedPtr( new FlatTableDBOutputter( filename, analysis_name, analysis_description, metadata ) ) ;
	}

	FlatTableDBOutputter::FlatTableDBOutputter( std::string const& filename, std::string const& analysis_name, std::string const& analysis_description, Metadata const& metadata ):
		m_outputter( filename, analysis_name, analysis_description, metadata ),
		m_table_name( "Analysis" + genfile::string_utils::to_string( m_outputter.analysis_id() ) ),
		m_max_snps_per_block( 1000 )
	{}

	FlatTableDBOutputter::~FlatTableDBOutputter() {
		store_block() ;
		m_snps.clear() ;
		m_values.clear() ;
		m_outputter.finalise() ;
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
		m_snps.clear() ;
		m_values.clear() ;

		if( options & eCreateIndices ) {
			db::Connection::ScopedTransactionPtr transaction = m_outputter.connection().open_transaction( 600 ) ;
			m_outputter.connection().run_statement( "CREATE INDEX IF NOT EXISTS " + m_table_name + "_index ON " + m_table_name + "( variant_id )" ) ;
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
				m_variables.left.insert( VariableMap::left_value_type( variable, m_variables.size() ) ).first ;
			}
		}
	}

	void FlatTableDBOutputter::create_new_variant( genfile::SNPIdentifyingData2 const& snp ) {
		if( m_snps.size() == m_max_snps_per_block ) {
			store_block() ;
			m_snps.clear() ;
			m_values.clear() ;
		}
		m_snps.push_back( snp ) ;
	}

	void FlatTableDBOutputter::store_per_variant_data(
		genfile::SNPIdentifyingData2 const& snp,
		std::string const& variable,
		genfile::VariantEntry const& value
	) {
		bool const new_snp = m_snps.empty() || snp != m_snps.back() ;
		if( new_snp ) {
			// If we have a whole block's worth of data, store it now.
			if( m_snps.size() == m_max_snps_per_block ) {
				store_block() ;
				m_snps.clear() ;
				m_values.clear() ;
			}
			m_snps.push_back( snp ) ;
		}

		VariableMap::left_const_iterator where = m_variables.left.find( variable ) ;
		if( where == m_variables.left.end() ) {
			if( m_insert_data_sql.get() ) {
				// Uh-oh, table columns are already fixed.
				// TODO: alter to put this data in SummaryData table?
				throw genfile::BadArgumentError( "qcdb::FlatTableDBOutputter::store_per_variant_data()", "variable=\"" + variable + "\"" ) ;
			}
			else {
				// Still have time to add the variable to our list of variables, retaining the order of addition.
				where = m_variables.left.insert( VariableMap::left_value_type( variable, m_variables.size() ) ).first ;
			}
		}

		// Store the value of this variable
		m_values[ std::make_pair( m_snps.size() - 1, where->second ) ] = value ;
	}

	void FlatTableDBOutputter::store_block() {
		db::Connection::ScopedTransactionPtr transaction = m_outputter.connection().open_transaction( 3600 ) ;

		if( !m_insert_data_sql.get() ) {
			create_schema() ;
			create_variables() ;
		}
		for( std::size_t i = 0; i < m_snps.size(); ++i ) {
			db::Connection::RowId const variant_id = m_outputter.get_or_create_variant( m_snps[i] ) ;
			store_data_for_variant( i, m_snps[i], m_outputter.analysis_id(), variant_id ) ;
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
			"variant_id INT NOT NULL REFERENCES Variant( id ), "
			"chromosome TEXT, "
			"position INTEGER"
		;
		
		insert_data_sql << "INSERT INTO "
			<< table_name ;
		insert_data_sql_columns << "( analysis_id, variant_id, chromosome, position" ;
		insert_data_sql_values << "VALUES( ?1, ?2, ?3, ?4" ;
		
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;

		for( std::size_t bind_i = 5; var_i != end_var_i; ++var_i, ++bind_i ) {
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
			<< "CREATE VIEW IF NOT EXISTS \"" << table_name << "View\" AS "
			<< "SELECT V.rsid, V.alleleA, V.alleleB, A.name AS analysis, T.* FROM \""
			<< table_name << "\" T "
			<< "INNER JOIN Variant V ON V.id = T.variant_id "
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

#if DEBUG_FLATTABLEDBOUTPUTTER
		std::cerr << "Getting \"table\" entry...\n"  ;
		db::Connection::RowId table_id = m_outputter.get_or_create_entity( "table", "Table holding results of an analysis" ) ;
		std::cerr << "ok.\n" ;
#endif

		m_outputter.get_or_create_entity_data(
			m_outputter.analysis_id(),
			m_outputter.get_or_create_entity( "table", "Table holding results of an analysis" ),
			table_name + "View"
		) ;
	}

	void FlatTableDBOutputter::create_variables() {
		db::Connection::RowId const variable_class_id = m_outputter.get_or_create_entity( "per-variant variable", "per-variant variable values" ) ;
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;
		for( ; var_i != end_var_i; ++var_i ) {
			m_outputter.get_or_create_entity( var_i->second, var_i->second, variable_class_id ) ;
		}
	}
	
	void FlatTableDBOutputter::store_data_for_variant(
		std::size_t const snp_i,
		genfile::SNPIdentifyingData2 const& snp,
		db::Connection::RowId const analysis_id,
		db::Connection::RowId const variant_id
	) {
		m_insert_data_sql->bind( 1, analysis_id ) ;
		m_insert_data_sql->bind( 2, variant_id ) ;
		m_insert_data_sql->bind( 3, std::string( snp.get_position().chromosome() ) ) ;
		m_insert_data_sql->bind( 4, snp.get_position().position() ) ;
		
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;
		
		for( std::size_t bind_i = 5; var_i != end_var_i; ++var_i, ++bind_i ) {
			ValueMap::const_iterator where = m_values.find( std::make_pair( snp_i, var_i->first )) ;
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
