
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread/thread.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/db/Connection.hpp"
#include "genfile/db/SQLStatement.hpp"
#include "components/SampleSummaryComponent/FlatTableDBOutputter.hpp"
#include "qcdb/DBOutputter.hpp"

//#define DEBUG_FLATTABLEDBOUTPUTTER 1

namespace sample_stats {
	FlatTableDBOutputter::UniquePtr FlatTableDBOutputter::create(
		std::string const& filename,
		std::string const& analysis_name,
		std::string const& analysis_description,
		qcdb::DBOutputter::Metadata const& metadata,
		genfile::CohortIndividualSource const& samples,
		boost::optional< genfile::db::Connection::RowId > analysis_id
	) {
		return UniquePtr( new FlatTableDBOutputter( filename, analysis_name, analysis_description, metadata, samples, analysis_id ) ) ;
	}

	FlatTableDBOutputter::SharedPtr FlatTableDBOutputter::create_shared(
		std::string const& filename,
		std::string const& analysis_name,
		std::string const& analysis_description,
		qcdb::DBOutputter::Metadata const& metadata,
		genfile::CohortIndividualSource const& samples,
		boost::optional< genfile::db::Connection::RowId > analysis_id
	) {
		return SharedPtr( new FlatTableDBOutputter( filename, analysis_name, analysis_description, metadata, samples, analysis_id ) ) ;
	}

	FlatTableDBOutputter::FlatTableDBOutputter(
		std::string const& filename,
		std::string const& analysis_name,
		std::string const& analysis_description,
		qcdb::DBOutputter::Metadata const& metadata,
		genfile::CohortIndividualSource const& samples, // unused at present
		boost::optional< genfile::db::Connection::RowId > analysis_id
	):
		m_outputter( filename, analysis_name, analysis_description, metadata, analysis_id ),
		m_table_name( "Analysis" + genfile::string_utils::to_string( m_outputter.analysis_id() ) + "SampleData" )
	{
		genfile::db::Connection::ScopedTransactionPtr transaction = m_outputter.connection().open_transaction( 240 ) ;
		m_outputter.connection().run_statement(
				"CREATE TABLE IF NOT EXISTS Sample ( "
				"analysis_id INTEGER NOT NULL REFERENCES Entity( id ), "
				"id TEXT NOT NULL, "
				"index_in_data INTEGER NOT NULL, "
				"PRIMARY KEY( analysis_id, id )"
			")"
		) ;

		construct_statements() ;
		
		store_samples( samples ) ;
	}
	
	void FlatTableDBOutputter::set_table_name( std::string const& table_name ) {
		if( table_name.find_first_not_of( "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_" ) != std::string::npos ) {
			throw genfile::BadArgumentError( "qcdb::FlatTableDBOutputter::set_table_name()", "table_name=\"" + table_name + "\"" ) ;
		}
		m_table_name = table_name ;
	}
	
	
	SampleStorage::AnalysisId FlatTableDBOutputter::analysis_id() const {
		return 0 ;
	}
	
	void FlatTableDBOutputter::store_per_sample_data(
		std::string const& computation_name,
		std::size_t sample,
		std::string const& variable,
		std::string const& description,
		genfile::VariantEntry const& value
	) {
		VariableMap::left_const_iterator where = add_variable_impl( variable ) ;
		m_values[ std::make_pair( sample, where->second ) ] = value ;
	}

	void FlatTableDBOutputter::add_variable( std::string const& variable ) {
		// add a variable but do not return the iterator.
		add_variable_impl( variable ) ;
	}

	FlatTableDBOutputter::VariableMap::left_const_iterator FlatTableDBOutputter::add_variable_impl( std::string const& variable ) {
		VariableMap::left_const_iterator where = m_variables.left.find( variable ) ;
		if( where == m_variables.left.end() ) {
			if( m_insert_data_sql.get() ) {
				// Uh-oh, table columns are already fixed.
				// TODO: alter to put this data in SummaryData table?
				throw genfile::BadArgumentError( "sample_stats::FlatTableDBOutputter::add_variable()", "variable=\"" + variable + "\"" ) ;
			}
			else {
				// Still have time to add the variable to our list of variables, retaining the order of addition.
				where = m_variables.left.insert( VariableMap::left_value_type( variable, m_variables.size() ) ).first ;
			}
		}
		return where ;
	}
	
	void FlatTableDBOutputter::store_block() {
		genfile::db::Connection::ScopedTransactionPtr transaction = m_outputter.connection().open_transaction( 3600 ) ;

		if( !m_insert_data_sql.get() ) {
			create_schema() ;
			create_variables() ;
		}
		for( std::size_t i = 0; i < m_sample_ids.size(); ++i ) {
			store_data_for_sample( m_outputter.analysis_id(), i, m_sample_ids[i] ) ;
		}
	}
	
	void FlatTableDBOutputter::create_schema() {
		using genfile::string_utils::to_string ;
		std::ostringstream schema_sql ;
		std::ostringstream index_sql ;
		std::ostringstream view_sql ;
		std::ostringstream insert_data_sql ;
		std::ostringstream insert_data_sql_columns ;
		std::ostringstream insert_data_sql_values ;
		std::string const& table_name = m_table_name ;
		schema_sql << "CREATE TABLE IF NOT EXISTS "
			<< table_name
			<< " ( "
			"analysis_id INT NOT NULL REFERENCES Entity( id ), "
			"sample_id TEXT NOT NULL REFERENCES Sample( id )"
		;

		view_sql << "CREATE VIEW IF NOT EXISTS \""
			<< ( table_name + "View" )
			<< "\" AS "
			<< "SELECT T.analysis_id AS analysis_id, E.name AS analysis, T.sample_id AS sample_id, S.index_in_data AS index_in_data" ;

		insert_data_sql << "INSERT INTO "
			<< table_name ;
		insert_data_sql_columns << "( analysis_id, sample_id" ;
		insert_data_sql_values << "VALUES( ?1, ?2" ;

		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;

		for( std::size_t bind_i = 3; var_i != end_var_i; ++var_i, ++bind_i ) {
			schema_sql << ", " ;
			insert_data_sql_columns << ", " ;
			insert_data_sql_values << ", " ;
			schema_sql
				<< '"'
				<< var_i->second
				<< '"'
				<< " NULL" ;

			view_sql << ", " << '"' << var_i->second << '"' ;
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

		view_sql
			<< " FROM \""
			<< table_name << "\" T "
			<< "INNER JOIN Sample S ON S.analysis_id == T.analysis_id AND S.id = T.sample_id "
			<< "INNER JOIN Analysis E ON E.id == S.analysis_id "
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

	void FlatTableDBOutputter::create_variables() {
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;
		for( ; var_i != end_var_i; ++var_i ) {
			m_outputter.create_variable( m_table_name, var_i->second ) ;
		}
	}
	
	void FlatTableDBOutputter::construct_statements() {
		m_find_sample_statement = m_outputter.connection().get_statement( "SELECT id FROM Sample WHERE analysis_id == ?1 AND id == ?2" ) ;
		m_insert_sample_statement = m_outputter.connection().get_statement( "INSERT INTO Sample (analysis_id, id, index_in_data ) VALUES( ?1, ?2, ?3 )" ) ;
	}

	FlatTableDBOutputter::~FlatTableDBOutputter() {
	}

	void FlatTableDBOutputter::finalise( long options ) {
		if( !m_values.empty() ) {
			store_block() ;
			m_values.clear() ;
		}
		if( options & qcdb::eCreateIndices ) {
			genfile::db::Connection::ScopedTransactionPtr transaction = m_outputter.connection().open_transaction( 600 ) ;
			m_outputter.connection().run_statement( "CREATE INDEX IF NOT EXISTS SampleIndex ON Sample( id )" ) ;
			m_outputter.connection().run_statement( "CREATE INDEX IF NOT EXISTS " + m_table_name + "_index ON " + m_table_name + "( sample_id )" ) ;
			m_outputter.connection().run_statement( "CREATE INDEX IF NOT EXISTS " + m_table_name + "_index2 ON " + m_table_name + "( analysis_id, sample_id )" ) ;
		}

		m_outputter.finalise( options ) ;
	}

	void FlatTableDBOutputter::store_samples( genfile::CohortIndividualSource const& samples ) {
		m_sample_ids.resize( samples.get_number_of_individuals() ) ;
		for( std::size_t index = 0; index < samples.get_number_of_individuals(); ++index ) {
			genfile::VariantEntry const id = samples.get_entry( index, "ID_1" ) ;
			get_or_create_sample( id, index ) ;
			m_sample_ids[index] = id ;
		}
	}

	genfile::db::Connection::RowId FlatTableDBOutputter::get_or_create_sample( genfile::VariantEntry const& identifier, std::size_t index ) const {
		genfile::db::Connection::RowId result ;

		if( identifier.is_missing() ) {
			throw genfile::BadArgumentError( "impl::FlatTableDBOutputter::get_or_create_sample()", "identifier=NA" ) ;
		}

		m_find_sample_statement
			->bind( 1, m_outputter.analysis_id() )
			.bind( 2, identifier )
			.step() ;

		if( m_find_sample_statement->empty() ) {
			m_insert_sample_statement
				->bind( 1, m_outputter.analysis_id() )
				.bind( 2, identifier )
				.bind( 3, int64_t( index ) )
				.step() ;
				
			result = m_outputter.connection().get_last_insert_row_id() ;

			m_insert_sample_statement->reset() ;
		}
		else {
			result = m_find_sample_statement->get< genfile::db::Connection::RowId >( 0 ) ;
		}

		m_find_sample_statement->reset() ;

		return result ;
	}
	
	void FlatTableDBOutputter::store_data_for_sample( genfile::db::Connection::RowId analysis_id, std::size_t sample_index, genfile::VariantEntry const& sample_id ) {
		m_insert_data_sql->bind( 1, analysis_id ) ;
		m_insert_data_sql->bind( 2, sample_id ) ;
		
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;
		
		for( std::size_t bind_i = 3; var_i != end_var_i; ++var_i, ++bind_i ) {
			ValueMap::const_iterator where = m_values.find( std::make_pair( sample_index, var_i->first )) ;
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
