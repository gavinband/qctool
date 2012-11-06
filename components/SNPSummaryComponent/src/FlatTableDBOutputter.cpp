
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

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
#include "components/SNPSummaryComponent/FlatTableDBOutputter.hpp"
#include "../../qctool_version_autogenerated.hpp"

namespace snp_summary_component {
	FlatTableDBOutputter::UniquePtr FlatTableDBOutputter::create( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) {
		return UniquePtr( new FlatTableDBOutputter( filename, cohort_name, metadata ) ) ;
	}
	FlatTableDBOutputter::SharedPtr FlatTableDBOutputter::create_shared( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ) {
		return SharedPtr( new FlatTableDBOutputter( filename, cohort_name, metadata ) ) ;
	}

	FlatTableDBOutputter::FlatTableDBOutputter( std::string const& filename, std::string const& cohort_name, Metadata const& metadata ):
		m_outputter( filename, cohort_name, metadata ),
		m_max_snps_per_block( 5 )
	{}

	FlatTableDBOutputter::~FlatTableDBOutputter() {
		store_block() ;
		m_outputter.finalise() ;
	}
	
	void FlatTableDBOutputter::finalise() {
		store_block() ;
		m_outputter.finalise() ;
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
				throw genfile::BadArgumentError( "snp_summary_component::FlatTableDBOutputter::store_per_variant_data()", "variable=\"" + variable + "\"" ) ;
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
		db::Connection::ScopedTransactionPtr transaction = m_outputter.connection().open_transaction( 240 ) ; // wait 4 minutes if we have to.

		if( !transaction.get() ) {
			throw genfile::OperationFailedError( "SNPSummaryComponent::FlatTableDBOutputter::write_data()", m_outputter.connection().get_spec(), "Opening transaction." ) ;
		}
		
		if( !m_insert_data_sql.get() ) {
			create_schema() ;
		}
		for( std::size_t i = 0; i < m_snps.size(); ++i ) {
			db::Connection::RowId const variant_id = m_outputter.get_or_create_variant( m_snps[i] ) ;
			store_data_for_variant( i, m_outputter.analysis_id(), variant_id ) ;
		}
	}

	void FlatTableDBOutputter::create_schema() {
		using genfile::string_utils::to_string ;
		std::ostringstream schema_sql ;
		std::ostringstream index_sql ;
		std::ostringstream insert_data_sql ;
		std::string const table_name = "Analysis" + to_string( m_outputter.analysis_id() )  ;
		schema_sql << "CREATE TABLE "
			<< table_name
			<< " ( "
			"analysis_id INT NOT NULL REFERENCES Entity( id ), "
			"variant_id INT NOT NULL REFERENCES Variant( id )"
		;
		
		insert_data_sql << "INSERT INTO "
			<< table_name
			<< " VALUES( ?1, ?2" ;
		
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;

		for( std::size_t bind_i = 3; var_i != end_var_i; ++var_i, ++bind_i ) {
			schema_sql << ", " ;
			insert_data_sql << ", " ;
			schema_sql
				<< '"'
				<< var_i->second
				<< '"'
				<< " NULL" ;
				
			insert_data_sql << "?" << to_string( bind_i ) ;
		}

		schema_sql << " ) ; " ;

		insert_data_sql << ") ; " ;

		std::cerr << "Creating table " << table_name << " using this SQL:\n"
			<< schema_sql.str()
			<< "\n" ;
		
		m_outputter.connection().run_statement( schema_sql.str() ) ;
		m_outputter.connection().run_statement( "CREATE INDEX " + table_name + "_index ON " + table_name + "( variant_id )" ) ;

		std::ostringstream view_sql ;	
		view_sql
			<< "CREATE VIEW \"" << table_name << "View\" AS "
			<< "SELECT V.chromosome, V.position, V.rsid, V.alleleA, V.alleleB, A.name AS analysis_id, T.* FROM \""
			<< table_name << "\" T "
			<< "INNER JOIN Variant V ON V.id = T.variant_id "
			<< "INNER JOIN Entity A ON A.id = T.analysis_id"
		;

		std::cerr << "Creating view using SQL:\n" << view_sql.str() << "\n" ;
		m_outputter.connection().run_statement( view_sql.str() ) ;

		std::cerr << "Inserts will use this SQL:\n"
			<< insert_data_sql.str()
			<< "\n" ;
		
		m_insert_data_sql = m_outputter.connection().get_statement(
			insert_data_sql.str()
		) ;
	}

	void FlatTableDBOutputter::store_data_for_variant(
		std::size_t const snp_i,
		db::Connection::RowId const analysis_id,
		db::Connection::RowId const variant_id
	) {
		m_insert_data_sql->bind( 1, analysis_id ) ;
		m_insert_data_sql->bind( 2, variant_id ) ;
		
		VariableMap::right_const_iterator
			var_i = m_variables.right.begin(),
			end_var_i = m_variables.right.end() ;
		
		for( std::size_t bind_i = 3; var_i != end_var_i; ++var_i, ++bind_i ) {
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
