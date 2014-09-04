
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SNP_OUTPUT_COMPONENT_STORE_SAMPLES_HPP
#define SNP_OUTPUT_COMPONENT_STORE_SAMPLES_HPP

#include <iostream>
#include <string>
#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/string_utils.hpp"
#include "qcdb/DBOutputter.hpp"
#include "qcdb/store_samples_in_db.hpp"

namespace qcdb {
	void store_samples_in_db(
		std::size_t number_of_samples,
		genfile::SNPDataSink::SampleNameGetter getter,
		genfile::CohortIndividualSource const& sample_data,
		qcdb::DBOutputter& outputter,
		std::string const& table_name
	) {
		if( number_of_samples != sample_data.get_number_of_individuals() ) {
			throw genfile::BadArgumentError(
				"snp_output_component::store_samples_in_db()",
				"number_of_samples",
				"number of samples ("
				+ genfile::string_utils::to_string( number_of_samples )
				+ ") does not match that in sample files ("
				+ genfile::string_utils::to_string( sample_data.get_number_of_individuals() )
				+ ")."
			) ;
		}
		genfile::CohortIndividualSource::ColumnSpec const spec = sample_data.get_column_spec() ;

		db::Connection::ScopedTransactionPtr transaction = outputter.connection().open_transaction( 2400 ) ;

		std::string sample_sql = 
			"CREATE TABLE IF NOT EXISTS Sample ( "
			"analysis_id INTEGER NOT NULL REFERENCES Analysis( id ), "
			"index_in_data INTEGER NOT NULL" ;

		std::string insert_sample_sql = 
			"INSERT INTO Sample ( analysis_id, index_in_data" ;

		for( std::size_t i = 0; i < spec.size(); ++i ) {
			sample_sql += ", " + spec[i].name() + " " + ( spec[i].is_continuous() ? "FLOAT" : "TEXT" ) ;
			insert_sample_sql += ", " + spec[i].name() ;
		}
		sample_sql += ", UNIQUE( analysis_id, \"" + spec[0].name() + "\" ), "
			+ "UNIQUE( analysis_id, index_in_data )"
			+ ")" ;

		insert_sample_sql += " ) VALUES ( ?, ?" ;
		for( std::size_t i = 0; i < spec.size(); ++i ) {
			insert_sample_sql += ", ?" ;
		}
		insert_sample_sql += " ) ;" ;

		outputter.connection().run_statement( sample_sql ) ;
	
		db::Connection::StatementPtr insert_sample_stmnt = outputter.connection().get_statement( insert_sample_sql ) ;

		for( std::size_t sample_i = 0; sample_i < number_of_samples; ++sample_i ) {
			if( getter(sample_i) != sample_data.get_entry( sample_i, "ID_1" )) {
				throw genfile::BadArgumentError(
					"snp_output_component::store_samples_in_db()",
					"getter",
					"getter(" + genfile::string_utils::to_string( sample_i )  + " = \""
						+ getter(sample_i).as< std::string >() + "\", but expected \""
					+ sample_data.get_entry( sample_i, "ID_1" ).as< std::string >()
					+ "\"."
				) ;
			}
			insert_sample_stmnt
				->bind( 1, outputter.analysis_id() )
				.bind( 2, int64_t( sample_i ) ) ;
			for( std::size_t column_i = 0; column_i < spec.size(); ++column_i ) {
				insert_sample_stmnt->bind( column_i+3, sample_data.get_entry( sample_i, spec[column_i].name() )) ;
			}
			insert_sample_stmnt->step() ;
			insert_sample_stmnt->reset() ;
		}		
	}
}

#endif
