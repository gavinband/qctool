
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef SQLITEHAPLOTYPESSNPDATASINK_HPP
#define SQLITEHAPLOTYPESSNPDATASINK_HPP

#include <iostream>
#include <string>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/GenLikeSNPDataSink.hpp"
#include "qcdb/DBOutputter.hpp"

	// This class represents a SNPDataSink which writes its data
	// to a plain GEN file.
	class SQLiteHaplotypesSNPDataSink: public genfile::SNPDataSink
	{
	public:
		SQLiteHaplotypesSNPDataSink( qcdb::DBOutputter::UniquePtr outputter ) ;
		
		operator bool() const { return true ; }
		std::string get_spec() const ;
		
	private:
		qcdb::DBOutputter::UniquePtr m_outputter ;
		std::size_t m_number_of_samples ;
		db::Connection::StatementPtr m_insert_sample_stmnt ;
		db::Connection::StatementPtr m_insert_data_stmnt ;

		std::vector< genfile::SNPIdentifyingData2 > m_snps ;
		std::vector< std::vector< char > > m_data ;
		std::size_t m_data_i ;

	private:
		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) ;

		void write_variant_data_impl(
			genfile::SNPIdentifyingData const& id_data,
			genfile::VariantDataReader& data_reader,
			Info const& info = Info()
		) ;
		
		void finalise_impl() ;
		
		void flush_data( std::size_t const data_count ) ;
	} ;

#endif