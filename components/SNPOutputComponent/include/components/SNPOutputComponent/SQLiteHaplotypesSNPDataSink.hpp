
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
		SQLiteHaplotypesSNPDataSink(
			qcdb::DBOutputter::UniquePtr outputter,
			genfile::CohortIndividualSource const& samples
		) ;
		
		operator bool() const { return true ; }
		std::string get_spec() const ;
		
        void set_genotype_field( std::string field ) ;
		
	private:
		std::string m_genotype_field ; 
		qcdb::DBOutputter::UniquePtr m_outputter ;
		genfile::CohortIndividualSource const& m_samples ;
		db::Connection::StatementPtr m_insert_data_stmnt ;

		std::vector< genfile::VariantIdentifyingData > m_snps ;
		std::vector< std::vector< uint8_t > > m_data ;
		std::size_t m_data_i ;

	private:
		void set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) ;

		void write_variant_data_impl(
			genfile::VariantIdentifyingData const& id_data,
			genfile::VariantDataReader& data_reader,
			Info const& info = Info()
		) ;
		
		void finalise_impl() ;
		
		void flush_data( std::size_t const data_count ) ;
	} ;

#endif
