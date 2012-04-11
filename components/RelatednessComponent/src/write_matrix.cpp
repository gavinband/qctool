
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef RELATEDNESS_COMPONENT_WRITE_MATRIX_HPP
#define RELATEDNESS_COMPONENT_WRITE_MATRIX_HPP

#include <string>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <Eigen/Core>
#include "genfile/VariantEntry.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/string_utils.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "components/RelatednessComponent/write_matrix.hpp"
#include "components/RelatednessComponent/names.hpp"

namespace pca {
	std::string get_metadata( std::string const& source, std::string const& description ) {
		return "Created by qctool ("
			+ source
			+ "), " + appcontext::get_current_time_as_string() + "\n"
			+ description ;
	}

	void write_matrix_to_sink(
		statfile::BuiltInTypeStatSink::UniquePtr sink,
		Eigen::MatrixXd const& matrix,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names = 0,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names = 0
	) {
		sink->write_metadata( get_metadata( source, description ) ) ;

		using genfile::string_utils::to_string ;
		
		if( !get_column_names ) {
			get_column_names = boost::bind( &pca::string_and_number, "sample_", _1 ) ;
		}
		if( get_row_names ) {
			(*sink) | "id" ;
		}
		for( int j = 0; j < matrix.cols(); ++j ) {
			(*sink) | to_string( get_column_names( j ) ) ;
		}

		for( int i = 0; i < matrix.rows(); ++i ) {
			if( get_row_names ) {
				(*sink) << to_string( get_row_names(i) ) ;
			}
			for( int j = 0; j < matrix.cols(); ++j ) {
				double const& value = matrix( i, j ) ;
				if( value == value ) {
					(*sink) << value ;
				}
				else {
					(*sink) << "NA" ;
				}
			}
			(*sink) << statfile::end_row() ;
		}
	}
	
	void write_matrix(
		std::string const& filename,
		Eigen::MatrixXd const& matrix,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names = 0,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names = 0
	) {
		return write_matrix_to_sink(
			statfile::BuiltInTypeStatSink::open( filename ),
			matrix, source, description, get_row_names, get_column_names
		) ;
	}

	void write_snp_and_vector_to_sink(
		boost::shared_ptr< statfile::BuiltInTypeStatSink > sink,
		genfile::SNPIdentifyingData snp,
		Eigen::VectorXd const& vector,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_names
	) {
		if( sink->number_of_rows_written() == 0 && sink->current_column() == 0 ) {
			(*sink) | "SNPID" | "rsid" | "chromosome" | "position" | "allele_A" | "allele_B" ;
			if( get_names ) {
				for( int i = 0; i < vector.size(); ++i ) {
					(*sink) | get_names( i ).as< std::string >() ;
				}
			}
			else {
				for( int i = 0; i < vector.size(); ++i ) {
					(*sink) | ( "v" + genfile::string_utils::to_string( i )) ;
				}
			}
		}
		(*sink)
			<< snp.get_SNPID()
			<< snp.get_rsid()
			<< std::string( snp.get_position().chromosome() )
			<< snp.get_position().position()
			<< snp.get_first_allele()
			<< snp.get_second_allele() ;
		for( int i = 0; i < vector.size(); ++i ) {
			(*sink) << vector(i) ;
		}
		(*sink) << statfile::end_row() ;
	}
}
#endif
