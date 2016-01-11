
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
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/string_utils.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "statfile/DelimitedStatSink.hpp"
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

	void write_sample_file(
		std::string const& filename,
		genfile::CohortIndividualSource const& samples,
		Eigen::MatrixXd const& matrix,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names = 0,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names = 0
	) {
		assert( samples.get_number_of_individuals() == matrix.rows() ) ;
		statfile::BuiltInTypeStatSink::UniquePtr sink = statfile::BuiltInTypeStatSink::open( filename ) ;
		genfile::CohortIndividualSource::ColumnSpec const spec = samples.get_column_spec() ;
		if( !get_column_names ) {
			get_column_names = boost::bind( &pca::string_and_number, "sample_", _1 ) ;
		}
		{
			std::stringstream str ;
			str << "metadata: { " ;
			for( std::size_t j = 0; j < spec.size(); ++j ) {
				str << ( ( j > 0 ) ? ", " : "" ) ;
				str << "\"" << spec[j].name() << "\": "
					<< "{ \"type\": \"" << spec[j].type() << "\" }" ;
			}
			for( int j = 0; j < matrix.cols(); ++j ) {
				str << ", \"" << get_column_names( j ) << "\": "
					<< "{ \"type\": \"C\" }" ;
			}
			str << " }" ;
			sink->write_metadata( get_metadata( source, description + "\n" + str.str() ) ) ;
		}
		{
			for( std::size_t j = 0; j < spec.size(); ++j ) {
				(*sink) | spec[j].name() ;
			}
			for( int j = 0; j < matrix.cols(); ++j ) {
				(*sink) | get_column_names( j ).as< std::string >() ;
			}
		}
		for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
			for( std::size_t j = 0; j < spec.size(); ++j ) {
				(*sink) << samples.get_entry( i, spec[j].name() ) ;
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

	void write_matrix_lower_diagonals_in_long_form(
		std::string const& filename,
		Eigen::MatrixXd const& matrix1,
		Eigen::MatrixXd const& matrix2,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names
	) {
		assert(
			matrix1.cols() == matrix1.rows()
			&& matrix1.rows() == matrix2.rows()
			&& matrix1.cols() == matrix2.cols()
		) ;
		using genfile::string_utils::to_string ;
		statfile::BuiltInTypeStatSink::UniquePtr sink = statfile::BuiltInTypeStatSink::open( filename ) ;
		sink->write_metadata( get_metadata( source, description ) ) ;
		(*sink) | "sample_1" | "sample_2" | "pairwise.complete.obs" | "value" ;
		for( int i = 0; i < matrix1.rows(); ++i ) {
			genfile::VariantEntry row_name = get_row_names ? get_row_names( i ) : ( "row_" + to_string( i ) ) ;
			for( int j = 0; j <= i; ++j ) {
				(*sink)
					<< row_name
					<< ( get_column_names ? get_column_names( j ) : ( "column_" + to_string( j ) ) )
					<< matrix1( i, j )
					<< matrix2( i, j )
				;
				(*sink) << statfile::end_row() ;
			}
		}
	}

	void write_loadings_to_sink(
		boost::shared_ptr< statfile::BuiltInTypeStatSink > sink,
		genfile::VariantIdentifyingData snp,
		double const non_missingness,
		double const allele_frequency,
		Eigen::VectorXd const& vector,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_names
	) {
		if( sink->number_of_rows_written() == 0 && sink->current_column() == 0 ) {
			(*sink) | "SNPID" | "rsid" | "chromosome" | "position" | "allele_A" | "allele_B" | "N" | "B_allele_frequency" ;
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
			<< snp.get_identifiers_as_string( ",", 1, snp.number_of_identifiers() )
			<< snp.get_primary_id()
			<< std::string( snp.get_position().chromosome() )
			<< snp.get_position().position()
			<< snp.get_allele(0)
			<< snp.get_allele(1)
			<< non_missingness
			<< allele_frequency
		;
		for( int i = 0; i < vector.size(); ++i ) {
			(*sink) << vector(i) ;
		}
		(*sink) << statfile::end_row() ;
	}
}
#endif
