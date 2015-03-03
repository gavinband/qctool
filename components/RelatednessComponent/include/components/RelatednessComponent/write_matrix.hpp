
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef RELATEDNESS_COMPONENT_WRITE_MATRIX_HPP
#define RELATEDNESS_COMPONENT_WRITE_MATRIX_HPP

#include <string>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <Eigen/Core>
#include "genfile/VariantEntry.hpp"
#include "genfile/string_utils.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"

namespace pca {
	std::string get_metadata( std::string const& source, std::string const& description ) ;
	
	void write_matrix_to_sink(
		statfile::BuiltInTypeStatSink::UniquePtr sink,
		Eigen::MatrixXd const& matrix,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names = 0,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names = 0
	) ;

	void write_matrix(
		std::string const& filename,
		Eigen::MatrixXd const& matrix,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names = 0,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names = 0
	) ;

	void write_sample_file(
		std::string const& filename,
		genfile::CohortIndividualSource const& samples,
		Eigen::MatrixXd const& matrix,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names = 0,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names = 0
	) ;

	// Write lower diagonal of matrix with 4 columns:
	// 
	void write_matrix_lower_diagonals_in_long_form(
		std::string const& filename,
		Eigen::MatrixXd const& matrix1,
		Eigen::MatrixXd const& matrix2,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names = 0,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names = 0
	) ;

	void write_loadings_to_sink(
		boost::shared_ptr< statfile::BuiltInTypeStatSink> sink,
		genfile::SNPIdentifyingData snp,
		double const non_missingness,
		double const allele_frequency,
		Eigen::VectorXd const& vector,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_names
	) ;
}
#endif
