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

	void write_snp_and_vector_to_sink(
		boost::shared_ptr< statfile::BuiltInTypeStatSink> sink,
		genfile::SNPIdentifyingData snp,
		Eigen::VectorXd const& vector,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_names
	) ;
}
#endif
