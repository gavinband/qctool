#include <iostream>
#include <iomanip>
#include <string>
#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/Error.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "components/RelatednessComponent/UDUTDecompositionLoader.hpp"

namespace relatedness {
	UDUTDecompositionLoader::UDUTDecompositionLoader( genfile::CohortIndividualSource const& samples ):
	 	m_samples( samples )
	{}

	void UDUTDecompositionLoader::send_UDUT_to( ResultCallback callback ) {
		m_result_signal.connect( callback ) ;
	}

	void UDUTDecompositionLoader::load_matrix( std::string const& filename ) const {
		Eigen::MatrixXd matrix ;
		std::size_t number_of_snps ;
		load_matrix_impl( filename, &matrix, &number_of_snps ) ;
		m_result_signal( matrix, m_samples.get_number_of_individuals(), number_of_snps ) ;
	}

	void UDUTDecompositionLoader::load_matrix_impl( std::string const& filename, Eigen::MatrixXd* matrix, std::size_t* number_of_snps ) const {
		statfile::BuiltInTypeStatSource::UniquePtr source = statfile::BuiltInTypeStatSource::open( filename ) ;
		assert( matrix ) ;
		assert( number_of_snps ) ;
		
		using namespace genfile::string_utils ;
		// Read the metadata
		std::size_t number_of_samples = 0 ;
		std::vector< std::string > metadata = split( source->get_descriptive_text(), "\n" ) ;
		for( std::size_t i = 0; i < metadata.size(); ++i ) {
			std::vector< std::string > elts = split_and_strip( metadata[i], ":", " " ) ;
			if( elts.size() == 2 && elts[0] == "Number of SNPs" ) {
				*number_of_snps = to_repr< std::size_t >( elts[1] ) ;
			}
			else if( elts.size() == 2 && elts[0] == "Number of samples" ) {
				number_of_samples = to_repr< std::size_t >( elts[1] ) ;
			}
		}
		if( number_of_samples != m_samples.get_number_of_individuals() ) {
			throw genfile::MismatchError(
				"UDUTDecompositionLoader::load_matrix()",
				filename,
				"number of samples: " + to_string( number_of_samples ),
				"expected number: " + to_string( m_samples.get_number_of_individuals() )
			) ;
		}

		if( source->number_of_columns() != number_of_samples + 1 ) {
			throw genfile::MalformedInputError( source->get_source_spec(), 0, std::min( source->number_of_columns(), number_of_samples + 1 )) ;
		}
		if( source->number_of_rows() != number_of_samples ) {
			throw genfile::MalformedInputError( source->get_source_spec(), 0, std::min( source->number_of_rows(), number_of_samples )) ;
		}
		// Read the matrix, making sure the samples come in the same order as in the sample file.
		matrix->resize( number_of_samples, number_of_samples+1 ) ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			for( std::size_t j = 0; j < ( number_of_samples + 1 ); ++j ) {
				(*source) >> (*matrix)(i,j) ;
			}
			(*source) >> statfile::end_row() ;
		}
	}
}
