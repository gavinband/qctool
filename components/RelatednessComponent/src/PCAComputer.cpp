#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/function.hpp>
#include "../config.hpp"
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "worker/Worker.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "components/RelatednessComponent/PCAComputer.hpp"
#include "components/RelatednessComponent/LapackEigenDecomposition.hpp"
#include "components/RelatednessComponent/names.hpp"

PCAComputer::PCAComputer(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	worker::Worker* worker,
	appcontext::UIContext& ui_context
):
	m_options( options ),
	m_ui_context( ui_context ),
	m_samples( samples ),
	m_number_of_samples( samples.get_number_of_individuals() ),
	m_number_of_PCAs_to_compute( 0 ),
	m_threshhold( 0.9 )
{
	assert( m_options.check( "-load-kinship" )) ;
	m_filename = m_options.get< std::string >( "-load-kinship" ) ;
	// remove .csv extension if present.
	std::string const extension = ".csv" ;
	if(
		m_filename.size() >= extension.size()
		&& ( m_filename.compare(
			m_filename.size() - extension.size(),
			extension.size(),
			extension
		) == 0 )
	) {
		m_filename.resize( m_filename.size() - extension.size() ) ;
	}
	load_matrix( m_filename + ".csv", &m_kinship_matrix, &m_number_of_snps ) ;

	if( m_options.check( "-PCAs" ) ) {	
		m_number_of_PCAs_to_compute = std::min( m_options.get< std::size_t >( "-PCAs" ), samples.get_number_of_individuals() ) ;
	}
}

namespace {
	genfile::VariantEntry get_eigendecomposition_header( std::string const& prefix, int value ) {
		if( value == 0 ) {
			return genfile::VariantEntry( "eigenvalue" ) ;
		}
		else {
			return genfile::VariantEntry( prefix + genfile::string_utils::to_string( value ) ) ;
		}
	}
	
	
	template< typename Stream, typename Matrix >
	void write_matrix_to_stream( Stream& stream, Matrix const& matrix ) {
		for( int i = 0; i < matrix.rows(); ++i ) {
			for( int j = 0; j < matrix.rows(); ++j ) {
				if( matrix( i,j ) != matrix( i,j )) {
					stream << std::setw( 10 ) << "NA" ;
				} else if( std::abs( matrix( i, j ) ) < 0.00001 ) {
					stream << std::setw( 10 ) << "0" ;
				}
				else {
					stream << std::setw( 10 ) << matrix(i,j) ;
				}
			}
			stream << "\n" ;
		}
		stream << "\n" ;
	}
}

void PCAComputer::load_matrix( std::string const& filename, Eigen::MatrixXd* matrix, std::size_t* number_of_snps ) const {
	statfile::BuiltInTypeStatSource::UniquePtr source = statfile::BuiltInTypeStatSource::open( filename ) ;

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
			"PCAComputer::load_matrix()",
			filename,
			"number of samples: " + to_string( number_of_samples ),
			"expected number: " + to_string( m_samples.get_number_of_individuals() )
		) ;
	}

	// Read the matrix, making sure the samples come in the same order as in the sample file.
	matrix->resize( m_samples.get_number_of_individuals(), m_samples.get_number_of_individuals() ) ;
	std::vector< std::size_t > sample_column_indices( m_samples.get_number_of_individuals() ) ;
	for( std::size_t sample_i = 0; sample_i < m_samples.get_number_of_individuals(); ++sample_i ) {
		sample_column_indices[sample_i] = source->index_of_column( pca::get_concatenated_sample_ids( &m_samples, sample_i )) ;
		if( sample_i > 0 ) {
			if( sample_column_indices[sample_i] <= sample_column_indices[sample_i-1] ) {
				m_ui_context.logger() << "!! Error in PCAComputer::load_matrix(): column order does not match samples.\n" ;
				throw genfile::MalformedInputError( source->get_source_spec(), 0, sample_column_indices[ sample_i ] ) ;
			}
		}
	}

	// Now apply sample exclusions.
	for( std::size_t sample_i = 0; sample_i < m_samples.get_number_of_individuals(); ++sample_i ) {
		std::string id ;
		// find row corresponding to next sample.
		for( (*source) >> id; (*source) && id != pca::get_concatenated_sample_ids( &m_samples, sample_i ); (*source) >> statfile::ignore_all() >> id ) ;
		if( !(*source) || source->number_of_rows_read() == source->number_of_rows() ) {
			throw genfile::MalformedInputError( source->get_source_spec(), source->number_of_rows_read(), 0 ) ;
		}
		if( source->number_of_rows_read() != ( sample_column_indices[ sample_i ] - 1 ) ) {
			m_ui_context.logger()
				<< "!! Error( PCAComputer::load_matrix ): sample " << m_samples.get_entry( sample_i, "id_1" )
				<< " is on row "
				<< source->number_of_rows_read()
				<< " but column "
				<< sample_column_indices[ sample_i ]
				<< ".\n" ;
			throw genfile::MalformedInputError( source->get_source_spec(), source->number_of_rows_read(), 0 ) ;
		}
		for( std::size_t sample_j = 0; sample_j < m_samples.get_number_of_individuals(); ++sample_j ) {
			(*source)
				>> statfile::ignore( sample_column_indices[ sample_j ] - source->current_column() )
				>> (*matrix)( sample_i, sample_j )
			;
		}
		(*source) >> statfile::ignore_all() ;
	}
}

void PCAComputer::compute_PCA() {
	assert( m_options.check_if_option_was_supplied( "-PCAs" ) ) ;
	Eigen::MatrixXd kinship_eigendecomposition( m_number_of_samples, m_number_of_samples + 1 ) ;
#if HAVE_LAPACK
	if( m_options.check( "-no-lapack" ))
#endif
	{
		m_ui_context.logger() << "========================================================================\n" ;
		m_ui_context.logger() << "PCAComputer: Computing eigenvalue decomposition of kinship matrix using Eigen...\n" ;
		Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > solver( m_kinship_matrix ) ;
		if( solver.info() == Eigen::Success ) {
			m_ui_context.logger() << "PCAComputer: Done, writing results...\n" ;
		} else if( solver.info() == Eigen::NumericalIssue ) {
			m_ui_context.logger() << "PCAComputer: Oh dear, numerical issue, writing results...\n" ;
		} else if( solver.info() == Eigen::NoConvergence ) {
			m_ui_context.logger() << "PCAComputer: Oh dear, no convergence, writing results...\n" ;
		} else {
			m_ui_context.logger() << "PCAComputer: Oh dear, unknown error, writing results...\n" ;
		}
		kinship_eigendecomposition.block( 0, 0, m_number_of_samples, 1 ) = solver.eigenvalues().reverse() ;
		kinship_eigendecomposition.block( 0, 1, m_number_of_samples, m_number_of_samples ) = Eigen::Reverse< Eigen::MatrixXd, Eigen::Horizontal >( solver.eigenvectors() ) ;
	}
#if HAVE_LAPACK
	else { // -no-lapack not specified.
		Eigen::VectorXd eigenvalues( m_number_of_samples ) ;
		Eigen::MatrixXd eigenvectors( m_number_of_samples, m_number_of_samples ) ;
		m_ui_context.logger() << "PCAComputer: Computing eigenvalue decomposition of kinship matrix using lapack...\n" ;
		lapack::compute_eigendecomposition( m_kinship_matrix, &eigenvalues, &eigenvectors ) ;
		//lapack::compute_partial_eigendecomposition( m_kinship_matrix, &eigenvalues, &eigenvectors, m_number_of_PCAs_to_compute ) ;
		kinship_eigendecomposition.block( 0, 0, m_number_of_samples, 1 )  = eigenvalues.reverse() ;
		kinship_eigendecomposition.block( 0, 1, m_number_of_samples, m_number_of_samples ) = Eigen::Reverse< Eigen::MatrixXd, Eigen::Horizontal >( eigenvectors ) ;
	}
#endif

	// Verify the decomposition.
	{
		std::size_t size = std::min( std::size_t(10), m_number_of_samples ) ;
		m_ui_context.logger() << "Top-left of U U^t is:\n" ;
		{
			
			Eigen::MatrixXd UUT = kinship_eigendecomposition.block( 0, 1, size, m_number_of_samples ) ;
			UUT *= kinship_eigendecomposition.block( 0, 1, m_number_of_samples, size ).transpose() ;
			write_matrix_to_stream( m_ui_context.logger(), UUT ) ;
		}
		m_ui_context.logger() << "Top-left of original kinship matrix is:\n" ;
		write_matrix_to_stream( m_ui_context.logger(), m_kinship_matrix.block( 0, 0, size, size )) ;
		
		Eigen::VectorXd v = kinship_eigendecomposition.block( 0, 0, m_number_of_samples, 1 ) ;
		Eigen::MatrixXd reconstructed_kinship_matrix = kinship_eigendecomposition.block( 0, 1, size, m_number_of_samples ) ;
		reconstructed_kinship_matrix *= v.asDiagonal() ;
		reconstructed_kinship_matrix *= kinship_eigendecomposition.block( 0, 1, m_number_of_samples, size ).transpose() ;
		// Knock out the upper triangle for comparison with the original kinship matrix.
		for( std::size_t i = 0; i < size; ++i ) {
			for( std::size_t j = i+1; j < size; ++j ) {
				reconstructed_kinship_matrix(i,j) = std::numeric_limits< double >::quiet_NaN() ;
			}
		}

		m_ui_context.logger() << "Top-left of reconstructed kinship matrix is:\n" ;
		write_matrix_to_stream( m_ui_context.logger(), reconstructed_kinship_matrix.block( 0, 0, size, size )) ;

		m_ui_context.logger() << "Verifying the decomposition...\n" ;
		double diff = 0.0 ;
		for( std::size_t i = 0; i < size; ++i ) {
			for( std::size_t j = 0; j <= i; ++j ) {
				diff = std::max( diff, std::abs( reconstructed_kinship_matrix(i,j) - m_kinship_matrix( i, j ))) ;
			}
		}
		m_ui_context.logger() << "...maximum discrepancy in lower diagonal is " << diff << ".\n" ;
		if( diff != diff ) {
			m_ui_context.logger() << "...yikes, there were NaNs.\n" ;
		} else if( diff > 0.01 ) {
			m_ui_context.logger() << "...warning: diff > 0.01.\n" ;
		}
		
		
	}
	
	using genfile::string_utils::to_string ;

	send_UDUT(
		"Number of SNPs: " + to_string( m_number_of_snps ) + "\n" +
			"Number of samples: " + to_string( m_number_of_samples ),
		m_number_of_snps,
		kinship_eigendecomposition,
		boost::function< genfile::VariantEntry ( int ) >(),
		boost::bind(
			&get_eigendecomposition_header,
			"v",
			_1
		)
	) ;

	if( m_number_of_PCAs_to_compute > 0 ) {
		//
		// Let X  be the L\times n matrix (L SNPs, n samples) of (mean-centred, scaled) genotypes.  
		// We want the projection of columns of X onto the matrix S whose columns are unit eigenvectors of (1/L) X X^t
		// where L is the number of SNPs.
		// to compute the row of the matrix S of unit eigenvectors of X X^t that corresponds to the current SNP.
		// The matrix S is given by
		//               1 
		//       S = ------- X U D^{-1/2}
		//           \sqrt(L)
		// where
		//       (1/L) X^t X = U D U^t
		// is the eigenvalue decomposition of (1/L) X^t X that we are passed in via set_UDUT (and L is the number of SNPs).
		//
		// The projections are then given by S^t X = \sqrt(L) U D^{1/2}
		// We also wish to compute the correlation between the SNP and the PCA component.
		// With S as above, the PCA components are the projections of columns of X onto columns of S.
		// If we want samples to correspond to columns, this is
		//   S^t X 
		// which can be re-written
		//   sqrt(L) U D^{1/2}
		// i.e. we may as well compute the correlation with columns of U.
		
		m_number_of_PCAs_to_compute = std::min( m_number_of_PCAs_to_compute, m_number_of_samples ) ;
		m_number_of_PCAs_to_compute = std::min( m_number_of_PCAs_to_compute, m_number_of_snps ) ;

		// PCAs are given by U D^{1/2}.
		Eigen::VectorXd PCA_eigenvalues = kinship_eigendecomposition.block( 0, 0, m_number_of_PCAs_to_compute, 1 ) ;
		Eigen::VectorXd v = PCA_eigenvalues.array().sqrt() ;
		Eigen::MatrixXd PCAs
			= kinship_eigendecomposition.block( 0, 1, m_number_of_samples, m_number_of_PCAs_to_compute ) * v.asDiagonal() ;

		send_PCAs(
			"Number of SNPs: " + to_string( m_number_of_snps ) + "\n" +
			"Number of samples: " + to_string( m_number_of_samples ) + "\n" +
			"Note: these PCAs are 1/sqrt(L) times the projection of samples onto unit eigenvectors of the variance-covariance matrix\n"
			"    1/(L-1) X X^t,\n"
			"where X is the L x N matrix of genotypes, L is the number of SNPs, and N the number of samples.\n"
			"The constant 1/sqrt(L) ensures that the PCAs do not grow with the number of SNPs.",
			PCA_eigenvalues,
			PCAs,
			boost::bind(
				&pca::get_concatenated_sample_ids,
				&m_samples,
				_1
			),
			boost::bind(
				&pca::string_and_number,
				"PCA_",
				_1
			)
		) ;
	}
	
	m_ui_context.logger() << "========================================================================\n" ;
}

void PCAComputer::send_UDUT_to( UDUTCallback callback ) {
	m_UDUT_signal.connect( callback ) ;
}

void PCAComputer::send_PCAs_to( PCACallback callback ) {
	m_PCA_signal.connect( callback ) ;
}

void PCAComputer::send_UDUT( std::string description, std::size_t number_of_snps, Eigen::MatrixXd const& UDUT, GetNames row_names, GetNames column_names ) {
	m_UDUT_signal( description, number_of_snps, UDUT, row_names, column_names ) ;
}

void PCAComputer::send_PCAs( std::string description, Eigen::VectorXd const& eigenvalues, Eigen::MatrixXd const& PCAs, GetNames pca_row_names, GetNames pca_column_names ) {
	m_PCA_signal( description, eigenvalues, PCAs, pca_row_names, pca_column_names ) ;
}

