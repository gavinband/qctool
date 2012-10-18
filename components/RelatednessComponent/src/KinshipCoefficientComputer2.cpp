
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/signals2.hpp>
#include "Eigen/Core"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "worker/Worker.hpp"
#include "worker/Task.hpp"
#include "appcontext/UIContext.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/FileUtil.hpp"
#include "components/RelatednessComponent/KinshipCoefficientManager.hpp"
#include "components/RelatednessComponent/KinshipCoefficientComputer2.hpp"
#include "components/RelatednessComponent/KinshipCoefficientBlockTask.hpp"
#include "components/RelatednessComponent/mean_centre_genotypes.hpp"

// #define DEBUG_KINSHIP_COEFFICIENT_COMPUTER2 1

namespace impl {
	// subdivide the lower triangle of the given matrix into square pieces of about the given size.
	std::vector< KinshipCoefficientComputer2::BlockExtent > subdivide_matrix_lower_triangle_into_squares(
		int const n,
		double square_size_requested
	) {
		std::vector< KinshipCoefficientComputer2::BlockExtent > result ;
		//
		// We wish to find the number a such that a (a-1) / 2 is close to the number of squares, n
		// That is
		// a ( a - 1 ) = 2 n
		// a^2 - a - 2n = 0
		// a = ( 1 + sqrt( 1 + 8n )) / 2
		//
		int const a = n / square_size_requested ;
		double square_size = double( n ) / a ;

		for( int i = 0; i < a; ++i ) {
			int const x = i * square_size ;
			for( int j = 0; j <= i; ++j ) {
				int const y = j * square_size ;
				result.push_back(
					KinshipCoefficientComputer2::BlockExtent(
						x, y,
						((i+1) == a) ? n : x + square_size,
						((j+1) == a) ? n : y + square_size 
					)
				) ;
			}
		}

		return result ;
	}
	
	std::vector< KinshipCoefficientComputer2::BlockExtent > subdivide_matrix_lower_triangle_into_blocks(
		int const n,
		int const off_diagonal_size
	) {
		// Make two diagonal blocks and a number of smaller blocks
		std::vector< KinshipCoefficientComputer2::BlockExtent > result ;
		int const k = n / 2;
		result.push_back( KinshipCoefficientComputer2::BlockExtent( 0, 0, k, k ) ) ;
		result.push_back( KinshipCoefficientComputer2::BlockExtent( k, k, n, n ) ) ;

		int const a = 2 ;
		double square_size = double( k ) / a ;
		
		for( int i = 0; i < a; ++i ) {
			int const x = k + ( i * square_size ) ;
			for( int j = 0; j < a; ++j ) {
				int const y = j * square_size ;
				result.push_back(
					KinshipCoefficientComputer2::BlockExtent(
						x, y,
						((i+1) == a) ? n : x + square_size,
						((j+1) == a) ? k : y + square_size 
					)
				) ;
			}
		}
		return result ;
	}
	
#if HAVE_CBLAS
	struct NormaliseGenotypesAndComputeXXtSymmetricBlockUsingCblasTask: public worker::Task {
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::MatrixXd Matrix ;
		typedef Eigen::Block< Matrix > MatrixBlock ;

		NormaliseGenotypesAndComputeXXtSymmetricBlockUsingCblasTask(
			MatrixBlock result,
			MatrixBlock non_missingness,
			Vector& genotypes,
			Vector& genotype_non_missingness,
			genfile::VariantDataReader::SharedPtr data_reader
		):
			m_result( result ),
			m_non_missingness( non_missingness ),
			m_genotypes( genotypes ),
			m_genotype_non_missingness( genotype_non_missingness ),
			m_data_reader( data_reader ),
			m_frequency_threshhold( 0.01 )
		{}
		
		void operator()() {
			genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( m_genotypes, m_genotype_non_missingness, 0.9, 0, 0, 1, 2 ) ;
			m_data_reader->get( "genotypes", setter ) ;
			double const allele_frequency = m_genotypes.sum() / ( 2.0 * m_genotype_non_missingness.sum() ) ;
			if( allele_frequency > m_frequency_threshhold ) {
				pca::mean_centre_genotypes( &m_genotypes, m_genotype_non_missingness, allele_frequency ) ;
				int const N = m_genotypes.size() ;
				cblas_dsyr(
					CblasColMajor,
					CblasLower,
					N,
					1.0 / ( 2.0 * allele_frequency * ( 1.0 - allele_frequency ) ),
					m_genotypes.data(),
					m_genotypes.innerStride(),
					m_result.data(),
					m_result.outerStride()
				) ;

				cblas_dsyr(
					CblasColMajor,
					CblasLower,
					N,
					1.0,
					m_genotype_non_missingness.data(),
					m_genotype_non_missingness.innerStride(),
					m_non_missingness.data(),
					m_non_missingness.outerStride()
				) ;
			}
		}

	private:
		MatrixBlock m_result ;
		MatrixBlock m_non_missingness ;
		Vector& m_genotypes ;
		Vector& m_genotype_non_missingness ;
		genfile::VariantDataReader::SharedPtr m_data_reader ;
		double const m_frequency_threshhold ;
	} ;
#endif

	struct NormaliseGenotypesAndComputeXXtSymmetricBlockUsingEigenTask: public worker::Task {
		typedef Eigen::VectorXd Vector ;
		typedef Eigen::VectorBlock< Vector > VectorBlock ;
		typedef Eigen::MatrixXd Matrix ;
		typedef Eigen::Block< Matrix > MatrixBlock ;

		NormaliseGenotypesAndComputeXXtSymmetricBlockUsingEigenTask(
			MatrixBlock result,
			MatrixBlock non_missingness,
			Vector& genotypes,
			Vector& genotype_non_missingness,
			genfile::VariantDataReader::SharedPtr data_reader
		):
			m_result( result ),
			m_non_missingness( non_missingness ),
			m_genotypes( genotypes ),
			m_genotype_non_missingness( genotype_non_missingness ),
			m_data_reader( data_reader ),
			m_frequency_threshhold( 0.01 )
		{}
		
		void operator()() {
			genfile::vcf::ThreshholdingGenotypeSetter< Eigen::VectorXd > setter( m_genotypes, m_genotype_non_missingness, 0.9, 0, 0, 1, 2 ) ;
			m_data_reader->get( "genotypes", setter ) ;
			double const allele_frequency = m_genotypes.sum() / ( 2.0 * m_genotype_non_missingness.sum() ) ;
			if( allele_frequency > m_frequency_threshhold ) {
				pca::mean_centre_genotypes( &m_genotypes, m_genotype_non_missingness, allele_frequency ) ;
				m_result
					.selfadjointView< Eigen::Lower >()
					.rankUpdate( m_genotypes, 1.0 / ( 2.0 * allele_frequency * ( 1.0 - allele_frequency ) ) ) ;
				m_non_missingness
					.selfadjointView< Eigen::Lower >()
					.rankUpdate( m_genotype_non_missingness, 1.0 ) ;
			}
		}

	private:
		MatrixBlock m_result ;
		MatrixBlock m_non_missingness ;
		Vector& m_genotypes ;
		Vector& m_genotype_non_missingness ;
		genfile::VariantDataReader::SharedPtr m_data_reader ;
		double const m_frequency_threshhold ;
	} ;
}

std::ostream& operator<<( std::ostream& out, KinshipCoefficientComputer2::BlockExtent const& block ) {
	return out << "[ " << block.x() << "-" << block.x_end() << " ), [ " << block.y() << "-" << block.y_end() << " )" ;
}


KinshipCoefficientComputer2::~KinshipCoefficientComputer2() throw() {}

KinshipCoefficientComputer2::KinshipCoefficientComputer2(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	worker::Worker* worker,
	appcontext::UIContext& ui_context
):
	m_options( options ),
	m_ui_context( ui_context ),
	m_worker( worker )
{	
}

void KinshipCoefficientComputer2::begin_processing_snps( std::size_t number_of_samples ) {
	m_result.resize( 7 ) ; //m_worker->get_number_of_worker_threads() + 1 ) ;
	m_non_missing_count.resize( m_result.size() ) ;
	for( std::size_t i = 0; i < m_result.size(); ++i ) {
		m_result[i].resize( number_of_samples, number_of_samples ) ;
		m_non_missing_count[i].resize( number_of_samples, number_of_samples ) ;
		m_result[i].setZero() ;
		m_non_missing_count[i].setZero() ;
	}

/*	m_subdivision.clear() ;
	if( 1 ) { //m_worker->get_number_of_worker_threads() == 1 ) {
		// Just do the whole matrix in one go, as this is faster.
		m_subdivision.push_back(
			BlockExtent( 0, 0, number_of_samples, number_of_samples )
		) ;
	}
	else {
		// Subdivide into two diagonal blocks (which run quickly) and a number of smaller off-diagonals.
		m_subdivision = impl::subdivide_matrix_lower_triangle_into_blocks(
			number_of_samples,
			std::sqrt( number_of_samples ) * 5
		) ;
	}
#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER2 
	std::cerr << "Subdivided matrix into " << m_subdivision.size() << "pieces.\n" ;
	for( std::size_t i = 0; i < m_subdivision.size(); ++i ) {
		std::cerr << m_subdivision[i] << "\n" ;
	}
#endif
*/
	m_number_of_snps_processed = 0 ;
	
	m_genotypes.resize( m_result.size() ) ;
	m_genotype_non_missingness.resize( m_result.size() ) ;
	m_data_index = 0 ;
}

void KinshipCoefficientComputer2::processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader::SharedPtr data_reader ) {
	Eigen::VectorXd& genotypes = m_genotypes[ m_data_index ] ;
	Eigen::VectorXd& genotype_non_missingness = m_genotype_non_missingness[ m_data_index ] ;
	Eigen::MatrixXd& result = m_result[ m_data_index ] ;
	Eigen::MatrixXd& non_missing_count = m_non_missing_count[ m_data_index ] ;
	std::size_t const task_index = m_data_index ;
	
	std::auto_ptr< worker::Task > task(
		new impl::NormaliseGenotypesAndComputeXXtSymmetricBlockUsingCblasTask(
			result.block( 0, 0, result.rows(), result.cols() ),
			non_missing_count.block( 0, 0, non_missing_count.rows(), non_missing_count.cols() ),
			genotypes, genotype_non_missingness,
			data_reader
		)
	) ;
	;
	if( task_index < m_tasks.size() ) {
		// we must not clobber the running task.
		m_tasks[ task_index ].wait_until_complete() ;
		m_tasks.replace( task_index, task ) ;
	} else {
		m_tasks.push_back( task ) ;
	}
	m_worker->tell_to_perform_task( m_tasks[ task_index ] ) ;
	m_data_index = ( m_data_index + 1 ) % m_result.size() ;
	++m_number_of_snps_processed ;
}

void KinshipCoefficientComputer2::end_processing_snps() {
	while( !m_tasks.empty() ) {
		m_tasks.front().wait_until_complete() ;
		m_tasks.pop_front() ;
	}

	for( std::size_t i = 1; i < m_result.size(); ++i ) {
		m_result[0] += m_result[i] ;
		m_non_missing_count[0] += m_non_missing_count[i] ;
	}

	m_result[0].array() /= m_non_missing_count[0].array() ;

#if DEBUG_KINSHIP_COEFFICIENT_COMPUTER2
	std::cerr << "result is:\n" << m_result[0].block( 0, 0, 10, 10 ) << ".\n" ;
	std::cerr << "non-missingness is:\n" << m_non_missing_count[0].block( 0, 0, 10, 10 ) << ".\n" ;
#endif
	
	std::string description = "Number of SNPs: "
		+ genfile::string_utils::to_string( m_number_of_snps_processed )
		+ "\nNumber of samples: "
		+ genfile::string_utils::to_string( m_result[0].cols() ) ;
	
	send_results(
		m_result[0],
		"KinshipCoefficientComputer2",
		description
	) ;
}
