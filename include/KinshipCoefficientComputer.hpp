#ifndef KinshipCoefficientComputation2_HPP
#define KinshipCoefficientComputation2_HPP

#include <boost/ptr_container/ptr_vector.hpp>
#include "Eigen/Core"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "worker/Task.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"

namespace impl {
	struct KinshipCoefficientComputerTask ;

	struct KinshipCoefficientComputerTask: public worker::Task {
		KinshipCoefficientComputerTask(
			std::size_t number_of_samples,
			genfile::SNPIdentifyingData const& id_data,
			genfile::VariantDataReader& data_reader,
			Eigen::MatrixXd* result,
			Eigen::MatrixXd* missing_count
		):
			m_number_of_samples( number_of_samples ),
			m_id_data( id_data ),
			m_result( result ),
			m_non_missing_count( missing_count ),
			m_threshhold( 0.9 ),
			m_genotypes( m_number_of_samples )
		{
			data_reader.get( "genotypes", m_genotypes ) ;
		}

		void operator()() {
			Eigen::VectorXd data = Eigen::VectorXd::Constant( m_number_of_samples, std::numeric_limits< double >::quiet_NaN() ) ;
			Eigen::VectorXd non_missingness_matrix = Eigen::VectorXd::Zero( m_number_of_samples ) ;
			double allele_sum = 0.0 ;
			for( std::size_t sample_i = 0; sample_i < m_number_of_samples; ++sample_i ) {
				for( std::size_t g = 0; g < 3; ++g ) {
					if( m_genotypes( sample_i, g ) >= m_threshhold ) {
						data( sample_i ) = double( g ) ;
						non_missingness_matrix( sample_i ) = 1.0 ;
						allele_sum += g ;
						break ;
					}
				}
			}
			double const allele_freq = allele_sum / ( 2.0 * non_missingness_matrix.sum() ) ;
			for( std::size_t sample_i = 0; sample_i < m_number_of_samples; ++sample_i ) {
				if( non_missingness_matrix( sample_i ) ) {
					data( sample_i ) -= 2.0 * allele_freq ;
				}
				else {
					data( sample_i ) = 0.0 ; // this sample does not contribute for this SNP.
				}
			}
			m_result->noalias() += ( data * data.transpose() ) / ( 2.0 * allele_freq * ( 1.0 - allele_freq ) ) ;
			m_non_missing_count->noalias() += non_missingness_matrix * non_missingness_matrix.transpose() ;
		}
		
	private:
		std::size_t const m_number_of_samples ;
		genfile::SNPIdentifyingData const& m_id_data ;
		Eigen::MatrixXd* m_result ;
		Eigen::MatrixXd* m_non_missing_count ;
		double const m_threshhold ;
		genfile::SingleSNPGenotypeProbabilities m_genotypes ;
	} ;
	
}

struct KinshipCoefficientComputer: public genfile::SNPDataSourceProcessor::Callback
{
public:
	KinshipCoefficientComputer(
		std::string const& filename,
		genfile::CohortIndividualSource const& samples,
		worker::Worker* worker
	):
		m_filename( filename ),
		m_samples( samples ),
		m_worker( worker )
	{}

	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
		assert( m_samples.get_number_of_individuals() == number_of_samples ) ;
		m_number_of_samples = number_of_samples ;
		m_number_of_snps = number_of_snps ;
		m_current_task = 0 ;
		m_number_of_tasks = m_worker->get_number_of_worker_threads() ;
		m_result.resize( m_number_of_tasks ) ;
		m_non_missing_count.resize( m_number_of_tasks ) ;
		for( std::size_t i = 0; i < m_result.size(); ++i ) {
			m_result[i] = Eigen::MatrixXd::Zero( m_number_of_samples, m_number_of_samples ) ;
			m_non_missing_count[i] = Eigen::MatrixXd::Zero( m_number_of_samples, m_number_of_samples ) ;
		}
	}

	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& data_reader ) {
		if( m_tasks.size() <= m_current_task ) {
			assert( m_current_task == m_tasks.size() ) ;
			m_tasks.push_back(
				new impl::KinshipCoefficientComputerTask(
					m_number_of_samples,
					id_data,
					data_reader,
					&m_result[ m_current_task ],
					&m_non_missing_count[ m_current_task ]
				)
			) ;
			m_worker->tell_to_perform_task( m_tasks.back() ) ;
		}
		else {
			m_tasks[ m_current_task ].wait_until_complete() ;
			m_tasks.replace(
				m_current_task,
				new impl::KinshipCoefficientComputerTask(
					m_number_of_samples,
					id_data,
					data_reader,
					&m_result[ m_current_task ],
					&m_non_missing_count[ m_current_task ]
				)
			) ;
			m_worker->tell_to_perform_task( m_tasks[ m_current_task ] ) ;
		}
		m_current_task = ( m_current_task + 1 ) % m_number_of_tasks ;
	}
	
	void end_processing_snps() {
		for( std::size_t i = 1; i < m_result.size(); ++i ) {
			m_result[0].noalias() += m_result[i] ;
			m_non_missing_count[0].noalias() += m_non_missing_count[i] ;
		}
		m_result[0].array() /= m_non_missing_count[0].array() ;
		write_output() ;
	}

private:
	void write_output() {
		appcontext::OUTPUT_FILE_PTR file = appcontext::open_file_for_output( m_filename ) ;
		(*file) << "# Written by KinshipCoefficientComputer, " << appcontext::get_current_time_as_string() << "\n" ;
		(*file) << "# Description: " << "" << "\n" ;
		(*file) << "# Number of SNPs: " << m_number_of_snps << ".\n" ;
		(*file) << "id," ;
		for( std::size_t sample_i = 0; sample_i < m_number_of_samples; ++sample_i ) {
			if( sample_i > 0 ) {
				(*file) << "," ;
			}
			(*file) << m_samples.get_entry( sample_i, "id_1" ).as< std::string >() ;
		}
		(*file) << "\n" ;

		for( std::size_t sample_i = 0; sample_i < m_number_of_samples; ++sample_i ) {
			(*file) << m_samples.get_entry( sample_i, "id_1" ).as< std::string >() << "," ;
			for( std::size_t sample_j = 0; sample_j < m_number_of_samples; ++sample_j ) {
				if( sample_j > 0 ) {
					(*file) << "," ;
				}
				(*file) << m_result[0]( sample_i, sample_j ) ;
			}
			(*file) << "\n" ;
		}
	}

private:
	std::string const m_filename ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps ;
	genfile::CohortIndividualSource const& m_samples ;
	worker::Worker* m_worker ;
	boost::ptr_vector< impl::KinshipCoefficientComputerTask > m_tasks ;
	std::vector< Eigen::MatrixXd > m_result ;
	std::vector< Eigen::MatrixXd > m_non_missing_count ;
	std::size_t m_current_task ;
	std::size_t m_number_of_tasks ;
} ;

#endif
