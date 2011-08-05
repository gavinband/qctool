#ifndef KinshipCoefficientComputation2_HPP
#define KinshipCoefficientComputation2_HPP

#include <boost/ptr_container/ptr_vector.hpp>
#include "Eigen/Core"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "worker/Worker.hpp"
#include "worker/Task.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"

namespace impl {
	struct KinshipCoefficientComputerTask ;

	struct KinshipCoefficientComputerTask: public worker::Task {
		KinshipCoefficientComputerTask(
			std::size_t number_of_samples,
			Eigen::MatrixXd* result,
			Eigen::MatrixXd* missing_count
		) ;

		void add_snp(
			genfile::SNPIdentifyingData const& id_data,
			genfile::VariantDataReader& data_reader
		) ;
		
		void finalise() { m_finalised = true ;}		
		bool is_finalised() const { return m_finalised ;}

		std::size_t number_of_snps() const { return m_id_data.size() ; }

		void operator()() ;
		
	private:
		std::size_t const m_number_of_samples ;
		Eigen::MatrixXd* m_result ;
		Eigen::MatrixXd* m_non_missing_count ;
		double const m_threshhold ;
		std::vector< genfile::SNPIdentifyingData > m_id_data ;
		std::vector< genfile::SingleSNPGenotypeProbabilities > m_genotypes ;
		
		Eigen::VectorXd m_data ;
		Eigen::VectorXd m_non_missingness_matrix ;
		bool m_finalised ;
	} ;
	
}

struct KinshipCoefficientComputer: public genfile::SNPDataSourceProcessor::Callback
{
public:
	KinshipCoefficientComputer(
		std::string const& filename,
		genfile::CohortIndividualSource const& samples,
		worker::Worker* worker
	) ;

	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

private:
	void write_output() ;
private:
	std::string const m_filename ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps ;
	std::size_t m_number_of_snps_processed ;
	genfile::CohortIndividualSource const& m_samples ;
	worker::Worker* m_worker ;
	boost::ptr_vector< impl::KinshipCoefficientComputerTask > m_tasks ;
	std::vector< Eigen::MatrixXd > m_result ;
	std::vector< Eigen::MatrixXd > m_non_missing_count ;
	std::size_t m_number_of_tasks ;
	std::size_t m_number_of_snps_per_task ;
	std::size_t m_current_task ;
} ;

#endif
