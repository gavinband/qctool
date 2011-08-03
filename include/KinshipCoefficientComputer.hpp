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
			genfile::SNPIdentifyingData const& id_data,
			genfile::VariantDataReader& data_reader,
			Eigen::MatrixXd* result,
			Eigen::MatrixXd* missing_count
		) ;

		void operator()() ;
		
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
	genfile::CohortIndividualSource const& m_samples ;
	worker::Worker* m_worker ;
	boost::ptr_vector< impl::KinshipCoefficientComputerTask > m_tasks ;
	std::vector< Eigen::MatrixXd > m_result ;
	std::vector< Eigen::MatrixXd > m_non_missing_count ;
	std::size_t m_current_task ;
	std::size_t m_number_of_tasks ;
} ;

#endif
