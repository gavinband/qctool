
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef COMPONENTS_RELATEDNESS_COMPONENT_KINSHIPCOEFFICIENTCOMPUTATION_HPP
#define COMPONENTS_RELATEDNESS_COMPONENT_KINSHIPCOEFFICIENTCOMPUTATION_HPP

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/signals2.hpp>
#include "Eigen/Core"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "worker/Worker.hpp"
#include "worker/Task.hpp"
#include "appcontext/UIContext.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/FileUtil.hpp"
#include "components/RelatednessComponent/KinshipCoefficientManager.hpp"
//#include "components/RelatednessComponent/KinshipCoefficientBlockTask.hpp"

namespace impl {
	struct KinshipCoefficientComputerTask: public worker::Task {
		typedef std::auto_ptr< KinshipCoefficientComputerTask > UniquePtr ;
		typedef boost::function< void ( Eigen::VectorXd* data, Eigen::MatrixXd* result, double const scale ) > AccumulateXXt ;
		
		virtual void add_snp(
			genfile::VariantDataReader::SharedPtr data_reader,
			Eigen::VectorXd& genotypes,
			Eigen::VectorXd& genotype_non_missingness
		) = 0 ;
		
		virtual void finalise() = 0 ;
		virtual bool is_finalised() const = 0 ;
		virtual std::size_t number_of_snps() const = 0 ;
	} ;
}


struct KinshipCoefficientComputer: public KinshipCoefficientManager, public genfile::SNPDataSourceProcessor::Callback
{
public:
	typedef std::auto_ptr< KinshipCoefficientComputer > UniquePtr ;
	typedef boost::function<
		impl::KinshipCoefficientComputerTask::UniquePtr ( std::size_t, Eigen::MatrixXd*, Eigen::MatrixXd*, impl::KinshipCoefficientComputerTask::AccumulateXXt, bool )
	> ComputationFactory ;
	
	static impl::KinshipCoefficientComputerTask::UniquePtr compute_kinship( std::size_t const, Eigen::MatrixXd*, Eigen::MatrixXd*, impl::KinshipCoefficientComputerTask::AccumulateXXt, bool const ) ;
	static impl::KinshipCoefficientComputerTask::UniquePtr compute_concordance( std::size_t const, Eigen::MatrixXd*, Eigen::MatrixXd*, impl::KinshipCoefficientComputerTask::AccumulateXXt, bool const ) ;
	static impl::KinshipCoefficientComputerTask::UniquePtr compute_intensity_covariance( std::size_t const, Eigen::MatrixXd*, Eigen::MatrixXd*, impl::KinshipCoefficientComputerTask::AccumulateXXt, bool const ) ;

public:
	~KinshipCoefficientComputer() throw() {}

	KinshipCoefficientComputer(
		appcontext::OptionProcessor const& options,
		genfile::CohortIndividualSource const& samples,
		worker::Worker* worker,
		appcontext::UIContext& ui_context,
		ComputationFactory computation_factory
	) ;

	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader::SharedPtr data_reader ) ;
	void end_processing_snps() ;

private:
	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps ;
	std::size_t m_number_of_snps_processed ;
	genfile::CohortIndividualSource const& m_samples ;
	worker::Worker* m_worker ;
	boost::ptr_vector< impl::KinshipCoefficientComputerTask > m_tasks ;
	std::vector< Eigen::MatrixXd > m_result ;
	std::vector< Eigen::MatrixXd > m_non_missing_count ;
	std::vector< Eigen::VectorXd > m_genotypes ;
	std::vector< Eigen::VectorXd > m_genotype_non_missingness ;
	std::size_t m_number_of_tasks ;
	std::size_t m_number_of_snps_per_task ;
	std::size_t m_current_task ;
	ComputationFactory m_computation_factory ;
	void (*m_accumulate_xxt)( Eigen::VectorXd*, Eigen::MatrixXd*, double const ) ;
} ;

#endif
