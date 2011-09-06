#ifndef KinshipCoefficientComputation2_HPP
#define KinshipCoefficientComputation2_HPP

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

struct KinshipCoefficientManager: public genfile::SNPDataSourceProcessor::Callback
{
public:
	static void declare_options( appcontext::OptionProcessor& options ) ;
	typedef std::auto_ptr< KinshipCoefficientManager > UniquePtr ;
	static UniquePtr create(
		appcontext::OptionProcessor const& options,
		genfile::CohortIndividualSource const& samples,
		worker::Worker* worker,
		appcontext::UIContext& ui_context
	) ;

public:
	virtual ~KinshipCoefficientManager() throw() {}

	typedef boost::function< genfile::VariantEntry ( std::size_t ) > GetNames ;
	typedef boost::function< void( std::string const& name, Eigen::MatrixXd const&, std::string const& source, std::string const& description, GetNames, GetNames ) > ResultsCallback ;

	void send_results_to( ResultsCallback callback ) ;
	void send_results( std::string const& name, Eigen::MatrixXd const&, std::string const& source, std::string const& description, GetNames, GetNames ) ;

private:
	typedef boost::signals2::signal< void( std::string const& name, Eigen::MatrixXd const&, std::string const& source, std::string const& description, GetNames, GetNames ) > ResultSignal ;
	ResultSignal m_result_signal ;
} ;

struct KinshipCoefficientComputer: public KinshipCoefficientManager
{
public:
	~KinshipCoefficientComputer() throw() {}

	KinshipCoefficientComputer(
		appcontext::OptionProcessor const& options,
		genfile::CohortIndividualSource const& samples,
		worker::Worker* worker,
		appcontext::UIContext& ui_context
	) ;

	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

private:
	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ;
	std::string m_filename ;
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

struct PCAComputer: public KinshipCoefficientManager
{
	~PCAComputer() throw() {}
	PCAComputer(
		appcontext::OptionProcessor const& options,
		genfile::CohortIndividualSource const& samples,
		worker::Worker* worker,
		appcontext::UIContext& ui_context
	) ;

	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& ) {}
	void end_processing_snps() ;
private:
	appcontext::OptionProcessor const& m_options ;
	appcontext::UIContext& m_ui_context ;
	genfile::CohortIndividualSource const& m_samples ;
	std::string m_filename ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps ;
	Eigen::MatrixXd m_matrix ;
	
	void load_matrix( std::string const& filename, Eigen::MatrixXd* matrix ) const ;
} ;


#endif
