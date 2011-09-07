#ifndef KINSHIP_COEFFICIENT_MANAGER_HPP
#define KINSHIP_COEFFICIENT_MANAGER_HPP

#include <string>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include <Eigen/Core>

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
	typedef boost::function< void( std::string const& name, genfile::SNPIdentifyingData const& snp, Eigen::VectorXd const&, GetNames ) > PerVariantResultsCallback ;

	void send_results_to( ResultsCallback callback ) ;
	void send_results( std::string const& name, Eigen::MatrixXd const&, std::string const& source, std::string const& description, GetNames, GetNames ) ;
	void send_per_variant_results_to( PerVariantResultsCallback callback ) ;
	void send_per_variant_results( std::string const& name, genfile::SNPIdentifyingData const& snp, Eigen::VectorXd const& data, GetNames ) ;

private:
	typedef boost::signals2::signal< void( std::string const&, Eigen::MatrixXd const&, std::string const&, std::string const&, GetNames, GetNames ) > ResultSignal ;
	ResultSignal m_result_signal ;
	typedef boost::signals2::signal< void( std::string const&, genfile::SNPIdentifyingData const&, Eigen::VectorXd const&, GetNames ) > PerVariantResultSignal ;
	PerVariantResultSignal m_per_variant_result_signal ;
} ;

#endif