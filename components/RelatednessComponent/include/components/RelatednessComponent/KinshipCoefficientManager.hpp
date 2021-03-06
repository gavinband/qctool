
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef COMPONENTS_RELATEDNESS_COMPONENT_KINSHIP_COEFFICIENT_MANAGER_HPP
#define COMPONENTS_RELATEDNESS_COMPONENT_KINSHIP_COEFFICIENT_MANAGER_HPP

#include <string>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include <Eigen/Core>
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantIdentifyingData.hpp"

struct KinshipCoefficientManager
{
public:
	typedef std::auto_ptr< KinshipCoefficientManager > UniquePtr ;
public:
	virtual ~KinshipCoefficientManager() throw() {}

	typedef boost::function< genfile::VariantEntry ( std::size_t ) > GetNames ;
	typedef boost::signals2::signal< void( std::size_t, Eigen::MatrixXd const&, Eigen::MatrixXd const&, std::string const&, std::string const& ) > ResultSignal ;
	typedef ResultSignal::slot_type ResultsCallback ;
	typedef boost::function< void( std::string const& name, genfile::VariantIdentifyingData const& snp, Eigen::VectorXd const&, GetNames ) > PerVariantResultsCallback ;

	void send_results_to( ResultsCallback callback ) ;
	void send_results( std::size_t const, Eigen::MatrixXd const&, Eigen::MatrixXd const&, std::string const& source, std::string const& description ) ;
private:
	ResultSignal m_result_signal ;
} ;

#endif
