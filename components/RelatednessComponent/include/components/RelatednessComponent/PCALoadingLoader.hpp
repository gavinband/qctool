
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_RELATEDNESS_COMPONENT_LOADINGS_LOADER_HPP
#define QCTOOL_RELATEDNESS_COMPONENT_LOADINGS_LOADER_HPP

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/VariantIdentifyingData.hpp"

namespace pca {
	struct PCALoadingLoader: public boost::noncopyable {
		PCALoadingLoader( genfile::CohortIndividualSource const& samples ) ;

		typedef boost::signals2::signal< void (
			std::vector< genfile::VariantIdentifyingData > const&,
			Eigen::VectorXd const&,
			Eigen::VectorXd const&,
			Eigen::MatrixXd const&,
			std::size_t number_of_snps,
			std::size_t number_of_loadings,
			std::vector< std::string > )
		> ResultSignal ;
		void send_loadings_to( ResultSignal::slot_type ) ;
		void load_loadings( std::string const& filename ) const ;
	private:
		genfile::CohortIndividualSource const& m_samples ;
		ResultSignal m_result_signal ;

		void load_loadings_impl(
			std::string const& filename,
			std::vector< genfile::VariantIdentifyingData >* snps,
			Eigen::VectorXd* counts,
			Eigen::VectorXd* frequencies,
			Eigen::MatrixXd* loadings,
			std::vector< std::string >* column_names
		) const ;
	} ;
}

#endif
