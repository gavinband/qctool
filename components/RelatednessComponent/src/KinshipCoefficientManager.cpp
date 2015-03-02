
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include "Eigen/Core"
#include "components/RelatednessComponent/KinshipCoefficientManager.hpp"

void KinshipCoefficientManager::send_results_to( ResultsCallback callback ) {
	m_result_signal.connect( callback ) ;
}

void KinshipCoefficientManager::send_results(
	std::size_t const number_of_snps,
	Eigen::MatrixXd const& matrix1,
	Eigen::MatrixXd const& matrix2,
	std::string const& source,
	std::string const& description
) {
	m_result_signal( number_of_snps, matrix1, matrix2, source, description ) ;
}
