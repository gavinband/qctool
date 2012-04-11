
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
	std::string const& name,
	Eigen::MatrixXd const& matrix,
	std::string const& source,
	std::string const& description,
	GetNames get_row_names,
	GetNames get_column_names
) {
	m_result_signal( name, matrix, source, description, get_row_names, get_column_names ) ;
}
