
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_RELATEDNESS_COMPONENT_UDUT_LOADER_HPP
#define QCTOOL_RELATEDNESS_COMPONENT_UDUT_LOADER_HPP

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"

namespace relatedness {
	struct UDUTDecompositionLoader: public boost::noncopyable {
		UDUTDecompositionLoader( genfile::CohortIndividualSource const& samples ) ;
		
		typedef boost::function< void ( Eigen::MatrixXd const&, std::size_t number_of_samples, std::size_t number_of_snps ) > ResultCallback ;
		void send_UDUT_to( ResultCallback ) ;
		
		void load_matrix( std::string const& filename ) const ;
	private:
		genfile::CohortIndividualSource const& m_samples ;
		typedef boost::signals2::signal< void ( Eigen::MatrixXd const&, std::size_t number_of_samples, std::size_t number_of_snps ) > ResultSignal ;
		ResultSignal m_result_signal ;

		void load_matrix_impl( std::string const& filename, Eigen::MatrixXd* matrix, std::size_t* number_of_snps ) const ;
	} ;
}

#endif
