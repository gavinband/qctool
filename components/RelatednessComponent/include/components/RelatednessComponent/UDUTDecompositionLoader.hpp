#ifndef QCTOOL_RELATEDNESS_COMPONENT_UDUT_LOADER_HPP
#define QCTOOL_RELATEDNESS_COMPONENT_UDUT_LOADER_HPP

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <boost/signals2/signal.hpp>
#include <Eigen/Core>

namespace relatedness {
	struct UDUTDecompositionLoader: public boost::noncopyable {
		UDUTDecompositionLoader() ;
		
		typedef boost::function< void ( Eigen::MatrixXd const& ) > ResultCallback ;
		void send_results_to( ResultCallback ) ;
		void load_matrix( std::string const& filename ) const ;

	private:
		typedef boost::signals2::signal< void ( Eigen::MatrixXd const& ) > ResultSignal ;
		ResultSignal m_result_signal ;
	} ;
}

#endif
