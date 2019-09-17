
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_REGRESSION_THREADED_LOG_LIKELIHOOD_HPP
#define METRO_REGRESSION_THREADED_LOG_LIKELIHOOD_HPP

#include <vector>
#include <memory>
#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "Eigen/Core"
#include "metro/regression/LogLikelihood.hpp"
#include "metro/regression/Design.hpp"
#include "boost/threadpool.hpp"
#include "metro/concurrency/threadpool.hpp"

namespace metro {
	namespace regression {
		/*
		* Base class for classes implementing regression log-likelihoods
		*/
		struct ThreadedLogLikelihood: public LogLikelihood
		{
		public:
			typedef boost::function< LogLikelihood::UniquePtr ( regression::Design&, std::vector< metro::SampleRange > ) > CreateLL ;
			//typedef std::auto_ptr< ThreadedLogLikelihood > UniquePtr ;

			static UniquePtr create( Design& design, CreateLL, int number_of_threads ) ;
			static UniquePtr create( Design::UniquePtr design, CreateLL, int number_of_threads ) ;
			static UniquePtr create( Design& design, CreateLL, metro::concurrency::threadpool& pool ) ;
			static UniquePtr create( Design::UniquePtr design, CreateLL, metro::concurrency::threadpool& pool ) ;

		public:
			~ThreadedLogLikelihood() ;
			ThreadedLogLikelihood( Design& design, CreateLL, int number_of_threads ) ;
			ThreadedLogLikelihood( Design::UniquePtr design, CreateLL, int number_of_threads ) ;
			ThreadedLogLikelihood( Design& design, CreateLL, metro::concurrency::threadpool& pool ) ;
			ThreadedLogLikelihood( Design::UniquePtr design, CreateLL, metro::concurrency::threadpool& pool ) ;

			regression::Design& design() const { return *m_design ; } ;
			std::string get_parameter_name( std::size_t i ) const ;
			std::vector< std::string > get_parameter_names() const ;

			// Return a lx2 matrix identifying the l parameters.
			// For each parameter the row contains:
			// column 0: the index of the outcome column.
			// column 1: the index of the design matrix column.
			int number_of_parameters() const ;
			int number_of_outcomes() const ;
			IntegerMatrix identify_parameters() const ;

			void evaluate_at( Point const& parameters, int const numberOfDerivatives = 2 ) ;
			void evaluate( int const numberOfDerivatives = 2 ) ;

			Point const& parameters() const ;
			double get_value_of_function() const ;
			Vector get_value_of_first_derivative() const ;
			Matrix get_value_of_second_derivative() const ;
			
			std::string get_summary() const ;
		private:
			regression::Design* m_design ;
			bool const m_design_owned ;
			metro::concurrency::threadpool* m_pool ;
			bool m_pool_owned ;

			boost::ptr_vector< metro::regression::LogLikelihood > m_lls ;

		private:
			void create_lls( CreateLL create_ll ) ;
		} ;
	}
}

#endif