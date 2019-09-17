
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include "metro/regression/LogLikelihood.hpp"
#include "metro/SmoothFunction.hpp"
#include "metro/regression/LogUnnormalisedPosterior.hpp"

namespace metro {
	namespace regression {
		LogUnnormalisedPosterior::LogUnnormalisedPosterior( LogLikelihood& loglikelihood, SmoothFunction& posterior ):
			m_ll( loglikelihood ),
			m_prior( posterior )
		{}

		LogLikelihood const& LogUnnormalisedPosterior::ll() const {
			return m_ll ;
		}

		SmoothFunction const& LogUnnormalisedPosterior::prior() const {
			return m_prior ;
		}
		
		std::string LogUnnormalisedPosterior::get_summary() const {
			std::string result = "LogUnnormalisedPosterior:\n" ;
			result +=
				    "  loglikelihood: " + m_ll.get_summary()
				+ "\n  prior: " + m_prior.get_summary()
				+ "\n" ;
			return result ;
		}
		
		void LogUnnormalisedPosterior::evaluate_at(
			LogUnnormalisedPosterior::Vector const& parameters,
			int const numberOfDerivatives
		) {
			m_parameters = parameters ;
			evaluate( numberOfDerivatives ) ;
		}

		void LogUnnormalisedPosterior::evaluate(
			int const numberOfDerivatives
		) {
			m_ll.evaluate_at( m_parameters, numberOfDerivatives ) ;
			m_prior.evaluate_at( m_parameters, numberOfDerivatives ) ;
		}
	}
}

