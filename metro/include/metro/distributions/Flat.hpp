
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_DISTRIBUTIONS_FLAT_HPP
#define METRO_DISTRIBUTIONS_FLAT_HPP

#include <memory>
#include "metro/UnivariateLogDensity.hpp"

namespace metro {
	namespace distributions {
		struct Flat: public UnivariateLogDensity {
		public:
			static UniquePtr create() ;
			
		public:
			Flat() {} ;
			Flat( Flat const& other ) {}
			Flat& operator=( Flat const& other ) { return *this ; }

		public:
			void evaluate_at( double ) {} ;

			double get_value_of_function() const {
				return 0.0 ;
			}

			double get_value_of_first_derivative() const {
				return 0.0 ;
			}

			double get_value_of_second_derivative() const {
				return 0.0 ;
			}
			
			std::string get_summary() const ;

		private:
		} ;
	}
}

#endif
