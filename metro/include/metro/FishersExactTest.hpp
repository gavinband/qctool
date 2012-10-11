
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef METRO_FISHERS_EXACT_TEST
#define METRO_FISHERS_EXACT_TEST

#include <Eigen/Core>
#include <boost/math/distributions/hypergeometric.hpp>

namespace metro {
	//
	// Fisher's exact test between categories
	//               red           black
	//     good |     a      |       b        |
	//      bad |     c      |       d        |
	//
	// Under the null that both rows are distributed the same, a is hypergeometrically
	// distributed.  Fishers' p-value is then the mass under the hypergeometric distribution
	// of all values of a greater than or equal to the one given that are consistent with
	// the given margins of the table.
	struct FishersExactTest {
		FishersExactTest( Eigen::Matrix2d const& matrix ) ;
		double get_OR() const ;
		std::pair< double, double > get_confidence_interval() const ;

		enum Alternative { eLess = 0, eGreater = 1, eTwoSided = 2 } ;
		double get_pvalue( Alternative const = eGreater ) const ;

		private:
			Eigen::Matrix2d const m_matrix ;
			boost::math::hypergeometric m_distribution ;
	} ;
}

#endif
