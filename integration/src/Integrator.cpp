
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "integration/Integrator.hpp"

namespace integration {
	Integrator::Integrator( double desired_error ):
		m_desired_error( desired_error )
	{}
}
