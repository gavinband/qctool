#ifndef FPUTILS_NAN_HPP
#define FPUTILS_NAN_HPP

namespace fputils {
	bool is_NaN( double value ) {
		return value != value ;
	}
}

#endif
