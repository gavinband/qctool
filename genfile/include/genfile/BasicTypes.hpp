#ifndef GENFILE_BASICTYPES_HPP
#define GENFILE_BASICTYPES_HPP

#include <stdint.h>

namespace genfile {
	typedef int64_t Integer ;
	namespace Eigen {
		class MatrixXd ;
		class VectorXd ;
	}
	typedef Eigen::MatrixXd Matrix ;
	typedef Eigen::VectorXd Vector ;
}

#endif
