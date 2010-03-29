#include <stddef.h>
#include <utility>
#include "genfile/GenomePosition.hpp"

namespace genfile {
	bool operator<( GenomePosition const& left, GenomePosition const& right ) {
		return left.m_data < right.m_data ;
	}

	bool operator<=( GenomePosition const& left, GenomePosition const& right ) {
		return left.m_data <= right.m_data ;
	}

	bool operator==( GenomePosition const& left, GenomePosition const& right ) {
		return left.m_data == right.m_data ;
	}

	bool operator!=( GenomePosition const& left, GenomePosition const& right ) {
		return left.m_data != right.m_data ;
	}

	std::ostream& operator<<( std::ostream& oStream, GenomePosition const& pos ) {
		return oStream << "(" << pos.chromosome() << ":" << pos.position() << ")" ;
	}

}
