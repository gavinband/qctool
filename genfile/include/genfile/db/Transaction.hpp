
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_DB_TRANSACTION_HPP
#define GENFILE_DB_TRANSACTION_HPP

#include <memory>

namespace genfile {
	namespace db {
		class Transaction {
		public:
			typedef std::auto_ptr< Transaction > UniquePtr ;
			virtual ~Transaction() {}
		} ;
	}
}

#endif
