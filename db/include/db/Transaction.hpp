#ifndef DB_TRANSACTION_HPP
#define DB_TRANSACTION_HPP

#include <memory>

namespace db {
	class Transaction {
	public:
		typedef std::auto_ptr< Transaction > UniquePtr ;
		virtual ~Transaction() {}
	} ;
}

#endif
