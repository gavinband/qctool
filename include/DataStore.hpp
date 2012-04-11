
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef APPCONTEXT_DATASTORE_HPP
#define APPCONTEXT_DATASTORE_HPP

#include <memory>
#include <string>
#include "stdint.h"
#include <boost/function.hpp>
#include "genfile/SNPIdentifyingData.hpp"

// Model of a store of data.
struct DataStore {
public:
	typedef std::auto_ptr< DataStore > UniquePtr ;
	static UniquePtr create( std::string const& spec ) ;
public:
	virtual ~DataStore() throw() {}
	typedef int64_t SnpId ;
	typedef int64_t EntityId ;
	virtual std::string get_spec() const = 0 ;
	virtual SnpId get_or_create_SNP( genfile::SNPIdentifyingData const& ) = 0 ;
	virtual SnpId get_or_create_entity( std::string name, std::string description ) = 0 ;
	virtual void set_relationship( std::string const& left, std::string const& relation, std::string const& right ) const = 0 ;
	virtual void store_per_variant_data( SnpId snp_id, int64_t field_id, int64_t cohort_id, int64_t storage_id, char const* buffer, char const* const end ) = 0 ;
	virtual void get_entities_by_relation( std::string const& relationship, std::string const& related_entity, boost::function< void ( db::Connection::RowId, std::string const& ) > callback ) = 0 ;
	struct Transaction {
		typedef std::auto_ptr< Transaction > UniquePtr ;
		virtual ~Transaction() {}
	} ;
	typedef Transaction::UniquePtr TransactionPtr ;
	virtual TransactionPtr open_transaction() = 0 ;
} ;

#endif
