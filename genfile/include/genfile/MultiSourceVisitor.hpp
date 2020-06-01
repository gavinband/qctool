
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_MULTISOURCEVISITOR_HPP
#define GENFILE_MULTISOURCEVISITOR_HPP 1

#include <vector>
#include "boost/optional.hpp"
#include "boost/noncopyable.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	struct MultiSourceVisitor: public boost::noncopyable {
	public:
		typedef boost::function<
			void (
				std::vector< int > const& changed,
				std::vector< genfile::VariantIdentifyingData > const& variants,
				std::vector< genfile::VariantDataReader::SharedPtr > const& readers
		) > Callback ;
	public:
		~MultiSourceVisitor() {} ;
		virtual bool step( Callback callback ) = 0 ;
		virtual boost::optional< std::size_t > count() const = 0 ;
	} ;
}

#endif
