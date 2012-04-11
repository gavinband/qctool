
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef STATFILE_BLOB_HPP
#define STATFILE_BLOB_HPP

#include <vector>

namespace statfile {
	struct Blob: private std::vector< unsigned char >
	{
	private:
		typedef std::vector< unsigned char > base_t ;
	public:
		Blob( std::size_t size ) ;
		std::size_t size() const { return base_t::size() ; }

		unsigned char const* begin() const { return &front() ; }
		unsigned char const* end() const { return (&front()) + size() ; }
		unsigned char* begin() { return &front() ; }
		unsigned char* end() { return (&front()) + size() ; }

		template< typename I > void assign( I begin, I const& end ) {
			base_t::assign( begin, end ) ;
		}
	} ;
}

#endif
