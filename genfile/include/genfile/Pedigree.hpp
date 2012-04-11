
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_PEDIGREE_HPP
#define GENFILE_PEDIGREE_HPP

#include <memory>
#include <string>
#include <vector>

namespace genfile {
	class Pedigree
	{
	public:
		typedef std::auto_ptr< Pedigree > UniquePtr ;
		enum Sex {eMale = 0, eFemale = 1, eUnknown = 2 } ;
	public:
		static UniquePtr create( std::string const& pedigree_spec ) ;
		virtual ~Pedigree() {}

		// Return the number of individuals in this pedigree.
		virtual std::size_t get_number_of_individuals() const = 0 ;
		// Return the identifier of the ith individual, in lexicographic order
		virtual std::string get_id_of( std::size_t i ) const = 0 ;
		// Return the family of the given individual
		virtual std::string get_family_of( std::string const& individual ) const = 0 ;
		// Return the parents of the given individual, in lexicographic order of id
		virtual std::vector< std::string > get_parents_of( std::string const& individual ) const = 0 ;
		// Return the children of the given individual, in lexicographic order of id
		virtual std::vector< std::string > get_children_of( std::string const& individual ) const = 0 ;
		// Return the sex of the given individual
		virtual Sex get_sex_of( std::string const& individual ) const = 0 ;
	} ;
}

#endif
