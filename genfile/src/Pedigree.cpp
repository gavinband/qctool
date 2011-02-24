#include <memory>
#include <string>
#include <vector>
#include "genfile/Pedigree.hpp"
#include "genfile/FromPedFilePedigree.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	Pedigree::UniquePtr Pedigree::create( std::string const& pedigree_spec ) {
		if( !pedigree_spec.size() > 5 && pedigree_spec.substr(0, 5) == "file:" ) {
			throw BadArgumentError( "genfile::Pedigree::create()", "pedigree_spec = \"" + pedigree_spec + "\"." ) ;
		}
		std::string spec = pedigree_spec.substr( 5, pedigree_spec.size() ) ;
		return Pedigree::UniquePtr(
			new FromPedFilePedigree(
				string_utils::strip( pedigree_spec.substr( 5, pedigree_spec.size() ), " \t" )
			)
		) ;
	}
}
