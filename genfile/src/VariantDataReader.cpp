#include <vector>
#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"

namespace genfile {
	VariantDataReader& VariantDataReader::get( std::string const& spec, std::vector< std::vector< Entry > >& data ) {
		vcf::VectorSetter setter( data ) ;
		get( spec, setter ) ;
		return *this ;
	}

	VariantDataReader& VariantDataReader::get( std::string const& spec, genfile::SingleSNPGenotypeProbabilities& data ) {
		vcf::GenotypeSetter< genfile::SingleSNPGenotypeProbabilities > setter( data ) ;
		get( spec, setter ) ;
		return *this ;
	}
}
