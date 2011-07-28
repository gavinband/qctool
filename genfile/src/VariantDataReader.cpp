#include "genfile/VariantDataReader.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"

namespace genfile {
	VariantDataReader& VariantDataReader::get( std::string const& spec, genfile::SingleSNPGenotypeProbabilities& data ) {
		vcf::GenotypeProbabilitySetter< genfile::GenotypeSetter< SingleSNPGenotypeProbabilities > > setter( set_genotypes( data ) ) ;
		get( spec, setter ) ;
		return *this ;
	}
}
