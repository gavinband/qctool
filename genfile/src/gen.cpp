#include "genfile/gen.hpp"

namespace genfile {
	namespace gen {
		namespace impl {
	        void read_snp_identifying_data(
	            std::istream& aStream,
	            std::string* SNPID,
	            std::string* RSID,
	            uint32_t* SNP_position,
	            char* first_allele,
	            char* second_allele
			) {
				aStream >> *SNPID >> *RSID >> *SNP_position >> *first_allele >> *second_allele ;
			}
		}
	}
}
