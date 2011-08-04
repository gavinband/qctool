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

	        void read_snp_identifying_data(
	            std::istream& aStream,
	            std::string* SNPID,
	            std::string* RSID,
				Chromosome* chromosome,
	            uint32_t* SNP_position,
	            char* first_allele,
	            char* second_allele
			) {
				std::string chromosome_string ;
				aStream >> *SNPID >> *RSID >> chromosome_string >> *SNP_position >> *first_allele >> *second_allele ;
				*chromosome = Chromosome( chromosome_string ) ;
			}
		}
	}
}
