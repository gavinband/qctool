#include "genfile/gen.hpp"
#include "genfile/Chromosome.hpp"

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
				Chromosome* chromosome,
	            std::string* SNPID,
	            std::string* RSID,
	            uint32_t* SNP_position,
	            char* first_allele,
	            char* second_allele
			) {
				aStream >> *chromosome >> *SNPID >> *RSID >> *SNP_position >> *first_allele >> *second_allele ;
			}
		}

		uint32_t count_snp_blocks(
			std::istream& aStream
		) {
			uint32_t number_of_snp_blocks = 0 ;
			std::vector< char > buffer( 1024*1024 ) ;
			do {
				aStream.read( &(buffer[0]), 1024*1024 ) ;
				number_of_snp_blocks += std::count( buffer.begin(), buffer.begin() + aStream.gcount(), '\n' ) ;
				// A gen file can't contain a blank line.
				// Because popular editors (vim, nano, ..., but not emacs) typically add a trailing newline,
				// we might get in the situation where the GEN file has two trailing newlines thus messing
				// with our count.
				// Therefore we check here for the special case where what we've read ends in two newlines.
				if( (aStream.gcount() > 1) && (buffer[ aStream.gcount() - 1] == '\n') && (buffer[ aStream.gcount() - 2] == '\n') ) {
					throw FileHasTwoTrailingNewlinesError( "(unknown)", number_of_snp_blocks ) ;
				}
			}
			while( aStream ) ;

			// Most editors (vim, nano, but not emacs) automatically add a newline to the end of the file.
			// If the file has a trailing newline, we already have the correct count.
			// But if not, we've undercounted by one.
			if( aStream.gcount() > 0 ) {
				std::size_t pos = aStream.gcount() - 1 ;
				if( buffer[pos] != '\n' ) {
					++number_of_snp_blocks ;
				}
			}

			// We should have reached eof.
			// If so, report the data now.
			if( !aStream.eof() ) {
				throw MalformedInputError( "(unknown)",  number_of_snp_blocks ) ;
			}

			return number_of_snp_blocks ;
		}
	}
}
