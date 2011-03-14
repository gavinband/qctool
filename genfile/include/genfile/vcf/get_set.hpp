#ifndef GENFILE_VCF_GET_SET_HPP
#define GENFILE_VCF_GET_SET_HPP

#include "genfile/VariantEntry.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	namespace vcf {
		template< typename Setter >
		struct GenotypeProbabilitySetter
		{
			GenotypeProbabilitySetter( Setter setter ): m_setter( setter ) {}
			void operator()( std::size_t i, std::vector< vcf::Entry > const& values ) {
				double AA = 0.0, AB = 0.0, BB = 0.0 ; // zero genotypes to represent missing call.
				if(
					values[0].is_int()
				) {
					if(
						values.size() == 2
					) {
						int A = values[0].as< int >(),
							B = values[1].as< int >() ;

						if( A >= 0 && A < 2 && B >= 0 && B < 2 ) {
							AA = ( A == 0 && B == 0 ) ;
							AB = ( A + B == 1 ) ;
							BB = ( A == 1 && B == 1 ) ;
						}
					}
				}
				else if( values[0].is_double() && values.size() == 3 ) {
					AA = values[0].is_missing() ? 0.0 : values[0].as< double >() ;
					AB = values[1].is_missing() ? 0.0 : values[1].as< double >() ;
					BB = values[2].is_missing() ? 0.0 : values[2].as< double >() ;
				}
				m_setter( i, AA, AB, BB ) ;
			}
		private:
			Setter m_setter ;
		} ;
		
		template< typename Setter >
		GenotypeProbabilitySetter< Setter > make_genotype_probability_setter( Setter setter ) {
			return GenotypeProbabilitySetter< Setter >( setter ) ;
		}
	}
}

#endif
