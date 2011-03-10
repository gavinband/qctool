#ifndef GENFILE_VCF_GET_SET_HPP
#define GENFILE_VCF_GET_SET_HPP

#include "genfile/VariantEntry.hpp"

namespace genfile {
	namespace vcf {
		template< typename Setter >
		struct GenotypeProbabilitySetter
		{
			GenotypeProbabilitySetter( Setter setter ): m_setter( setter ) {}
			void operator()( std::size_t i, std::vector< vcf::Entry > const& values ) {
				assert( values.size() >= 3 ) ;
				double AA = values[0].is_missing() ? 0.0 : values[0].as< double >() ;
				double AB = values[1].is_missing() ? 0.0 : values[1].as< double >() ;
				double BB = values[2].is_missing() ? 0.0 : values[2].as< double >() ;
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
