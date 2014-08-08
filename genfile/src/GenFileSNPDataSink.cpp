
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/GenLikeSNPDataSink.hpp"
#include "genfile/GenFileSNPDataSink.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/vcf/get_set_eigen.hpp"

namespace genfile {
	
	GenFileSNPDataSink::GenFileSNPDataSink( std::string const& filename, Chromosome chromosome ):
		GenLikeSNPDataSink( filename, chromosome, get_compression_type_indicated_by_filename( filename ) )
	{}

	GenFileSNPDataSink::GenFileSNPDataSink( std::string const& filename, Chromosome chromosome, CompressionType compression_type ):
		GenLikeSNPDataSink( filename, chromosome, compression_type )
	{}

	namespace {
		struct GenotypeWriter: public vcf::GenotypeSetterBase {
			GenotypeWriter( std::ostream& stream ):
				m_stream( stream )
			{}
			~GenotypeWriter() throw() {}
			
			void set( std::size_t i, double AA, double AB, double BB ) {
				m_stream << " " << AA << " " << AB << " " << BB ;
			}

		private:
			std::ostream& m_stream ;
		} ;
	}

	void GenFileSNPDataSink::write_variant_data_impl(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		write_variant( stream(), id_data ) ;
		GenotypeWriter writer( stream() ) ;
		if( data_reader.supports( "GT" )) {
			data_reader.get( "GT", writer ) ;
		} else if( data_reader.supports( "GP" )) {
			data_reader.get( "GP", writer ) ;
		} else {
			throw genfile::BadArgumentError(
				"genfile::GenFileSNPDataSink::write_variant_data_impl()",
				"data_reader",
				"Data source must support GT or GP fields."
			) ;
		}
		stream() << "\n" ;
	}
}

