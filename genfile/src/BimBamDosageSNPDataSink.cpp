
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
#include "genfile/BimBamDosageSNPDataSink.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/vcf/get_set_eigen.hpp"

namespace genfile {
	BimBamDosageSNPDataSink::BimBamDosageSNPDataSink( std::string const& filename ):
		GenLikeSNPDataSink( filename, get_compression_type_indicated_by_filename( filename ) )
	{}

	BimBamDosageSNPDataSink::BimBamDosageSNPDataSink( std::string const& filename, CompressionType compression_type ):
		GenLikeSNPDataSink( filename, compression_type )
	{}

	void BimBamDosageSNPDataSink::write_variant( std::ostream& out, genfile::VariantIdentifyingData const& variant ) {
		assert( variant.number_of_alleles() == 2 ) ;
		using genfile::string_utils::to_string ;
		out << variant.get_identifiers_as_string( "," ) << ":" << to_string( variant.get_position().chromosome() ) << ":" << to_string( variant.get_position().position() ) << " "
			<< variant.get_allele(0) << " "
			<< variant.get_allele(1) ;
	}

	namespace {
		struct GenotypeWriter: public vcf::GenotypeSetterBase {
			GenotypeWriter( Eigen::VectorXd& data ):
				m_data( data )
			{
				m_data.setConstant( -1 ) ;
			}
			~GenotypeWriter() throw() {}
			
			void set( std::size_t i, double AA, double AB, double BB ) {
				if( AA == 0 && AB == 0 && BB == 0 ) {
					m_data[i] = -1 ;
				} else {
					m_data[i] = ( ( 2 * BB ) + AB ) ;
				}
			}

			void write_to_stream( std::ostream& stream ) const {
				for( std::size_t i = 0; i < m_data.size(); ++i ) {
					if( m_data[i] == -1 ) {
						stream << " NA" ;
					} else {
						stream << " " << m_data[i] ;
					}
				}
			}

		private:
			Eigen::VectorXd& m_data ;
		} ;
	}
	
	void BimBamDosageSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) {
		m_data.resize( number_of_samples ) ;
		m_data.setZero() ;
	}

	void BimBamDosageSNPDataSink::write_variant_data_impl(
		VariantIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		write_variant( stream(), id_data ) ;
		GenotypeWriter writer( m_data ) ;
		data_reader.get( ":genotypes:", writer ) ;
		writer.write_to_stream( stream() ) ;
		stream() << "\n" ;
	}
}

