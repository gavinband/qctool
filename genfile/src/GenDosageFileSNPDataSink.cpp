
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
#include "genfile/GenDosageFileSNPDataSink.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/vcf/get_set_eigen.hpp"

namespace genfile {
	GenDosageFileSNPDataSink::GenDosageFileSNPDataSink( std::string const& filename, Chromosome chromosome ):
		GenLikeSNPDataSink( filename, chromosome, get_compression_type_indicated_by_filename( filename ) )
	{}

	GenDosageFileSNPDataSink::GenDosageFileSNPDataSink( std::string const& filename, Chromosome chromosome, CompressionType compression_type ):
		GenLikeSNPDataSink( filename, chromosome, compression_type )
	{}

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
	
	void GenDosageFileSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) {
		if( write_chromosome_column() ) {
			stream() << "chromosome " ;
		}
		stream() << "SNPID rsid position alleleA alleleB" ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			stream() << " " << getter(i) ;
		}
		stream() << "\n" ;
		m_data.resize( number_of_samples ) ;
		m_data.setZero() ;
	}

	void GenDosageFileSNPDataSink::write_variant_data_impl(
		SNPIdentifyingData const& id_data,
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

