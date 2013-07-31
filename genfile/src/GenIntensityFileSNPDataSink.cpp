
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
#include "genfile/GenIntensityFileSNPDataSink.hpp"
#include "genfile/vcf/get_set.hpp"
#include "genfile/vcf/get_set_eigen.hpp"

namespace genfile {
	GenIntensityFileSNPDataSink::GenIntensityFileSNPDataSink( std::string const& filename, Chromosome chromosome ):
		GenLikeSNPDataSink( filename, chromosome, get_compression_type_indicated_by_filename( filename ) )
	{}

	GenIntensityFileSNPDataSink::GenIntensityFileSNPDataSink( std::string const& filename, Chromosome chromosome, CompressionType compression_type ):
		GenLikeSNPDataSink( filename, chromosome, compression_type )
	{}

	namespace {
		struct IntensityWriter: public VariantDataReader::PerSampleSetter {
			IntensityWriter( std::ostream& stream, std::size_t number_of_samples ):
				m_stream( stream ),
				m_number_of_samples( number_of_samples )
			{}

			~IntensityWriter() throw() {}

			void set_number_of_samples( std::size_t n ) {
				if( n != m_number_of_samples ) {
					throw genfile::BadArgumentError( "genfile::IntensityWriter::set_number_of_samples()", "n=" + string_utils::to_string(n), "Number of samples does not match expected number (" + string_utils::to_string( m_number_of_samples ) + ")" ) ;
				}
			}

			void set_sample( std::size_t i ) {
				assert( i < m_number_of_samples ) ;
			}
			void set_number_of_entries( std::size_t n ) {
				if( n != 2 ) {
					throw genfile::BadArgumentError( "genfile::IntensityWriter::set_number_of_entries()", "n=" + string_utils::to_string(n), "Expected 2 entries per sample." ) ;
				}
			}
			void operator()( MissingValue const value ) {
				m_stream << " NA" ;
			}

			void operator()( std::string& value ) {
				m_stream << " " << value ;
			}

			void operator()( Integer const value ) {
				m_stream << " " << value ;
			}
			
			void operator()( double const value ) {
				m_stream << " " << value ;
			}
			
		private:
			std::ostream& m_stream ;
			std::size_t const m_number_of_samples ;
		} ;
	}
	
	void GenIntensityFileSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) {
		m_number_of_samples = number_of_samples ;
		if( write_chromosome_column() ) {
			stream() << "chromosome " ;
		}
		stream() << "SNPID rsid position alleleA alleleB" ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			stream() << " " << getter(i) << "_X" ;
			stream() << " " << getter(i) << "_Y" ;
		}
		stream() << "\n" ;
	}

	void GenIntensityFileSNPDataSink::write_variant_data_impl(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		write_variant( stream(), id_data ) ;
		IntensityWriter writer( stream(), m_number_of_samples ) ;
		data_reader.get( "intensities", writer ) ;
		stream() << "\n" ;
	}
}

