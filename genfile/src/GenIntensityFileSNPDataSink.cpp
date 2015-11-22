
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <Eigen/Core>
#include <boost/bind.hpp>
#include <boost/format.hpp>
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
			IntensityWriter( Eigen::VectorXd& data ):
				m_data( data ),
				m_number_of_samples( m_data.size() / 2 ),
				m_number_of_alleles( 0 ),
				m_sample_i( 0 ),
				m_entry_i( 0 )
			{
				assert( m_data.size() % 2 == 0 ) ;
				m_data.setConstant( -1 ) ;
			}

			~IntensityWriter() throw() {}

			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				if( (nSamples*2) != m_data.size() ) {
					throw genfile::BadArgumentError(
						"genfile::IntensityWriter::initialise()",
						"n=" + string_utils::to_string(nSamples),
						"Number of samples does not match expected number (" + string_utils::to_string( m_data.size() / 2 ) + ")"
					) ;
				}
				if( nAlleles != 2 ) {
					throw genfile::BadArgumentError(
						"genfile::IntensityWriter::initialise()",
						( boost::format( "n=%d" ) % nAlleles).str(),
						"Expected two alleles."
					) ;
				}
				m_number_of_alleles = nAlleles ;
			}

			bool set_sample( std::size_t i ) {
				assert( i < m_number_of_samples ) ;
				m_sample_i = i ;
				return true ;
			}

			void set_number_of_entries( std::size_t n, OrderType const order_type, ValueType const value_type ) {
				if( n != m_number_of_alleles || order_type != ePerAllele ) {
					throw genfile::BadArgumentError(
						"genfile::IntensityWriter::set_number_of_entries()",
						"n=" + string_utils::to_string(n),
						"Expected one entry per allele per sample."
					) ;
				}
				m_entry_i = 0 ;
			}
			void set_value( MissingValue const value ) {
				m_data( m_sample_i * 2 + m_entry_i++ ) = -1 ;
			}

			void set_value( std::string& value ) {
				assert(0) ;
			}

			void set_value( Integer const value ) {
				assert(0) ;
			}
			
			void set_value( double const value ) {
				m_data( m_sample_i * 2 + m_entry_i++ ) = value ;
			}
			
			void write_to_stream( std::ostream& ostr ) const {
				for( std::size_t i = 0; i < m_data.size(); ++i ) {
					if( m_data(i) == -1 ) {
						ostr << " NA" ;
					} else {
						ostr << " " << m_data(i) ;
					}
				}
			}
			
		private:
			Eigen::VectorXd& m_data ;
			std::size_t const m_number_of_samples ;
			std::size_t m_number_of_alleles ;
			std::size_t m_sample_i ;
			std::size_t m_entry_i ;
		} ;
	}
	
	void GenIntensityFileSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) {
		if( write_chromosome_column() ) {
			stream() << "chromosome " ;
		}
		stream() << "SNPID rsid position alleleA alleleB" ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			stream() << " " << getter(i) << "_X" ;
			stream() << " " << getter(i) << "_Y" ;
		}
		stream() << "\n" ;
		m_data.resize( number_of_samples * 2 ) ;
		m_data.setConstant( -1 ) ;
	}

	void GenIntensityFileSNPDataSink::write_variant_data_impl(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		write_variant( stream(), id_data ) ;
		IntensityWriter writer( m_data ) ;
		data_reader.get( "XY", writer ) ;
		writer.write_to_stream( stream() ) ;
		stream() << "\n" ;
	}
}

