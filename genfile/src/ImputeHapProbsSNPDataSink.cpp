
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <boost/spirit/include/karma.hpp>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/GenLikeSNPDataSink.hpp"
#include "genfile/ImputeHapProbsSNPDataSink.hpp"
//#include "genfile/vcf/get_set.hpp"
//#include "genfile/vcf/get_set_eigen.hpp"
#include "genfile/ToGP.hpp"
#include "genfile/format_float.hpp"

namespace genfile {
	
	ImputeHapProbsSNPDataSink::ImputeHapProbsSNPDataSink( std::string const& filename ):
		GenLikeSNPDataSink( filename, get_compression_type_indicated_by_filename( filename ) ),
		m_precision( 5 )
	{}

	ImputeHapProbsSNPDataSink::ImputeHapProbsSNPDataSink( std::string const& filename, CompressionType compression_type ):
		GenLikeSNPDataSink( filename, compression_type ),
		m_precision( 5 )
	{}

	void ImputeHapProbsSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) {
		m_data.resize( number_of_samples * 4 ) ;
		m_data.setConstant( -1 ) ;
	}

	void ImputeHapProbsSNPDataSink::set_precision( std::size_t precision ) {
        m_precision = precision ;
	}

	namespace {
		struct HapProbWriter: public VariantDataReader::PerSampleSetter {
			HapProbWriter( Eigen::VectorXd& data, std::size_t const precision = 5 ):
				m_data( data ),
				m_precision( precision ),
				m_number_of_samples( m_data.size() / 4 ),
				m_sample_i(0),
				m_entry_i(0)
			{
				m_data.setConstant( -1 ) ;
			}

			~HapProbWriter() throw() {}
			
			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				if( nSamples != m_number_of_samples ) {
					throw genfile::BadArgumentError(
						"genfile::HapProbWriter::initialise()",
						"n=" + string_utils::to_string( nSamples ),
						"Number of samples does not match expected number (" + string_utils::to_string( m_number_of_samples ) + ")"
					) ;
				}
				if( nAlleles != 2 ) {
					throw genfile::BadArgumentError(
						"genfile::HapProbWriter::initialise()",
						( boost::format( "n=%d" ) % nAlleles ).str(),
						"Expected two alleles."
					) ;
				}
			}

			bool set_sample( std::size_t i ) {
				assert( i < m_number_of_samples ) ;
				m_sample_i = i ;
				return true ;
			}

			void set_number_of_entries( uint32_t, std::size_t n, OrderType const order_type, ValueType const value_type ) {
				if( n != 4 || order_type != ePerPhasedHaplotypePerAllele || value_type != eProbability ) {
					throw genfile::BadArgumentError(
						"genfile::IntensityWriter::set_number_of_entries()",
						"n=" + string_utils::to_string(n),
						"Expected 4 genotype probabilities per sample."
					) ;
				}
				m_entry_i = 0 ;
			}
			void set_value( std::size_t j, MissingValue const value ) {
				// Value of -1 means missing data for this sample
				m_data[ 4 * m_sample_i + 0 ] = -1 ;
			}
			void set_value( std::size_t, std::string& value ) {
				assert(0) ;
			}
			void set_value( std::size_t, Integer const value ) {
				m_data( m_sample_i * 4 + m_entry_i++ ) = value ;
			}
			void set_value( std::size_t, double const value ) {
				m_data( m_sample_i * 4 + m_entry_i++ ) = value ;
			}
			void finalise() {} ;

			void write_to_stream( std::ostream& stream ) {
				for( std::size_t i = 0; i < m_number_of_samples; ++i ) {
					if( m_data[2*i+0] < 0 ) {
						stream << " 0 0" ;
					} else {
						if( m_precision == 5 ) {
							// Fast formatting at 5dps.
							impl::FormatFloat5dp formatter ;
							boost::spirit::karma::generate(
								&m_buffer[0],
								( " " << formatter << " " << formatter << "\0" ),
								m_data[4*i+1],
								m_data[4*i+3]
							) ;
							stream << m_buffer ;
						} else {
							// Slow but flexible formatting at any precision.
							stream 
								<< std::setprecision( m_precision )
								<< " " << m_data[4*i+1] 
								<< " " << m_data[4*i+3] 
							;
						}
					}
				}
			}

		private:
			Eigen::VectorXd& m_data ;
			std::size_t const m_precision ;
			std::size_t m_number_of_samples ;
			std::size_t m_sample_i ;
			std::size_t m_entry_i ;
			char m_buffer[100];
		} ;
	}

	void ImputeHapProbsSNPDataSink::write_variant_data_impl(
		VariantIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		write_variant( stream(), id_data ) ;
		HapProbWriter writer( m_data, m_precision ) ;
		if( data_reader.supports( ":genotypes:" )) {
			data_reader.get( ":genotypes:", to_GP_phased( writer ) ) ;
		} else {
			throw genfile::BadArgumentError(
				"genfile::ImputeHapProbsSNPDataSink::write_variant_data_impl()",
				"data_reader",
				"Data source must support :genotypes: field."
			) ;
		}
		writer.write_to_stream( stream() ) ;
		stream() << "\n" ;
	}
}

