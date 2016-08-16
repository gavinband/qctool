
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
#include "genfile/GenFileSNPDataSink.hpp"
//#include "genfile/vcf/get_set.hpp"
//#include "genfile/vcf/get_set_eigen.hpp"
#include "genfile/ToGP.hpp"

namespace genfile {
	
	GenFileSNPDataSink::GenFileSNPDataSink( std::string const& filename ):
		GenLikeSNPDataSink( filename, get_compression_type_indicated_by_filename( filename ) )
	{}

	GenFileSNPDataSink::GenFileSNPDataSink( std::string const& filename, CompressionType compression_type ):
		GenLikeSNPDataSink( filename, compression_type )
	{}

	void GenFileSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) {
		m_data.resize( number_of_samples * 3 ) ;
		m_data.setConstant( -1 ) ;
	}

	void GenFileSNPDataSink::set_precision( std::size_t precision ) {
		assert( precision == 5 ) ;
	}

	namespace {
		template< typename Num >
		struct float_policy: public boost::spirit::karma::real_policies< Num >
		{
			static bool trailing_zeros( Num n ) { return true ; }
			static int floatfield( Num n ) { return boost::spirit::karma::real_policies< Num >::fmtflags::fixed ; }
			static unsigned precision( Num n ) { return 5 ; }
		} ;

		typedef boost::spirit::karma::real_generator<double, float_policy<double> > FormatFloat ;

		struct GenotypeWriter: public VariantDataReader::PerSampleSetter {
			GenotypeWriter( Eigen::VectorXd& data, std::size_t const precision = 5 ):
				m_data( data ),
				m_number_of_samples( m_data.size() / 3 ),
				m_sample_i(0),
				m_entry_i(0)
			{
				m_data.setConstant( -1 ) ;
			}

			~GenotypeWriter() throw() {}
			
			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				if( nSamples != m_number_of_samples ) {
					throw genfile::BadArgumentError(
						"genfile::GenotypeWriter::initialise()",
						"n=" + string_utils::to_string( nSamples ),
						"Number of samples does not match expected number (" + string_utils::to_string( m_number_of_samples ) + ")"
					) ;
				}
				if( nAlleles != 2 ) {
					throw genfile::BadArgumentError(
						"genfile::GenotypeWriter::initialise()",
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
				if( n != 3 || order_type != ePerUnorderedGenotype || value_type != eProbability ) {
					throw genfile::BadArgumentError(
						"genfile::IntensityWriter::set_number_of_entries()",
						"n=" + string_utils::to_string(n),
						"Expected 3 genotype probabilities per sample."
					) ;
				}
				m_entry_i = 0 ;
			}
			void set_value( std::size_t j, MissingValue const value ) {
				// Value of -1 means missing data for this sample
				m_data[ 3 * m_sample_i + 0 ] = -1 ;
			}
			void set_value( std::size_t, std::string& value ) {
				assert(0) ;
			}
			void set_value( std::size_t, Integer const value ) {
				m_data( m_sample_i * 3 + m_entry_i++ ) = value ;
			}
			void set_value( std::size_t, double const value ) {
				m_data( m_sample_i * 3 + m_entry_i++ ) = value ;
			}
			void finalise() {} ;

			void write_to_stream( std::ostream& stream ) {
				for( std::size_t i = 0; i < m_number_of_samples; ++i ) {
					if( m_data[3*i+0] < 0 ) {
						stream << " 0 0 0 " ;
					} else {
						using boost::spirit::karma::double_ ;
						using boost::spirit::karma::generate ;
						FormatFloat formatter ;
						//sprintf( &m_buffer[0], " %.5f %.5f %.5f", m_data[3*i+0], m_data[3*i+1], m_data[3*i+2] ) ;
						generate( &m_buffer[0], ( " " << formatter << " " << formatter << " " << formatter << "\0" ), m_data[3*i+0], m_data[3*i+1], m_data[3*i+2] ) ;
						stream << m_buffer ;
						//generate( std::ostream_iterator<char>( stream ), ( " " << formatter << " " << formatter << " " << formatter ), m_data[3*i+0], m_data[3*i+1], m_data[3*i+2] ) ;
						/*
						stream 
							<< " " << m_data[3*i+0] 
							<< " " << m_data[3*i+1] 
							<< " " << m_data[3*i+2]
						;
						*/
					}
				}
			}

		private:
			Eigen::VectorXd& m_data ;
			std::size_t m_number_of_samples ;
			std::size_t m_sample_i ;
			std::size_t m_entry_i ;
			char m_buffer[100];
		} ;
	}

	void GenFileSNPDataSink::write_variant_data_impl(
		VariantIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		write_variant( stream(), id_data ) ;
		GenotypeWriter writer( m_data ) ;
		if( data_reader.supports( ":genotypes:" )) {
			data_reader.get( ":genotypes:", to_GP_unphased( writer ) ) ;
		} else {
			throw genfile::BadArgumentError(
				"genfile::GenFileSNPDataSink::write_variant_data_impl()",
				"data_reader",
				"Data source must support :genotypes: field."
			) ;
		}
		writer.write_to_stream( stream() ) ;
		stream() << "\n" ;
	}
}

