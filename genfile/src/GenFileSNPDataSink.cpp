
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
//#include "genfile/vcf/get_set.hpp"
//#include "genfile/vcf/get_set_eigen.hpp"
#include "genfile/ToGP.hpp"

namespace genfile {
	
	GenFileSNPDataSink::GenFileSNPDataSink( std::string const& filename, Chromosome chromosome ):
		GenLikeSNPDataSink( filename, chromosome, get_compression_type_indicated_by_filename( filename ) )
	{}

	GenFileSNPDataSink::GenFileSNPDataSink( std::string const& filename, Chromosome chromosome, CompressionType compression_type ):
		GenLikeSNPDataSink( filename, chromosome, compression_type )
	{}

	void GenFileSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) {
		m_data.resize( number_of_samples * 3 ) ;
		m_data.setConstant( -1 ) ;
	}

	namespace {
		struct GenotypeWriter: public VariantDataReader::PerSampleSetter {
			GenotypeWriter( Eigen::VectorXd& data ):
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
			void set_value( std::size_t, MissingValue const value ) {
				m_data( m_sample_i * 3 + m_entry_i++ ) = -1 ;
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
			std::size_t m_number_of_samples ;
			std::size_t m_sample_i ;
			std::size_t m_entry_i ;
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
			data_reader.get( ":genotypes:", to_GP( writer ) ) ;
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

