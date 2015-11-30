
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/ShapeITHaplotypesSNPDataSink.hpp"

namespace genfile {
	
	ShapeITHaplotypesSNPDataSink::ShapeITHaplotypesSNPDataSink( std::string const& filename, Chromosome chromosome ):
		GenLikeSNPDataSink( filename, chromosome, get_compression_type_indicated_by_filename( filename ) )
	{}

	ShapeITHaplotypesSNPDataSink::ShapeITHaplotypesSNPDataSink( std::string const& filename, Chromosome chromosome, CompressionType compression_type ):
		GenLikeSNPDataSink( filename, chromosome, compression_type )
	{}

	void ShapeITHaplotypesSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter ) {
		m_data.resize( number_of_samples * 2 ) ;
		m_data.setConstant( -1 ) ;
	}

	namespace {
		struct HaplotypeWriter: public VariantDataReader::PerSampleSetter {
			HaplotypeWriter( Eigen::VectorXd& data ):
				m_data( data ),
				m_sample_i(0),
				m_entry_i(0)
			{
				m_data.setConstant( -1 ) ;
			}

			~HaplotypeWriter() throw() {}
			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				assert( (2*nSamples) == m_data.size() ) ;
				assert( nAlleles == 2 ) ;
			}
			void set_number_of_alleles( std::size_t n ) {
				assert( n == 2 ) ;
			}
			bool set_sample( std::size_t i ) {
				m_sample_i = i ;
				return true ;
			}
			void set_number_of_entries( uint32_t, std::size_t n, OrderType const order_type, ValueType const value_type ) {
				if( n != 2 || order_type != ePerOrderedHaplotype || value_type != eAlleleIndex) {
					throw genfile::BadArgumentError(
						"genfile::HaplotypeWriter::set_number_of_entries()",
						"n != 2",
						"Expected two entries representing ordered haplotype calls"
					) ;
				}
				m_entry_i = 0 ;
			}
			virtual void set_order_type( OrderType order_type ) {
				if( order_type != eOrderedList ) {
					throw genfile::BadArgumentError(
						"genfile::HaplotypeWriter::set_order_type()",
						"order_type=eUnorderedList"
					) ;
				}
			}
			void set_value( MissingValue const value ) {
				m_data( m_sample_i * 2 + m_entry_i++ ) = -1 ;
			}
			void set_value( std::string& value ) {
				assert(0) ;
			}
			void set_value( Integer const value ) {
				m_data( m_sample_i * 2 + m_entry_i++ ) = value ;
			}
			void set_value( double const value ) {
				m_data( m_sample_i * 2 + m_entry_i++ ) = value ;
			}
			
			void write_to_stream( std::ostream& ostr ) const {
				for( std::size_t i = 0; i < m_data.size(); ++i ) {
					if( m_data(i) == -1 ) {
						ostr << " NA" ;
					} else {
						ostr << " " << m_data[i] ;
					}
				}
			}
			
		private:
			Eigen::VectorXd& m_data ;
			std::size_t m_sample_i ;
			std::size_t m_entry_i ;
		} ;
	}

	void ShapeITHaplotypesSNPDataSink::write_variant_data_impl(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		write_variant( stream(), id_data ) ;
		HaplotypeWriter writer( m_data ) ;
		data_reader.get( "genotypes", writer ) ;
		writer.write_to_stream( stream() ) ;
		stream() << "\n" ;
	}
}

