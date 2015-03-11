
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

	namespace {
		struct HaplotypeWriter: public VariantDataReader::PerSampleSetter {
			HaplotypeWriter( std::ostream& stream ):
				m_stream( stream )
			{}
			~HaplotypeWriter() throw() {}
			
			void set_number_of_samples( std::size_t n ) {
				// nothing to do.
			}
			void set_number_of_alleles( std::size_t n ) {
				assert( n == 2 ) ;
			}
			bool set_sample( std::size_t i ) {
				// nothing to do
				return true ;
			}
			void set_order_type( OrderType const order_type, ValueType const value_type ) {
				assert( order_type == ePerOrderedHaplotype && value_type == eAlleleIndex ) ;
			}

			virtual void set_number_of_entries( std::size_t n ) {
				if( n != 2 ) {
					throw genfile::BadArgumentError(
						"genfile::HaplotypeWriter::set_number_of_entries()",
						"n != 2"
					) ;
				}
			}
			virtual void set_order_type( OrderType order_type ) {
				if( order_type != eOrderedList ) {
					throw genfile::BadArgumentError(
						"genfile::HaplotypeWriter::set_order_type()",
						"order_type=eUnorderedList"
					) ;
				}
			}
			virtual void operator()( MissingValue const value ) {
				m_stream << " " << value ;
			}
			virtual void operator()( std::string& value ) {
				m_stream << " " << value ;
			}
			virtual void operator()( Integer const value ) {
				m_stream << " " << value ;
			}
			virtual void operator()( double const value ) {
				m_stream << " " << value ;
			}
			
		private:
			std::ostream& m_stream ;
		} ;
	}

	void ShapeITHaplotypesSNPDataSink::write_variant_data_impl(
		SNPIdentifyingData const& id_data,
		VariantDataReader& data_reader,
		Info const& info
	) {
		write_variant( stream(), id_data ) ;
		HaplotypeWriter writer( stream() ) ;
		data_reader.get( "genotypes", writer ) ;
		stream() << "\n" ;
	}
}

