
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "stdint.h"
#include <limits>
#include <boost/algorithm/string/replace.hpp>
#include <Eigen/Core>
#include "genfile/snp_data_utils.hpp"
#include "genfile/endianness_utils.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/zlib.hpp"
#include "qcdb/DBOutputter.hpp"
#include "components/SNPOutputComponent/SQLiteHaplotypesSNPDataSink.hpp"
#include "qcdb/store_samples_in_db.hpp"

// #define DEBUG_SQLITEHAPLOTYPESSNPDATASINK 1

SQLiteHaplotypesSNPDataSink::SQLiteHaplotypesSNPDataSink(
	qcdb::DBOutputter::UniquePtr outputter,
	genfile::CohortIndividualSource const& samples
):
	m_genotype_field( ":genotypes:" ),
	m_outputter( outputter ),
	m_samples( samples ),
	m_data_i( 0 )
{
	m_outputter->connection().run_statement(
		"CREATE TABLE IF NOT EXISTS Haplotype ( "
		"  analysis_id INT NOT NULL REFERENCES Analysis( id ), "
		"  variant_id INT NOT NULL REFERENCES Variant( id ), "
		"  N INT NOT NULL,"
		"  data BLOB NOT NULL"
		")"
	) ;
	m_outputter->connection().run_statement(
		"CREATE View IF NOT EXISTS HaplotypeView AS "
		"SELECT analysis_id, A.name AS analysis, variant_id, chromosome, position, rsid, alleleA, alleleB, quote( data ) "
		"FROM Haplotype H "
		"INNER JOIN Analysis A "
		"ON A.id == H.analysis_id "
		"INNER JOIN Variant V "
		"ON V.id == H.variant_id"
	) ;
	m_insert_data_stmnt = m_outputter->connection().get_statement(
		"INSERT INTO Haplotype ( analysis_id, variant_id, N, data ) VALUES( ?, ?, ?, ? ) ;"
	) ;
	m_snps.resize( 10000 ) ;
	m_data.resize( 10000 ) ;
}

std::string SQLiteHaplotypesSNPDataSink::get_spec() const {
	return "SQLiteHaplotypesSNPDataSink" ;
}

void SQLiteHaplotypesSNPDataSink::set_genotype_field( std::string field ) {
	m_genotype_field = field ;
}

namespace {
	template< typename Integer, typename Target >
	bool fits_in( Integer const& value ) {
		// This code adapted from code found on StackOverflow
		// here: http://stackoverflow.com/a/17251989/517251.

		// Return false if Integer can hold lower numbers than target and the value is less than the minimum of target
		// Or if Integer can hold larger numbers than Target and the value is larger than the target.
		// Otherwise return true.
		//
		// (I guess the point here is that the first part of each test can be done at compile time.)
		//
		const intmax_t IntegerMin = intmax_t( std::numeric_limits< Integer >::min() ) ;
		const uintmax_t IntegerMax = uintmax_t( std::numeric_limits< Integer >::max() ) ;
		const intmax_t TargetMin = intmax_t( std::numeric_limits< Target >::min() ) ;
		const uintmax_t TargetMax = uintmax_t( std::numeric_limits< Target >::max() ) ;
       

		bool result = !(
			( ( IntegerMin < TargetMin ) && value < static_cast< Integer >( TargetMin ))
			||
			( ( IntegerMax > TargetMax ) && value > static_cast< Integer >( TargetMax ))
		) ;
		
		return result ;
	}

	struct HaplotypeWriter: public genfile::VariantDataReader::PerSampleSetter {
		enum Encoding { eUBJSON = 0, eBITPACK = 1 } ; // must take contiguous values starting at zero.
		
		HaplotypeWriter( std::vector< uint8_t >* compression_buffer ):
			m_compression_buffer( compression_buffer ),
			m_buffer_validity( 2, true ),
			m_ploidy( 0 )
		{
			assert( compression_buffer != 0 ) ;
		}

		~HaplotypeWriter() throw() {}
		
		void initialise( std::size_t n, std::size_t nAlleles ) {
			assert( nAlleles == 2 ) ;
			m_buffers.resize(2) ;
			m_buf_p.resize(2) ;
			m_buf_end_p.resize(2) ;

			{
				m_buffers[ eUBJSON ].clear() ;
				m_buffers[ eUBJSON ].resize( 9 + ( n * 2 * 10 ) + 2 ) ; // name, plus two brackets per sample, two values up to 9 bytes per sample, two outside brackets.
				m_buf_p[ eUBJSON ] = &(m_buffers[ eUBJSON ][0]) ;
				m_buf_end_p[ eUBJSON ] = &(m_buffers[ eUBJSON ][0]) + m_buffers[ eUBJSON ].size() ;
				std::string const encoding_name = "ubjson" ;
				m_buffers[ eUBJSON ][0] = 's' ;
				++m_buf_p[ eUBJSON ] ;
				m_buf_p[ eUBJSON ] = genfile::write_big_endian_integer( m_buf_p[ eUBJSON ], m_buf_end_p[ eUBJSON ], static_cast< int16_t >( encoding_name.size() )) ;
				m_buf_p[ eUBJSON ] = std::copy( encoding_name.begin(), encoding_name.end(), m_buf_p[ eUBJSON ] ) ;
				// Write the UBJSON array marker
				*(m_buf_p[ eUBJSON ]++) = '[' ;
			}

			{
				m_buffers[ eBITPACK ].clear() ;
				m_buffers[ eBITPACK ].resize( 10 + ( n * 2 ) ) ; // name, plus 2 haplotypes per sample.
				m_buf_p[ eBITPACK ] = &(m_buffers[ eBITPACK ][0]) ;
				m_buf_end_p[ eBITPACK ] = &(m_buffers[ eBITPACK ][0]) + m_buffers[ eBITPACK ].size() ;
				std::string const encoding_name = "bitpack" ;
				m_buffers[ eBITPACK ][0] = 's' ;
				++m_buf_p[ eBITPACK ] ;
				m_buf_p[ eBITPACK ] = genfile::write_big_endian_integer( m_buf_p[ eBITPACK ], m_buf_end_p[ eBITPACK ], static_cast< int16_t >( encoding_name.size() )) ;
				m_buf_p[ eBITPACK ] = std::copy( encoding_name.begin(), encoding_name.end(), m_buf_p[ eBITPACK ] ) ;
				m_bitpack_index = 0 ;
			}
		}
		
		bool set_sample( std::size_t i ) {
			if( i > 0 ) {
				*(m_buf_p[ eUBJSON ]++) = ']' ;
			}
			*(m_buf_p[ eUBJSON ]++) = '[' ;
			return true ;
		}
		void set_number_of_entries( std::size_t n, OrderType const order_type, ValueType const value_type ) {
			if( n != 1 && n != 2 ) {
				m_buffer_validity[ eBITPACK ] = false ;
				throw genfile::BadArgumentError(
					"genfile::HaplotypeWriter::set_number_of_entries()",
					"n != 1 or 2"
				) ;
			}
			if( order_type != genfile::ePerOrderedHaplotype ) {
				throw genfile::BadArgumentError(
					"genfile::HaplotypeWriter::set_order_type()",
					"order_type",
					"Expected values to represent ordered haplotype calls."
				) ;
			}
			if( value_type != genfile::eAlleleIndex ) {
				throw genfile::BadArgumentError(
					"genfile::HaplotypeWriter::set_order_type()",
					"value_type",
					"Expected values to represent allele indices (i.e. GT field)."
				) ;
			}
			m_ploidy = n ;
		}

		void operator()( genfile::MissingValue const value ) {
			assert( ( m_buf_p[ eUBJSON ] + 1 ) <= m_buf_end_p[ eUBJSON ] ) ;
			*(m_buf_p[ eUBJSON ]++) = 'Z' ;
			if( m_buffer_validity[ eBITPACK ] ) {
				*(m_buf_p[ eBITPACK ]++) = 0xFF ;
			}
		}

		void operator()( std::string& value ) {
			assert( ( m_buf_p[ eUBJSON ] + 1 + value.size() ) <= m_buf_end_p[ eUBJSON ] ) ;
			m_buffer_validity[ eBITPACK ] = false ;
			*(m_buf_p[ eUBJSON ]++) = 'S' ;
			std::copy( value.begin(), value.end(), m_buf_p[ eUBJSON ] ) ;
			m_buf_p[ eUBJSON ] += value.size() ;
		}

		void operator()( Integer const value ) {
			// We allow a maximum value of 127 here to allow for reserved values.
	 		if( fits_in< Integer, uint8_t >( value ) && value < 127 ) {
				assert( ( m_buf_p[ eUBJSON ] + 2 ) <= m_buf_end_p[ eUBJSON ] ) ;
				*(m_buf_p[ eUBJSON ]++) = 'i' ;
				m_buf_p[ eUBJSON ] = genfile::write_big_endian_integer( m_buf_p[ eUBJSON ], m_buf_end_p[ eUBJSON ], static_cast< uint8_t >( value )) ;

				if( m_buffer_validity[ eBITPACK ] ) {
					*(m_buf_p[ eBITPACK ]++) = static_cast< uint8_t >( value ) ;
					if( m_ploidy == 1 ) {
						// store two haplotypes for males, the second filled wth missing values.
						*(m_buf_p[ eBITPACK ]++) = 0xFF ;
					}
				}
			} else {
				m_buffer_validity[ eBITPACK ] = false ;
				if( fits_in< Integer, int16_t >( value ) ) {
					assert( ( m_buf_p[ eUBJSON ] + 3 ) <= m_buf_end_p[ eUBJSON ] ) ;
					*(m_buf_p[ eUBJSON ]++) = 'I' ;
					m_buf_p[ eUBJSON ] = genfile::write_big_endian_integer( m_buf_p[ eUBJSON ], m_buf_end_p[ eUBJSON ], static_cast< int16_t >( value )) ;
				} else if( fits_in< Integer, int32_t >( value ) ) {
					assert( ( m_buf_p[ eUBJSON ] + 5 ) <= m_buf_end_p[ eUBJSON ] ) ;
					*(m_buf_p[ eUBJSON ]++) = 'l' ;
					m_buf_p[ eUBJSON ] = genfile::write_big_endian_integer( m_buf_p[ eUBJSON ], m_buf_end_p[ eUBJSON ], static_cast< int32_t >( value )) ;
				}
				else {
					assert( ( m_buf_p[ eUBJSON ] + 9 ) <= m_buf_end_p[ eUBJSON ] ) ;
					*(m_buf_p[ eUBJSON ]++) = 'L' ;
					m_buf_p[ eUBJSON ] = genfile::write_big_endian_integer( m_buf_p[ eUBJSON ], m_buf_end_p[ eUBJSON ], static_cast< int64_t >( value )) ;
				}
			}
		}

		void operator()( double const value ) {
			assert( ( m_buf_p[ eUBJSON ] + 5 ) <= m_buf_end_p[ eUBJSON ] ) ;
			m_buffer_validity[ eBITPACK ] = false ;
			float const float_value = value ;
			*(m_buf_p[ eUBJSON ]++) = 'd' ;
			m_buf_p[ eUBJSON ] = write_float( m_buf_p[ eUBJSON ], m_buf_end_p[ eUBJSON ], float_value ) ;
		}
		
		void finalise() {
			*(m_buf_p[ eUBJSON ]++) = ']' ;
			*(m_buf_p[ eUBJSON ]++) = ']' ;

			for( std::size_t i = 0; i < 2; ++i ) {
				m_buffers[ i ].resize( ( m_buf_p[ i ] ) - ( &( m_buffers[ i ][0] ) ) ) ;
			}
			
			if( m_buffer_validity[ eBITPACK ] ) {
				genfile::zlib_compress( &( m_buffers[ eBITPACK ][0] ), m_buf_p[eBITPACK], m_compression_buffer ) ;
				//std::cerr << "BITPACK: buffer size: " << m_buffers[ eBITPACK ].size() << " before compression, " << m_compression_buffer->size() << " after compression.\n" ;
			} else {
				genfile::zlib_compress( &( m_buffers[ eUBJSON ][0] ), m_buf_p[eUBJSON], m_compression_buffer ) ;
				//std::cerr << "UBJSON: buffer size: " << m_buffers[ eUBJSON ].size() << " before compression, " << m_compression_buffer->size() << " after compression.\n" ;
			}
		}
		
		std::string get_encoding_name() const {
			if( m_buffer_validity[ eBITPACK ] ) {
				return "bitpack" ;
			} else {
				return "ubjson" ;
			}
		}
		
	private:
		std::vector< uint8_t >* m_compression_buffer ;
		std::vector< std::vector< uint8_t > > m_buffers ;
		std::size_t m_bitpack_index ;
		std::vector< bool > m_buffer_validity ;
		std::vector< uint8_t* > m_buf_p ;
		std::vector< uint8_t* > m_buf_end_p ;
		std::size_t m_ploidy ;

	private:
		uint8_t* write_float( uint8_t* buf_p, uint8_t const* const buf_end_p, float const& value ) const {
			uint8_t const* v_p = reinterpret_cast< uint8_t const* >( &value ) ;
			*(buf_p++) = *v_p++ ;
			*(buf_p++) = *v_p++ ;
			*(buf_p++) = *v_p++ ;
			*(buf_p++) = *v_p++ ;
			return buf_p ;
		}
	} ;
}

void SQLiteHaplotypesSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) {
	if( number_of_samples != m_samples.get_number_of_individuals() ) {
		throw genfile::BadArgumentError(
			"SQLiteGenotypesSNPDataSink::set_sample_names_impl()",
			"number_of_samples",
			"number of samples ("
			+ genfile::string_utils::to_string( number_of_samples )
			+ ") does not match that in sample files ("
			+ genfile::string_utils::to_string( m_samples.get_number_of_individuals() )
			+ ")."
		) ;
	}
	qcdb::store_samples_in_db(
		number_of_samples,
		getter,
		m_samples,
		*m_outputter,
		"Sample"
	) ;
}

void SQLiteHaplotypesSNPDataSink::write_variant_data_impl(
	genfile::SNPIdentifyingData const& id_data,
	genfile::VariantDataReader& data_reader,
	Info const& info
) {
	if( m_data_i == m_snps.size() ) {
		flush_data( m_data_i ) ;
		m_data_i = 0 ;
	}
	m_snps[ m_data_i ] = id_data ;
	HaplotypeWriter writer( &(m_data[ m_data_i ]) ) ;
	try {
		data_reader.get( ":genotypes:", writer ) ;
		writer.finalise() ;
		++m_data_i ;
	} catch( genfile::InputError const& e ) {
		std::cerr << "SQLiteHaplotypesSNPDataSink:write_variant_data_impl(): at variant " << id_data << ": "
			<< e.format_message() << ".\n" ;
		std::cerr << "I will not store haplotypes at this SNP.\n" ;
	}
}

void SQLiteHaplotypesSNPDataSink::finalise_impl() {
	if( m_data_i > 0 ) {
		flush_data( m_data_i ) ;
		m_data_i = 0 ;
	}
	m_outputter->finalise() ;
}

void SQLiteHaplotypesSNPDataSink::flush_data( std::size_t const data_count ) {
	db::Connection::ScopedTransactionPtr transaction = m_outputter->connection().open_transaction( 600 ) ;

#if DEBUG_SQLITEHAPLOTYPESSNPDATASINK
	std::cerr << "Flushing " << data_count << " elements..." ;
	std::size_t max_data_size = 0 ;
#endif
	std::vector< db::Connection::RowId > variant_ids( data_count ) ;
	for( std::size_t i = 0; i < data_count; ++i ) {
		variant_ids[i] = m_outputter->get_or_create_variant( m_snps[i] ) ;
	}
#if DEBUG_SQLITEHAPLOTYPESSNPDATASINK
	std::cerr << "stored variants..." ;
#endif
	for( std::size_t i = 0; i < data_count; ++i ) {
		uint8_t const* buffer = &(m_data[i][0]) ;
		uint8_t const* end_buffer = &(m_data[i][0]) + m_data[i].size() ;

		m_insert_data_stmnt
			->bind( 1, m_outputter->analysis_id() )
			.bind( 2, variant_ids[i] )
			.bind( 3, int64_t( m_samples.get_number_of_individuals() ) )
			.bind( 4, buffer, end_buffer )
			.step() ;
		m_insert_data_stmnt->reset() ;
#if DEBUG_SQLITEHAPLOTYPESSNPDATASINK
		max_data_size = std::max( max_data_size, m_data[i].size() ) ;
#endif
	}
#if DEBUG_SQLITEHAPLOTYPESSNPDATASINK
	std::cerr << "Done, max data size was " << max_data_size << ".\n" ;
#endif
}


