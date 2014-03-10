
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

// #define DEBUG_SQLITEHAPLOTYPESSNPDATASINK 1

SQLiteHaplotypesSNPDataSink::SQLiteHaplotypesSNPDataSink( qcdb::DBOutputter::UniquePtr outputter ):
	m_outputter( outputter ),
	m_number_of_samples( 0 ),
	m_data_i( 0 )
{
	m_outputter->connection().run_statement(
		"CREATE TABLE IF NOT EXISTS Haplotype ( "
		"  analysis_id INT NOT NULL REFERENCES Entity( id ), "
		"  variant_id INT NOT NULL REFERENCES Variant( id ), "
		"  N INT NOT NULL,"
		"  data BLOB NOT NULL"
		")"
	) ;
	m_outputter->connection().run_statement(
		"CREATE View IF NOT EXISTS HaplotypeView AS "
		"SELECT analysis_id, E.name AS analysis, variant_id, chromosome, position, rsid, alleleA, alleleB, quote( data ) "
		"FROM Haplotype H "
		"INNER JOIN Entity E "
		"ON E.id == H.analysis_id "
		"INNER JOIN Variant V "
		"ON V.id == H.variant_id"
	) ;
	m_insert_data_stmnt = m_outputter->connection().get_statement(
		"INSERT INTO Haplotype ( analysis_id, variant_id, N, data ) VALUES( ?, ?, ?, ? ) ;"
	) ;
	m_outputter->connection().run_statement(
		"CREATE TABLE IF NOT EXISTS Sample ( "
		"  analysis_id INTEGER NOT NULL REFERENCES Entity( id ), "
		"  identifier TEXT NOT NULL, "
		"  index_in_data INTEGER NOT NULL, "
		"  UNIQUE( analysis_id, identifier ), "
		"  UNIQUE( analysis_id, index_in_data )"
		")"
	) ;
	m_insert_sample_stmnt = m_outputter->connection().get_statement(
		"INSERT INTO Sample ( analysis_id, identifier, index_in_data ) VALUES( ?, ?, ? ) ;"
	) ;
	
	m_snps.resize( 10000 ) ;
	m_data.resize( 10000 ) ;
}

std::string SQLiteHaplotypesSNPDataSink::get_spec() const {
	return "SQLiteHaplotypesSNPDataSink" ;
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
		
		HaplotypeWriter( std::vector< char >* compression_buffer ):
			m_compression_buffer( compression_buffer ),
			m_buffer_validity( 2, true ),
			m_ploidy( 0 )
		{
			assert( compression_buffer != 0 ) ;
		}

		~HaplotypeWriter() throw() {}
		
		void set_number_of_samples( std::size_t n ) {
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
		void set_sample( std::size_t i ) {
			if( i > 0 ) {
				*(m_buf_p[ eUBJSON ]++) = ']' ;
			}
			*(m_buf_p[ eUBJSON ]++) = '[' ;
		}
		void set_number_of_entries( std::size_t n ) {
			if( n != 1 && n != 2 ) {
				m_buffer_validity[ eBITPACK ] = false ;
				throw genfile::BadArgumentError(
					"genfile::HaplotypeWriter::set_number_of_entries()",
					"n != 1 or 2"
				) ;
			}
			m_ploidy = n ;
		}

		void set_order_type( OrderType order_type ) {
			if( order_type != eOrderedList ) {
				throw genfile::BadArgumentError(
					"genfile::HaplotypeWriter::set_order_type()",
					"order_type=eUnorderedList"
				) ;
			}
		}

		void operator()( genfile::MissingValue const value ) {
			assert( ( m_buf_p[ eUBJSON ] + 1 ) <= m_buf_end_p[ eUBJSON ] ) ;
			*(m_buf_p[ eUBJSON ]++) = 'Z' ;
			if( m_buffer_validity[ eBITPACK ] ) {
				*(m_buf_p[ eBITPACK ]++) = static_cast< char >( -1 ) ; // missing value encoded by -1.
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
			if( fits_in< Integer, char >( value ) ) {
				assert( ( m_buf_p[ eUBJSON ] + 2 ) <= m_buf_end_p[ eUBJSON ] ) ;
				*(m_buf_p[ eUBJSON ]++) = 'i' ;
				m_buf_p[ eUBJSON ] = genfile::write_big_endian_integer( m_buf_p[ eUBJSON ], m_buf_end_p[ eUBJSON ], static_cast< char >( value )) ;

				if( m_buffer_validity[ eBITPACK ] ) {
					*(m_buf_p[ eBITPACK ]++) = static_cast< char >( value ) ;
					if( m_ploidy == 1 ) {
						// store two haplotypes for males, the second completely filled with missing values.
						*(m_buf_p[ eBITPACK ]++) = static_cast< char >( -1 ) ;
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
		std::vector< char >* m_compression_buffer ;
		std::vector< std::vector< char > > m_buffers ;
		std::size_t m_bitpack_index ;
		std::vector< bool > m_buffer_validity ;
		std::vector< char* > m_buf_p ;
		std::vector< char* > m_buf_end_p ;
		std::size_t m_ploidy ;

	private:
		char* write_float( char* buf_p, char* buf_end_p, float const& value ) const {
			char const* v_p = reinterpret_cast< char const* >( &value ) ;
			*(buf_p++) = *v_p++ ;
			*(buf_p++) = *v_p++ ;
			*(buf_p++) = *v_p++ ;
			*(buf_p++) = *v_p++ ;
			return buf_p ;
		}
	} ;
}

void SQLiteHaplotypesSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) {
	db::Connection::ScopedTransactionPtr transaction = m_outputter->connection().open_transaction( 600 ) ;
	for( std::size_t i = 0; i < number_of_samples; ++i ) {
		m_insert_sample_stmnt
			->bind( 1, m_outputter->analysis_id() )
			.bind( 2, getter(i) )
			.bind( 3, int64_t( i ) )
			.step()
		;
		m_insert_sample_stmnt->reset() ;
	}
	m_number_of_samples = number_of_samples ;
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
	data_reader.get( "genotypes", writer ) ;
	writer.finalise() ;
	++m_data_i ;
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
		char const* buffer = &(m_data[i][0]) ;
		char const* end_buffer = &(m_data[i][0]) + m_data[i].size() ;

		m_insert_data_stmnt
			->bind( 1, m_outputter->analysis_id() )
			.bind( 2, variant_ids[i] )
			.bind( 3, int64_t( m_number_of_samples ) )
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


