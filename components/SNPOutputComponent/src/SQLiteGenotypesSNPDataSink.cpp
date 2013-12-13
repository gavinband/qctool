
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
#include "components/SNPOutputComponent/SQLiteGenotypesSNPDataSink.hpp"

// #define DEBUG_SQLITEGENOTYPESNPDATASINK 1

SQLiteGenotypesSNPDataSink::SQLiteGenotypesSNPDataSink( qcdb::DBOutputter::UniquePtr outputter ):
	m_outputter( outputter ),
	m_number_of_samples( 0 ),
	m_data_i( 0 )
{
	m_outputter->connection().run_statement(
		"CREATE TABLE IF NOT EXISTS Genotype ( "
		"  analysis_id INT NOT NULL REFERENCES Entity( id ), "
		"  variant_id INT NOT NULL REFERENCES Variant( id ), "
		"  N INT NOT NULL,"
		"  data BLOB NOT NULL"
		")"
	) ;
	m_outputter->connection().run_statement(
		"CREATE View IF NOT EXISTS GenotypeView AS "
		"SELECT analysis_id, E.name AS analysis, variant_id, chromosome, position, rsid, alleleA, alleleB, quote( data ) "
		"FROM Genotype H "
		"INNER JOIN Entity E "
		"ON E.id == H.analysis_id "
		"INNER JOIN Variant V "
		"ON V.id == H.variant_id"
	) ;
	m_insert_data_stmnt = m_outputter->connection().get_statement(
		"INSERT INTO Genotype ( analysis_id, variant_id, N, data ) VALUES( ?, ?, ?, ? ) ;"
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

std::string SQLiteGenotypesSNPDataSink::get_spec() const {
	return "SQLiteGenotypesSNPDataSink" ;
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

	struct GenotypeProbabilityWriter: public genfile::VariantDataReader::PerSampleSetter {
		enum Encoding { eBITPACK = 0 } ; // must take contiguous values starting at zero.
		
		GenotypeProbabilityWriter( std::vector< char >* compression_buffer ):
			m_compression_buffer( compression_buffer ),
			m_buffer_validity( 1, true )
		{
			assert( compression_buffer != 0 ) ;
		}

		~GenotypeProbabilityWriter() throw() {}
		
		void set_number_of_samples( std::size_t n ) {
			m_buffers.resize(1) ;
			m_buf_p.resize(1) ;
			m_buf_end_p.resize(1) ;

			{
				m_buffers[ eBITPACK ].clear() ;
				std::string const encoding_name = "3float" ;
				m_buffers[ eBITPACK ].resize( 3 + encoding_name.size() + 5 + 3 + ( n * 4 * 3 ) ) ; // name, plus 3 genotype probabilities per sample.
				m_buf_p[ eBITPACK ] = &(m_buffers[ eBITPACK ][0]) ;
				m_buf_end_p[ eBITPACK ] = &(m_buffers[ eBITPACK ][0]) + m_buffers[ eBITPACK ].size() ;
				m_buffers[ eBITPACK ][0] = 's' ; // marker that a string follows, UBJSON style.
				++m_buf_p[ eBITPACK ] ;
				m_buf_p[ eBITPACK ] = genfile::write_big_endian_integer( m_buf_p[ eBITPACK ], m_buf_end_p[ eBITPACK ], static_cast< int16_t >( encoding_name.size() )) ;
				m_buf_p[ eBITPACK ] = std::copy( encoding_name.begin(), encoding_name.end(), m_buf_p[ eBITPACK ] ) ;

				*(m_buf_p[ eBITPACK ]++) = 'I' ; // marker that an int32 follows, UBJSON style.
				m_buf_p[ eBITPACK ] = genfile::write_big_endian_integer( m_buf_p[ eBITPACK ], m_buf_end_p[ eBITPACK ], static_cast< int32_t >( n )) ; // number of samples

				*(m_buf_p[ eBITPACK ]++) = 'i' ; // marker that an int16 follows, UBJSON style.
				m_buf_p[ eBITPACK ] = genfile::write_big_endian_integer( m_buf_p[ eBITPACK ], m_buf_end_p[ eBITPACK ], static_cast< int16_t >( 3 )) ; // 3 probs per samples

				m_bitpack_index = 0 ;
			}
		}
		void set_sample( std::size_t i ) {
			// nothing to do.
		}
		void set_number_of_entries( std::size_t n ) {
			if( n != 3 ) {
				m_buffer_validity[ eBITPACK ] = false ;
				throw genfile::BadArgumentError(
					"genfile::GenotypeProbabilityWriter::set_number_of_entries()",
					"n != 3"
				) ;
			}
		}
		void set_order_type( OrderType order_type ) {
			if( order_type != eUnorderedList ) {
				throw genfile::BadArgumentError(
					"genfile::GenotypeProbabilityWriter::set_order_type()",
					"order_type != eUnorderedList"
				) ;
			}
		}
		void operator()( genfile::MissingValue const value ) {
			double float_value = -1 ; // encode missingness as -1.
			operator()( float_value ) ;
		}

		void operator()( std::string& value ) {
			assert(0) ; // not allowed in this encoding.
		}
		void operator()( Integer const value ) {
			assert(0) ; // not allowed in this encoding.
		}

		void operator()( double const value ) {
			assert( ( m_buf_p[ eBITPACK ] + 4 ) <= m_buf_end_p[ eBITPACK ] ) ;
			float float_value = value ;
			if( float_value != float_value ) {
				float_value = -1 ; // encode missingness as -1.
			}
			m_buf_p[ eBITPACK ] = write_float( m_buf_p[ eBITPACK ], m_buf_end_p[ eBITPACK ], float_value ) ;
		}
		
		void finalise() {
			for( std::size_t i = 0; i < m_buffers.size(); ++i ) {
				m_buffers[ i ].resize( ( m_buf_p[ i ] ) - ( &( m_buffers[ i ][0] ) ) ) ;
			}
			
			if( m_buffer_validity[ eBITPACK ] ) {
				genfile::zlib_compress( &( m_buffers[ eBITPACK ][0] ), m_buf_p[eBITPACK], m_compression_buffer ) ;
			} else {
				assert(0) ;
			}
		}
		
		std::string get_encoding_name() const {
			if( m_buffer_validity[ eBITPACK ] ) {
				return "3float" ;
			} else {
				assert(0) ;
			}
		}
		
	private:
		std::vector< char >* m_compression_buffer ;
		std::vector< std::vector< char > > m_buffers ;
		std::size_t m_bitpack_index ;
		std::vector< bool > m_buffer_validity ;
		std::vector< char* > m_buf_p ;
		std::vector< char* > m_buf_end_p ;
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

void SQLiteGenotypesSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) {
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

void SQLiteGenotypesSNPDataSink::write_variant_data_impl(
	genfile::SNPIdentifyingData const& id_data,
	genfile::VariantDataReader& data_reader,
	Info const& info
) {
	if( m_data_i == m_snps.size() ) {
		flush_data( m_data_i ) ;
		m_data_i = 0 ;
	}
	m_snps[ m_data_i ] = id_data ;
	GenotypeProbabilityWriter writer( &(m_data[ m_data_i ]) ) ;
	data_reader.get( "genotypes", writer ) ;
	writer.finalise() ;
	++m_data_i ;
}

void SQLiteGenotypesSNPDataSink::finalise_impl() {
	if( m_data_i > 0 ) {
		flush_data( m_data_i ) ;
		m_data_i = 0 ;
	}
	m_outputter->finalise() ;
}

void SQLiteGenotypesSNPDataSink::flush_data( std::size_t const data_count ) {
	db::Connection::ScopedTransactionPtr transaction = m_outputter->connection().open_transaction( 600 ) ;

#if DEBUG_SQLITEGENOTYPESNPDATASINK
	std::cerr << "Flushing " << data_count << " elements..." ;
	std::size_t max_data_size = 0 ;
#endif
	std::vector< db::Connection::RowId > variant_ids( data_count ) ;
	for( std::size_t i = 0; i < data_count; ++i ) {
		variant_ids[i] = m_outputter->get_or_create_variant( m_snps[i] ) ;
	}
#if DEBUG_SQLITEGENOTYPESNPDATASINK
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
#if DEBUG_SQLITEGENOTYPESNPDATASINK
		max_data_size = std::max( max_data_size, m_data[i].size() ) ;
#endif
	}
#if DEBUG_SQLITEGENOTYPESNPDATASINK
	std::cerr << "Done, max data size was " << max_data_size << ".\n" ;
#endif
}


