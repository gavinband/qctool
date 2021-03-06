
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
#include "genfile/bgen/bgen.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/zlib.hpp"
#include "genfile/endianness_utils.hpp"
#include "qcdb/DBOutputter.hpp"
#include "components/SNPOutputComponent/SQLiteGenotypesSNPDataSink.hpp"
#include "qcdb/store_samples_in_db.hpp"

// #define DEBUG_SQLITEGENOTYPESNPDATASINK 1

SQLiteGenotypesSNPDataSink::SQLiteGenotypesSNPDataSink(
	qcdb::DBOutputter::UniquePtr outputter,
	genfile::CohortIndividualSource const& samples
):
	m_genotype_field( ":genotypes:" ),
	m_intensity_field( ":intensities:" ),
	m_outputter( outputter ),
	m_samples( samples ),
	m_genotype_data_i( 0 ),
	m_intensity_data_i( 0 )
{
	genfile::db::Connection::ScopedTransactionPtr transaction = m_outputter->connection().open_transaction( 2400 ) ;
	m_outputter->connection().run_statement(
		"CREATE TABLE IF NOT EXISTS Genotype ( "
		"  analysis_id INT NOT NULL REFERENCES Analysis( id ), "
		"  variant_id INT NOT NULL REFERENCES Variant( id ), "
		"  N INT NOT NULL,"
		"  data BLOB NOT NULL"
		")"
	) ;
	m_outputter->connection().run_statement(
		"CREATE View IF NOT EXISTS GenotypeView AS "
		"SELECT analysis_id, A.name AS analysis, variant_id, chromosome, position, rsid, alleleA, alleleB, quote( data ) "
		"FROM Genotype H "
		"INNER JOIN Analysis A "
		"ON A.id == H.analysis_id "
		"INNER JOIN Variant V "
		"ON V.id == H.variant_id"
	) ;
	m_insert_genotype_stmnt = m_outputter->connection().get_statement(
		"INSERT INTO Genotype ( analysis_id, variant_id, N, data ) VALUES( ?, ?, ?, ? ) ;"
	) ;

	m_outputter->connection().run_statement(
		"CREATE TABLE IF NOT EXISTS Intensity ( "
		"  analysis_id INT NOT NULL REFERENCES Analysis( id ), "
		"  variant_id INT NOT NULL REFERENCES Variant( id ), "
		"  N INT NOT NULL,"
		"  data BLOB NOT NULL"
		")"
	) ;
	
	m_outputter->connection().run_statement(
		"CREATE View IF NOT EXISTS IntensityView AS "
		"SELECT analysis_id, A.name AS analysis, variant_id, chromosome, position, rsid, alleleA, alleleB, quote( data ) "
		"FROM Intensity H "
		"INNER JOIN Analysis A "
		"ON A.id == H.analysis_id "
		"INNER JOIN Variant V "
		"ON V.id == H.variant_id"
	) ;
	m_insert_intensity_stmnt = m_outputter->connection().get_statement(
		"INSERT INTO Intensity ( analysis_id, variant_id, N, data ) VALUES( ?, ?, ?, ? ) ;"
	) ;

	m_genotype_snps.resize( 1000 ) ;
	m_genotype_data.resize( 1000 ) ;

	m_intensity_snps.resize( 1000 ) ;
	m_intensity_data.resize( 1000 ) ;
}

std::string SQLiteGenotypesSNPDataSink::get_spec() const {
	return "SQLiteGenotypesSNPDataSink" ;
}

void SQLiteGenotypesSNPDataSink::set_genotype_field( std::string field ) {
	m_genotype_field = field ;
}

void SQLiteGenotypesSNPDataSink::set_intensity_field( std::string field ) {
	m_intensity_field = field ;
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

	struct PositiveFloatWriter: public genfile::VariantDataReader::PerSampleSetter {
		enum Encoding { eBITPACK = 0 } ; // must take contiguous values starting at zero.
		
		 PositiveFloatWriter( std::vector< uint8_t >* compression_buffer, std::size_t number_of_entries ):
			m_compression_buffer( compression_buffer ),
			m_expected_number_if_entries( number_of_entries ),
			m_buffer_validity( 1, true )
		{
			assert( compression_buffer != 0 ) ;
		}

		~PositiveFloatWriter() throw() {}
		
		void initialise( std::size_t n, std::size_t nAlleles ) {
			assert( nAlleles == 2 ) ;

			m_buffers.resize(1) ;
			m_buf_p.resize(1) ;
			m_buf_end_p.resize(1) ;

			{
				m_buffers[ eBITPACK ].clear() ;
				std::string const encoding_name = "floatarray" ;
				m_buffers[ eBITPACK ].resize( 3 + encoding_name.size() + 5 + 3 + ( n * 4 * m_expected_number_if_entries ) ) ; // name, plus m_expected_number_if_entries floats per sample.
				m_buf_p[ eBITPACK ] = &(m_buffers[ eBITPACK ][0]) ;
				m_buf_end_p[ eBITPACK ] = &(m_buffers[ eBITPACK ][0]) + m_buffers[ eBITPACK ].size() ;
				m_buffers[ eBITPACK ][0] = 's' ; // marker that a string follows, UBJSON style.
				++m_buf_p[ eBITPACK ] ;
				m_buf_p[ eBITPACK ] = genfile::write_big_endian_integer( m_buf_p[ eBITPACK ], m_buf_end_p[ eBITPACK ], static_cast< int16_t >( encoding_name.size() )) ;
				m_buf_p[ eBITPACK ] = std::copy( encoding_name.begin(), encoding_name.end(), m_buf_p[ eBITPACK ] ) ;

				*(m_buf_p[ eBITPACK ]++) = 'I' ; // marker that an int32 follows, UBJSON style.
				m_buf_p[ eBITPACK ] = genfile::write_big_endian_integer( m_buf_p[ eBITPACK ], m_buf_end_p[ eBITPACK ], static_cast< int32_t >( n )) ; // number of samples

				*(m_buf_p[ eBITPACK ]++) = 'i' ; // marker that an int16 follows, UBJSON style.
				m_buf_p[ eBITPACK ] = genfile::write_big_endian_integer( m_buf_p[ eBITPACK ], m_buf_end_p[ eBITPACK ], static_cast< int16_t >( m_expected_number_if_entries )) ;

				m_bitpack_index = 0 ;
			}
		}
		
		bool set_sample( std::size_t i ) {
			// nothing to do.
			return true ;
		}
		void set_number_of_entries( uint32_t, std::size_t n, OrderType const order_type, ValueType const value_type ) {
			if( n != m_expected_number_if_entries ) {
				m_buffer_validity[ eBITPACK ] = false ;
				throw genfile::BadArgumentError(
					"genfile:: PositiveFloatWriter::set_number_of_entries()",
					"n != " + genfile::string_utils::to_string( m_expected_number_if_entries )
				) ;
			}
			if( order_type != genfile::eOrderedList && order_type != genfile::ePerUnorderedGenotype && order_type != genfile::ePerAllele ) {
				throw genfile::BadArgumentError(
					"genfile:: PositiveFloatWriter::set_order_type()",
					"order_type",
					"Expected order_type to be eUnorderedList or ePerUnorderedGenotype"
				) ;
			}
			if( value_type != genfile::eProbability && value_type != genfile::eUnknownValueType ) {
				throw genfile::BadArgumentError(
					"genfile:: PositiveFloatWriter::set_order_type()",
					"value_type",
					"Expected value_type == eProbability"
				) ;
			}
		}
		void set_value( std::size_t i, genfile::MissingValue const value ) {
			double float_value = -1 ; // encode missingness as -1.
			set_value( i, float_value ) ;
		}

		void set_value( std::size_t i, double const value ) {
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
				return "floatarray" ;
			} else {
				assert(0) ;
			}
		}
		
	private:
		std::vector< uint8_t >* m_compression_buffer ;
		std::size_t const m_expected_number_if_entries ;
		std::vector< std::vector< uint8_t > > m_buffers ;
		std::size_t m_bitpack_index ;
		std::vector< bool > m_buffer_validity ;
		std::vector< uint8_t* > m_buf_p ;
		std::vector< uint8_t* > m_buf_end_p ;
	private:
		uint8_t* write_float( uint8_t* buf_p, uint8_t* buf_end_p, float const& value ) const {
			uint8_t const* v_p = reinterpret_cast< uint8_t const* >( &value ) ;
			*(buf_p++) = *v_p++ ;
			*(buf_p++) = *v_p++ ;
			*(buf_p++) = *v_p++ ;
			*(buf_p++) = *v_p++ ;
			return buf_p ;
		}
	} ;
	
	// This class handles genotypes set as single int (encoding dosage)
	// as a pair of ints (as in vcf GT field)
	// or as three probabilities.
	// In all cases we convert to three probabilities.
	// Unlike GEN files, we encode missing probabilities as missing data not zeroes.
	struct GenotypeMunger: public genfile::VariantDataReader::PerSampleSetter {
		GenotypeMunger( PositiveFloatWriter& writer ):
			m_writer( writer ),
			m_number_of_entries( 0 ),
			m_order_type( genfile::eUnknownOrderType ),
			m_value_type( genfile::eUnknownValueType ),
			m_allele_probs( 2 ),
			m_genotype_probs( 3 )
		{}
		
		~GenotypeMunger() throw() {}
		
		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
			m_writer.initialise( nSamples, nAlleles ) ;
		}
		bool set_sample( std::size_t i ) {
			m_writer.set_sample(i) ;
			return true ;
		}
		void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
			if(
				!(
					( order_type == genfile::ePerSample && n == 1 )
					|| (( order_type == genfile::ePerOrderedHaplotype || order_type == genfile::ePerUnorderedHaplotype ) && n == 2 )
					|| ( order_type == genfile::ePerUnorderedGenotype && n == 3 )
				)
			) {
				throw genfile::BadArgumentError(
					"genfile:: GenotypeMunger::set_number_of_entries()",
					"Expected one dosage, two genotype calls, or three probabilities."
				) ;
			}
			m_number_of_entries = n ;
			m_order_type = order_type ;
			m_value_type = value_type ;

			m_writer.set_number_of_entries( 2, 3, genfile::ePerUnorderedGenotype, genfile::eProbability ) ;
			m_allele_prob_i = 0 ;
			m_genotype_prob_i = 0 ;
			m_genotype_probs[0] = m_genotype_probs[1] = m_genotype_probs[2] = -1.0 ; // encode missingness as -1
			m_allele_probs[0] = m_allele_probs[1] = -1.0 ;  // encode missingness as -1
		}

		void set_value( std::size_t i, genfile::MissingValue const value ) {
			double float_value = -1 ; // encode missingness as -1.
			set_value( i, float_value ) ;
		}
		void set_value( std::size_t i, Integer const value ) {
			set_value( i, double( value )) ;
		}
		void set_value( std::size_t, double const value ) {
			// TODO: fix this to handle order_type and value_type properly.
			switch( m_order_type ) {
				case genfile::ePerSample:
					// genotype dosage information.  Just write the probabilities (two zeros, one one) directly.
					// If -1 comes in (via set_value( MissingValue )) then all genotype probs will remain as -1.
					if( value >= 0.0 ) {
						for( std::size_t g = 0; g < 3; ++g ) {
							m_genotype_probs[g] = ( ( value == g ) ? 1.0 : 0.0 )  ;
						}
					}
					go() ;
					break ;
				case genfile::ePerOrderedHaplotype:
				case genfile::ePerUnorderedHaplotype:
					// GT-style pair of genotypes.  Need to store 'em.
					m_allele_probs[ m_allele_prob_i++ ] = value ;
					if( m_allele_prob_i == 2 ) {
						// compute and store the genotype probs.
						if( m_allele_probs[0] != -1.0 && m_allele_probs[1] != -1.0 ) {
							for( std::size_t g = 0; g < 3; ++g ) {
								m_genotype_probs[g] = (( m_allele_probs[0] + m_allele_probs[1] ) == g ) ? 1.0 : 0.0 ;
							}
						}
						go() ;
					}
					break ;
				case genfile::ePerUnorderedGenotype:
					m_genotype_probs[ m_genotype_prob_i++ ] = value ;
					if( m_genotype_prob_i == 3 ) {
						// store the genotype probs
						go() ;
					}
					break ;
				default:
					assert(0) ;
			}
		}
		
		void finalise() {}
		
		void go() {
			if(
				( m_genotype_probs[0] == 0.0 && m_genotype_probs[1] == 0.0 && m_genotype_probs[2] == 0.0 ) // gen-style setting of all probs to 0
			 	|| ( m_genotype_probs[0] == -1.0 || m_genotype_probs[1] == -1.0 || m_genotype_probs[2] == -1.0 ) // all probs set as -1's.
			) {
				for( int g = 0; g < 3; ++g ) {
					m_writer.set_value( g, genfile::MissingValue() ) ;
				}
			} else {
				for( std::size_t g = 0; g < 3; ++g ) {
					m_writer.set_value( g, m_genotype_probs[g] ) ;
				}
			}
		}
	private:
		PositiveFloatWriter& m_writer ;
		std::size_t m_number_of_entries ;
		OrderType m_order_type ;
		ValueType m_value_type ;
		std::vector< double > m_allele_probs ;
		std::size_t m_allele_prob_i ;
		std::vector< double > m_genotype_probs ;
		std::size_t m_genotype_prob_i ;
	} ;

	
}

void SQLiteGenotypesSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) {
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

void SQLiteGenotypesSNPDataSink::write_variant_data_impl(
	genfile::VariantIdentifyingData const& id_data,
	genfile::VariantDataReader& data_reader,
	Info const& info
) {
	if( data_reader.supports( m_genotype_field ) ) {
		if( m_genotype_data_i == m_genotype_snps.size() ) {
			flush_genotype_data( m_genotype_data_i ) ;
			m_genotype_data_i = 0 ;
		}

		m_genotype_snps[ m_genotype_data_i ] = id_data ;
		PositiveFloatWriter writer( &(m_genotype_data[ m_genotype_data_i ]), 3 ) ;
		GenotypeMunger munger( writer ) ;
		data_reader.get( m_genotype_field, munger ) ;
		writer.finalise() ;
		++m_genotype_data_i ;
	}
	
	if( data_reader.supports( m_intensity_field )) {
		if( m_intensity_data_i == m_intensity_snps.size() ) {
			flush_intensity_data( m_intensity_data_i ) ;
			m_intensity_data_i = 0 ;
		}

		m_intensity_snps[ m_intensity_data_i ] = id_data ;
		PositiveFloatWriter writer( &(m_intensity_data[ m_intensity_data_i ]), 2 ) ;
		data_reader.get( ":intensities:", writer ) ;
		writer.finalise() ;
		++m_intensity_data_i ;
	}
}

void SQLiteGenotypesSNPDataSink::finalise_impl() {
	if( m_genotype_data_i > 0 ) {
		flush_genotype_data( m_genotype_data_i ) ;
		m_genotype_data_i = 0 ;
	}

	if( m_intensity_data_i > 0 ) {
		flush_intensity_data( m_intensity_data_i ) ;
		m_intensity_data_i = 0 ;
	}
	
	{
		genfile::db::Connection::ScopedTransactionPtr transaction = m_outputter->connection().open_transaction( 2400 ) ;
		m_outputter->connection().run_statement( "CREATE INDEX IF NOT EXISTS IntensityVariantIndex ON Intensity ( variant_id )" ) ;
		m_outputter->connection().run_statement( "CREATE INDEX IF NOT EXISTS IntensityAnalysisVariantIndex ON Intensity ( analysis_id, variant_id )" ) ;
		m_outputter->connection().run_statement( "CREATE INDEX IF NOT EXISTS GenotypeVariantIndex ON Genotype ( variant_id )" ) ;
		m_outputter->connection().run_statement( "CREATE INDEX IF NOT EXISTS GenotypeAnalysisVariantIndex ON Genotype ( analysis_id, variant_id )" ) ;
	}

	m_outputter->finalise() ;
}

void SQLiteGenotypesSNPDataSink::flush_genotype_data( std::size_t const data_count ) {
	genfile::db::Connection::ScopedTransactionPtr transaction = m_outputter->connection().open_transaction( 2400 ) ;

#if DEBUG_SQLITEGENOTYPESNPDATASINK
	std::cerr << "Flushing " << data_count << " genotypes..." ;
	std::size_t max_data_size = 0 ;
#endif
	std::vector< genfile::db::Connection::RowId > variant_ids( data_count ) ;
	for( std::size_t i = 0; i < data_count; ++i ) {
		variant_ids[i] = m_outputter->get_or_create_variant( m_genotype_snps[i] ) ;
	}
#if DEBUG_SQLITEGENOTYPESNPDATASINK
	std::cerr << "stored variants..." ;
#endif
	for( std::size_t i = 0; i < data_count; ++i ) {
		uint8_t const* buffer = &(m_genotype_data[i][0]) ;
		uint8_t const* end_buffer = &(m_genotype_data[i][0]) + m_genotype_data[i].size() ;

		m_insert_genotype_stmnt
			->bind( 1, m_outputter->analysis_id() )
			.bind( 2, variant_ids[i] )
			.bind( 3, int64_t( m_samples.get_number_of_individuals() ) )
			.bind( 4, buffer, end_buffer )
			.step() ;
		m_insert_genotype_stmnt->reset() ;
#if DEBUG_SQLITEGENOTYPESNPDATASINK
		max_data_size = std::max( max_data_size, m_genotype_data[i].size() ) ;
#endif
	}
#if DEBUG_SQLITEGENOTYPESNPDATASINK
	std::cerr << "Done, max data size was " << max_data_size << ".\n" ;
#endif
}

void SQLiteGenotypesSNPDataSink::flush_intensity_data( std::size_t const data_count ) {
	genfile::db::Connection::ScopedTransactionPtr transaction = m_outputter->connection().open_transaction( 1200 ) ;

#if DEBUG_SQLITEGENOTYPESNPDATASINK
	std::cerr << "Flushing " << data_count << " intensities..." ;
	std::size_t max_data_size = 0 ;
#endif
	std::vector< genfile::db::Connection::RowId > variant_ids( data_count ) ;
	for( std::size_t i = 0; i < data_count; ++i ) {
		variant_ids[i] = m_outputter->get_or_create_variant( m_intensity_snps[i] ) ;
	}
#if DEBUG_SQLITEGENOTYPESNPDATASINK
	std::cerr << "stored variants..." ;
#endif
	for( std::size_t i = 0; i < data_count; ++i ) {
		uint8_t const* buffer = &(m_intensity_data[i][0]) ;
		uint8_t const* end_buffer = &(m_intensity_data[i][0]) + m_intensity_data[i].size() ;

		m_insert_intensity_stmnt
			->bind( 1, m_outputter->analysis_id() )
			.bind( 2, variant_ids[i] )
			.bind( 3, int64_t( m_samples.get_number_of_individuals() ) )
			.bind( 4, buffer, end_buffer )
			.step() ;
		m_insert_intensity_stmnt->reset() ;
#if DEBUG_SQLITEGENOTYPESNPDATASINK
		max_data_size = std::max( max_data_size, m_intensity_data[i].size() ) ;
#endif
	}
#if DEBUG_SQLITEGENOTYPESNPDATASINK
	std::cerr << "Done, max data size was " << max_data_size << ".\n" ;
#endif
}



