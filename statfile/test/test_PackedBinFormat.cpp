
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <utility>
#include <map>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <limits>
#include <boost/tuple/tuple.hpp>
#include "test_case.hpp"
#include "statfile/RFormatStatSource.hpp"
#include "statfile/PackedBinFormatStatSink.hpp"
#include "statfile/PackedBinFormatStatSource.hpp"
#include "statfile/endianness_utils.hpp"

void write_header_block( std::ostream& stream, uint32_t number_of_rows, uint32_t number_of_columns, std::string const& version_text, std::string const& descriptive_text ) {
	uint32_t header_block_size = 5 + 5 + 5 + 5 + version_text.size() + 5 + descriptive_text.size() + 5 ;
	// we will write number_of_columns 4-byte strings
	header_block_size += (9 * number_of_columns) ;
	stream.put( 'u' ) ;
	statfile::write_little_endian_integer( stream, header_block_size ) ;
	stream.put( 's' ) ;
	uint32_t size = version_text.size() ;
	statfile::write_little_endian_integer( stream, size ) ;
	stream.write( version_text.data(), size ) ;
	stream.put( 'u' ) ;
	statfile::write_little_endian_integer( stream, number_of_rows ) ;
	stream.put( 'u' ) ;
	statfile::write_little_endian_integer( stream, number_of_columns ) ;
	stream.put( 'u' ) ;
	statfile::write_little_endian_integer( stream, number_of_columns ) ;
	stream.put( 's' ) ;
	size = descriptive_text.size() ;
	statfile::write_little_endian_integer( stream, size ) ;
	stream.write( descriptive_text.data(), size ) ;

	// Write column headers.
	for( std::size_t i = 0; i < number_of_columns; ++i ) {
		stream.put( 's' ) ;
		size = 4 ;
		std::ostringstream ostr ;
		ostr << std::setw(4) << i ;
		statfile::write_little_endian_integer( stream, size ) ;
		stream.write( ostr.str().data(), 4 ) ;		
	}
}

typedef boost::tuple< uint32_t, uint32_t, std::string > data_t ;

using statfile::write_little_endian_integer ;

void construct_test_data(
	std::ostream & stream,
	uint32_t number_of_rows,
	uint32_t number_of_columns,
	std::string const& descriptive_text
) {
	std::string version_text = "packedbin_v1" ;
	write_header_block( stream, number_of_rows, number_of_columns, version_text, descriptive_text ) ;

	for( std::size_t i = 0; i < number_of_rows; ++i ) {
		for( std::size_t j = 0; j < number_of_columns; ++j ) {
			int c = (i+j) % 4 ;
			int32_t ivalue = (i+j) % 9 ;
			uint32_t uvalue = (i+j) % 9 ;
			double dvalue = (i+j) % 9 ;
			std::ostringstream str ;
			if( ivalue != 0 ) {
				str << ivalue ;
			}
			uint32_t str_size = str.str().size() ;
			std::string svalue = str.str() ;
			switch( c ) {
				case 0:
					if( ivalue != int32_t(0)) {
						write_little_endian_integer( stream, uint32_t(i) ) ;
						write_little_endian_integer( stream, uint32_t(j) ) ;
						stream.put( 'i' ) ;
						write_little_endian_integer( stream, ivalue ) ;
					}
					break ;
				case 1:
					if( uvalue != uint32_t(0)) {
						write_little_endian_integer( stream, uint32_t(i) ) ;
						write_little_endian_integer( stream, uint32_t(j) ) ;
						stream.put( 'u' ) ;
						write_little_endian_integer( stream, uvalue ) ;
					}
					break ;
				case 2:
					if( dvalue != 0.0 ) {
						write_little_endian_integer( stream, uint32_t(i) ) ;
						write_little_endian_integer( stream, uint32_t(j) ) ;
						stream.put( 'd' ) ;
						stream.write( reinterpret_cast< char const* >( &dvalue ), sizeof( dvalue )) ;
					}
					break ;
				case 3:
					if( svalue != "" ) {
						write_little_endian_integer( stream, uint32_t(i) ) ;
						write_little_endian_integer( stream, uint32_t(j) ) ;
						stream.put( 's' ) ;
						write_little_endian_integer( stream, str_size ) ;
						stream.write( svalue.data(), str_size ) ;
					}
					break ;
				default:
					assert( 0 ) ;
			}
		}
	}
	write_little_endian_integer( stream, std::numeric_limits< uint32_t> ::max() ) ;
	write_little_endian_integer( stream, std::numeric_limits< uint32_t> ::max() ) ;
}

void write_test_data(
	statfile::PackedBinFormatStatSink& sink,
	data_t const& data
) {
	// set the descriptive text
	sink.set_descriptive_text( data.get<2>() ) ;

	// set up the column names.
	for( std::size_t i = 0; i < data.get<1>() ; ++i ) {
		std::ostringstream ostr ;
		ostr << std::setw(4) << i ;
		sink | ostr.str() ;
	}

	for( std::size_t i = 0; i < data.get<0>(); ++i ) {
		for( std::size_t j = 0; j < data.get<1>(); ++j ) {
			int value = (i+j) % 9 ;
			int type = (i+j) % 4 ;
			std::ostringstream str ;
			if( value != 0 ) {
				str << value ;
			}
			switch( type ) {
				case 0:
					sink << int64_t( value ) ;
					break ;
				case 1:
					sink << int64_t( value ) ;
					break ;
				case 2:
					sink << double( value ) ;
					break ;
				case 3:
					sink << str.str() ;
					break ;
				default:
					assert( 0 ) ;
			}
		}
		
		sink << statfile::end_row() ;
	}
}

AUTO_TEST_CASE( test_PackedBinFormatStatSink ) 
{
	std::vector< data_t > test_data ;
	test_data.push_back( boost::make_tuple( 0, 0, "" )) ;
	test_data.push_back( boost::make_tuple( 0, 10, "" )) ;
	test_data.push_back( boost::make_tuple( 5, 5, "A description" )) ;
	test_data.push_back( boost::make_tuple( 10, 10, "" )) ;
	test_data.push_back( boost::make_tuple( 2, 2, "This is some text" )) ;
	test_data.push_back( boost::make_tuple( 5, 2, "description" )) ;

	for( std::size_t di = 0; di < test_data.size(); ++di ) {
		std::cerr << "test_PackedBinFormatStatSink: data item " << di << ".\n" ;
		
		std::ostringstream str1 ;
		construct_test_data( str1, test_data[di].get<0>(), test_data[di].get<1>(), test_data[di].get<2>() ) ;
		std::string expected = str1.str() ;

		std::vector< char > buffer( 10000, char(0) ) ;
		std::auto_ptr< std::ostringstream > str2( new std::ostringstream ) ;
		str2->rdbuf()->pubsetbuf( &buffer[0], buffer.size() ) ;
		std::auto_ptr< std::ostream > str2b( str2.release() ) ;
		{
			statfile::PackedBinFormatStatSink sink( str2b ) ;
			write_test_data( sink, test_data[di] ) ;
		} // sink must be deconstructed before comparison.

		// Check output data is terminated by something != 0.
		std::string actual( buffer.begin(), buffer.end() );
		std::size_t pos = actual.find_last_not_of( std::string( 1, char(0) )) ;
		TEST_ASSERT( pos != std::string::npos ) ;
		actual.resize( pos + 1 ) ;

		if( expected.size() != actual.size() ) {
			std::cerr << "Expected size: " << expected.size() << ", actual size: " << actual.size() << ".\n" ;
		}
		TEST_ASSERT( expected.size() == actual.size() ) ;
		if( expected != actual ) {
			std::pair< std::string::const_iterator, std::string::const_iterator > mismatch_i = std::mismatch( expected.begin(), expected.end(), actual.begin() ) ;
			std::size_t location = (mismatch_i.first - expected.begin()) ;
			std::cerr << "First mismatch at position "
				<< location
				<< ": expected 0x"
				<< std::hex << int( *(mismatch_i.first))
				<< ", got 0x"
				<< std::hex << int( *(mismatch_i.second))
				<< ".\n" ;
		}
		TEST_ASSERT( expected == actual ) ;
	}
}

AUTO_TEST_CASE( test_PackedBinFormatStatSource ) {
	std::vector< data_t > test_data ;
	test_data.push_back( boost::make_tuple( 0, 0, "" )) ;
	test_data.push_back( boost::make_tuple( 0, 10, "" )) ;
	test_data.push_back( boost::make_tuple( 5, 5, "A description" )) ;
	test_data.push_back( boost::make_tuple( 10, 10, "" )) ;
	test_data.push_back( boost::make_tuple( 2, 2, "This is some text" )) ;
	test_data.push_back( boost::make_tuple( 5, 2, "description" )) ;


	for( std::size_t di = 0; di < test_data.size(); ++di ) {
		std::cerr << "test_PackedBinFormatStatSource: data item " << di << ".\n" ;
		
		// Write the test data via a stat sink into a buffer
		std::vector< char > buffer( 10000 ) ;
		std::auto_ptr< std::ostringstream > oStr( new std::ostringstream ) ;
		oStr->rdbuf()->pubsetbuf( &buffer[0], buffer.size() ) ;
		std::auto_ptr< std::ostream > oStr2( oStr.release() ) ;
		{
			statfile::PackedBinFormatStatSink sink( oStr2 ) ;
			write_test_data( sink, test_data[di] ) ;
		} // sink must be deconstructed before data is fully written.

		// Check output data is terminated by something != 0.
		std::string actual( buffer.begin(), buffer.end() );
		std::size_t pos = actual.rfind( char(0) ) ;
		TEST_ASSERT( pos > 1 ) ;
		actual.resize( pos - 1 ) ;

		// Now the buffer should hold the written data.
		std::auto_ptr< std::istringstream > instr( new std::istringstream() ) ;
		instr->rdbuf()->pubsetbuf( &buffer[0], actual.size() ) ;
		std::auto_ptr< std::istream > instr2( instr.release() ) ;
		statfile::PackedBinFormatStatSource source( instr2 ) ;

		TEST_ASSERT( source.number_of_rows() == test_data[di].get<0>() ) ;
		TEST_ASSERT( source.number_of_columns() == test_data[di].get<1>() ) ;
		TEST_ASSERT( source.get_descriptive_text() == test_data[di].get<2>() ) ;
		for( std::size_t i = 0; i < source.number_of_columns(); ++i ) {
			std::ostringstream oStream ;
			oStream << std::setw(4) << i ;
			TEST_ASSERT( source.column_name(i) == oStream.str() ) ;
		}
		
		for( std::size_t i = 0; i < source.number_of_rows(); ++i ) {
			for( std::size_t j = 0; j < source.number_of_columns(); ++j ) {
				int c = (i + j) % 4 ;
				int expected_value = (i+j) % 9 ;
				if( c == 0 ) {
					int32_t value ;
					source >> value ;
					TEST_ASSERT( value == int32_t( expected_value )) ;
				}
				else if( c == 1 ) {
					uint32_t value ;
					source >> value ;
					std::cerr << value << " : expected " << expected_value << ".\n" ; 
					TEST_ASSERT( value == uint32_t( expected_value )) ;
				}
				else if( c == 2 ) {
					double value ;
					source >> value ;
					std::cerr << "i = " << i << ", j = " << j << ", value = " << value << " expected = " << expected_value << ".\n" ; 
					TEST_ASSERT( value == double( expected_value )) ;
				}
				else if( c == 3 ) {
					std::string value ;
					source >> value ;
					std::ostringstream oStream ;
					if( expected_value != 0 ) {
						oStream << expected_value ;
					}
					TEST_ASSERT( value == oStream.str() ) ;
				}
			}
			
			source >> statfile::end_row() ;
		}
		
		TEST_ASSERT( source ) ;
		double dummy_value ;
		TEST_ASSERT( !( source >> dummy_value )) ;
	}
}

