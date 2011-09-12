#include <iostream>
#include <utility>
#include <map>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include "test_case.hpp"
#include "statfile/RFormatStatSource.hpp"
#include "statfile/BinFormatStatSink.hpp"
#include "statfile/BinFormatStatSource.hpp"
#include "statfile/endianness_utils.hpp"

namespace {
	void write_header_block( std::ostream& stream, uint32_t number_of_rows, uint32_t number_of_columns, uint32_t number_of_id_columns, std::string const& version_text, std::string const& descriptive_text ) {
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
		statfile::write_little_endian_integer( stream, number_of_id_columns ) ; // number of id columns.
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
		std::string version_text = "statfilebin_v1" ;
		std::size_t number_of_id_columns = number_of_columns / 2;
		write_header_block( stream, number_of_rows, number_of_columns, number_of_id_columns, version_text, descriptive_text ) ;

		for( std::size_t i = 0; i < number_of_rows; ++i ) {
			for( std::size_t j = 0; j < number_of_id_columns; ++j ) {
				int c = (i+j) % 4 ;
				double value = std::exp( i+j ) ;
				std::ostringstream str ;
				str << ( i+j ) ;
				uint32_t str_size = str.str().size() ;
				switch( c ) {
					case 0:
						stream.put( 'i' ) ;
						write_little_endian_integer( stream, int32_t( i+j ) ) ;
						break ;
					case 1:
						stream.put( 'u' ) ;
						write_little_endian_integer( stream, uint32_t( i+j ) ) ;
						break ;
					case 2:
						stream.put( 'd' ) ;
						stream.write( reinterpret_cast< char const* >( &value ), sizeof( value )) ;
						break ;
					case 3:
						stream.put( 's' ) ;
						write_little_endian_integer( stream, str_size ) ;
						stream.write( str.str().data(), str_size ) ;
						break ;
					default:
						assert( 0 ) ;
				}
			}
		
			for( std::size_t j = number_of_id_columns; j < number_of_columns; ++j ) {
				double value = std::exp( i+j ) ;
				stream.write( reinterpret_cast< char const* >( &value ), sizeof( value )) ;
			}
		} 
	}

	void write_test_data(
		statfile::BinFormatStatSink& sink,
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

		std::size_t number_of_id_columns = data.get<1>() / 2;

		for( std::size_t i = 0; i < data.get<0>(); ++i ) {
			for( std::size_t j = 0; j < number_of_id_columns; ++j ) {
				int c = (i+j) % 4 ;
				std::ostringstream str ;
				str << (i+j) ;
				switch( c ) {
					case 0:
						sink << int64_t( i+j ) ;
						break ;
					case 1:
						sink << int64_t( i+j ) ;
						break ;
					case 2:
						sink << double( std::exp( i+j ) ) ;
						break ;
					case 3:
						sink << str.str() ;
						break ;
					default:
						assert( 0 ) ;
				}
			}
		
			for( std::size_t j = number_of_id_columns; j < data.get<1>(); ++j ) {
				sink << double( std::exp( i+j ) ) ;
			}
			sink << statfile::end_row() ;
		}
	}
}

AUTO_TEST_CASE( test_BinFormatStatSink ) 
{
	std::vector< data_t > test_data ;
	test_data.push_back( boost::make_tuple( 0, 0, "" )) ;
	test_data.push_back( boost::make_tuple( 0, 10, "" )) ;
	test_data.push_back( boost::make_tuple( 5, 5, "A description" )) ;
	test_data.push_back( boost::make_tuple( 10, 10, "" )) ;
	test_data.push_back( boost::make_tuple( 2, 2, "This is some text" )) ;
	test_data.push_back( boost::make_tuple( 5, 2, "description" )) ;

	for( std::size_t di = 0; di < test_data.size(); ++di ) {
		std::cerr << "test_BinFormatStatSink: data item " << di << ".\n" ;
		
		std::ostringstream str1 ;
		construct_test_data( str1, test_data[di].get<0>(), test_data[di].get<1>(), test_data[di].get<2>() ) ;
		std::string expected = str1.str() ;

		std::vector< char > buffer( 10000 ) ;
		std::auto_ptr< std::ostringstream > str2( new std::ostringstream ) ;
		str2->rdbuf()->pubsetbuf( &buffer[0], buffer.size() ) ;
		std::ostringstream* dodgy_ptr = str2.get() ;
		std::auto_ptr< std::ostream > str2b( str2.release() ) ;
		std::size_t data_size ;
		{
			uint32_t number_of_id_columns = test_data[di].get<1>() / 2 ;
			statfile::BinFormatStatSink sink( str2b, number_of_id_columns ) ;
			write_test_data( sink, test_data[di] ) ;
			data_size = dodgy_ptr->tellp() ;
		} // sink must be deconstructed before comparison.

		std::string actual( buffer.begin(), buffer.begin() + data_size ) ;

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

AUTO_TEST_CASE( test_BinFormatStatSource ) {
	std::vector< data_t > test_data ;
	test_data.push_back( boost::make_tuple( 0, 0, "" )) ;
	test_data.push_back( boost::make_tuple( 0, 10, "" )) ;
	test_data.push_back( boost::make_tuple( 5, 5, "A description" )) ;
	test_data.push_back( boost::make_tuple( 10, 10, "" )) ;
	test_data.push_back( boost::make_tuple( 2, 2, "This is some text" )) ;
	test_data.push_back( boost::make_tuple( 5, 2, "description" )) ;


	for( std::size_t di = 0; di < test_data.size(); ++di ) {
		std::cerr << "test_BinFormatStatSource: data item " << di << ".\n" ;
		
		// Write the test data via a stat sink into a buffer
		std::vector< char > buffer( 10000 ) ;
		std::auto_ptr< std::ostringstream > oStr( new std::ostringstream ) ;
		oStr->rdbuf()->pubsetbuf( &buffer[0], buffer.size() ) ;
		std::ostringstream* dodgy_ptr = oStr.get() ;
		std::auto_ptr< std::ostream > oStr2( oStr.release() ) ;
		std::size_t data_size ;
		{
			uint32_t number_of_id_columns = test_data[di].get<1>() / 2 ;
			statfile::BinFormatStatSink sink( oStr2, number_of_id_columns ) ;
			write_test_data( sink, test_data[di] ) ;
			data_size = dodgy_ptr->tellp() ;
		} // sink must be deconstructed before data is fully written.

		// Now the buffer should hold the written data.
		std::auto_ptr< std::istringstream > instr( new std::istringstream() ) ;
		instr->rdbuf()->pubsetbuf( &buffer[0], data_size ) ;
		std::auto_ptr< std::istream > instr2( instr.release() ) ;
		statfile::BinFormatStatSource source( instr2 ) ;

		TEST_ASSERT( source.number_of_rows() == test_data[di].get<0>() ) ;
		TEST_ASSERT( source.number_of_columns() == test_data[di].get<1>() ) ;
		TEST_ASSERT( source.get_descriptive_text() == test_data[di].get<2>() ) ;
		for( std::size_t i = 0; i < source.number_of_columns(); ++i ) {
			std::ostringstream oStream ;
			oStream << std::setw(4) << i ;
			TEST_ASSERT( source.column_name(i) == oStream.str() ) ;
		}
		
		for( std::size_t i = 0; i < source.number_of_rows(); ++i ) {
			std::size_t number_of_id_columns = source.number_of_columns() / 2 ;
			for( std::size_t j = 0; j < number_of_id_columns; ++j ) {
				int c = (i + j) % 4 ;
				if( c == 0 ) {
					int32_t value ;
					source >> value ;
					TEST_ASSERT( value == int32_t( i+j )) ;
				}
				else if( c == 1 ) {
					uint32_t value ;
					source >> value ;
					TEST_ASSERT( value == uint32_t( i+j )) ;
				}
				else if( c == 2 ) {
					double value ;
					source >> value ;
					TEST_ASSERT( value == double( std::exp( i+j ))) ;
				}
				else if( c == 3 ) {
					std::string value ;
					source >> value ;
					std::ostringstream oStream ;
					oStream << (i+j) ;
					TEST_ASSERT( value == oStream.str() ) ;
				}
			}
			
			for( std::size_t j = number_of_id_columns; j < source.number_of_columns(); ++j ) {
				double value ;
				source >> value ;
				TEST_ASSERT( value == std::exp( i+j )) ;
			}

			source >> statfile::end_row() ;
			TEST_ASSERT( source ) ;
		}
	}
}
