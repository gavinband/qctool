#include <iostream>
#include <boost/bind.hpp>
#include "genfile/FileUtils.hpp"
#include "db/SQLite3Connection.hpp"
#include "db/SQLStatement.hpp"
#include "genfile/Error.hpp"
#include "genfile/zlib.hpp"
#include "genfile/endianness_utils.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "DataReadTest.hpp"


void DataReadTest::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Test / benchmark options" ) ;
	options[ "-read-test" ]
		.set_description( "Don't do anything; just read the data." ) ;
}

DataReadTest::DataReadTest() {}

void DataReadTest::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps = number_of_snps ;
	m_number_of_snps_read = 0 ;
	m_data.resize( m_number_of_samples ) ;
}

void DataReadTest::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	std::vector< std::string > fields ;
	data_reader.get_supported_specs( boost::bind( &std::vector< std::string >::push_back, &fields, _1 )) ;
	
	for( std::size_t field_i = 0; field_i < fields.size(); ++field_i ) {
		std::string const field = fields[ field_i ] ;
		// Make sure we've got these fields in Entity
		data_reader.get( field, genfile::VariantDataReader::set( m_data ) ) ;
	}
}

void DataReadTest::end_processing_snps() {
	// nothing to do.
}

