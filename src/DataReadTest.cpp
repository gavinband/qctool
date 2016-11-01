
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <boost/bind.hpp>
#include "genfile/FileUtils.hpp"
#include "genfile/Error.hpp"
#include "genfile/zlib.hpp"
#include "genfile/endianness_utils.hpp"
#include "genfile/vcf/get_set.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "DataReadTest.hpp"


void DataReadTest::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Test / benchmark options" ) ;
	options[ "-read-test" ]
		.set_description( "Don't do anything; just read the data." ) ;
}

DataReadTest::DataReadTest() {}

void DataReadTest::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps_read = 0 ;
	m_data.resize( m_number_of_samples ) ;
}

void DataReadTest::processed_snp( genfile::VariantIdentifyingData const& snp, genfile::VariantDataReader::SharedPtr data_reader ) {
	std::vector< std::string > fields ;
	void(std::vector<std::string>::*push_back)(std::string const&) = &std::vector<std::string>::push_back ;
	data_reader->get_supported_specs( boost::bind( push_back, &fields, _1 )) ;
	
	for( std::size_t field_i = 0; field_i < fields.size(); ++field_i ) {
		std::string const field = fields[ field_i ] ;
		// Make sure we've got these fields in Entity
		genfile::vcf::VectorSetter setter( m_data ) ;
		data_reader->get( field, setter ) ;
	}
}

void DataReadTest::end_processing_snps() {
	// nothing to do.
}

