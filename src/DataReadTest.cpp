
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <boost/bind.hpp>
#include "genfile/FileUtils.hpp"
#include "genfile/Error.hpp"
#include "genfile/VariantDataReader.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "DataReadTest.hpp"


namespace {
	struct DataSetter: public genfile::VariantDataReader::PerSampleSetter {
		DataSetter(
			std::vector< uint32_t >& ploidy,
			std::vector< std::vector< double > >& doubles,
			std::vector< std::vector< int64_t > >& ints,
			std::vector< std::vector< std::string > >& strings
		):
			m_sample_i(0),
			m_ploidy( ploidy ),
			m_doubles( doubles ),
			m_ints( ints ),
			m_strings( strings )
		{}

		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
		}

		bool set_sample( std::size_t n ) {
			m_sample_i = n ;
			return true ;
		}

		inline void set_number_of_entries(
			uint32_t ploidy, std::size_t n,
			genfile::OrderType const order_type,
			genfile::ValueType const value_type
		) {
			// data is preallocated, just check it will fit.
			assert( n <= 10 ) ;
		}

		void set_value( std::size_t entry_i, genfile::MissingValue const value ) {
			// nothing to do
		}

		inline void set_value( std::size_t entry_i, Integer const value ) {
			m_ints[ m_sample_i ][ entry_i ] = value ;
		}

		void set_value( std::size_t entry_i, double const value ) {
			m_doubles[ m_sample_i ][ entry_i ] = value ;
		}

		void finalise() {
			// do nothing
		}
		
	private:
		std::size_t m_sample_i ;
		std::vector< uint32_t >& m_ploidy ;
		std::vector< std::vector< double > >& m_doubles ;
		std::vector< std::vector< int64_t > >& m_ints ;
		std::vector< std::vector< std::string > >& m_strings ;
	} ;
}

void DataReadTest::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Test / benchmark options" ) ;
	options[ "-read-test" ]
		.set_description( "Don't do anything; just read the data." ) ;
}

DataReadTest::DataReadTest() {}

void DataReadTest::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps_read = 0 ;
	m_ploidy.resize( number_of_samples ) ;
	m_doubles.resize( number_of_samples ) ;
	m_ints.resize( number_of_samples ) ;
	m_strings.resize( number_of_samples ) ;
	for( std::size_t i = 0; i < number_of_samples; ++i ) {
		m_doubles[i].resize( 10 ) ;
		m_ints[i].resize( 10 ) ;
		m_strings[i].resize( 10 ) ;
	}
}

void DataReadTest::processed_snp( genfile::VariantIdentifyingData const& snp, genfile::VariantDataReader::SharedPtr data_reader ) {
	std::vector< std::string > fields ;
	void(std::vector<std::string>::*push_back)(std::string const&) = &std::vector<std::string>::push_back ;
	data_reader->get_supported_specs( boost::bind( push_back, &fields, _1 )) ;

	DataSetter setter( m_ploidy, m_doubles, m_ints, m_strings ) ;
	for( std::size_t field_i = 0; field_i < fields.size(); ++field_i ) {
		std::string const field = fields[ field_i ] ;
		data_reader->get( field, setter ) ;
	}
}

void DataReadTest::end_processing_snps() {
	// nothing to do.
}

