#ifndef QCTOOL_DATA_READ_TEST_HPP
#define QCTOOL_DATA_READ_TEST_HPP

#include <boost/noncopyable.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "db/SQLite3Connection.hpp"

class DataReadTest: public genfile::SNPDataSourceProcessor::Callback, public boost::noncopyable {
public:
	static void declare_options( appcontext::OptionProcessor& options ) ;

	DataReadTest() ;
	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const& , genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;
private:
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps_read ;
	std::vector< std::vector< genfile::VariantEntry > > m_data ;
} ;

#endif

