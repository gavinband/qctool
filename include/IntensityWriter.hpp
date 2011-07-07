#ifndef QCTOOL_INTENSITY_WRITER_HPP
#define QCTOOL_INTENSITY_WRITER_HPP

#include <boost/noncopyable.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "db/SQLite3Connection.hpp"

class IntensityWriter: public genfile::SNPDataSourceProcessor::Callback, public boost::noncopyable {
public:
	IntensityWriter( std::string const& filename ) ;
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( genfile::SNPIdentifyingData const& , genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;
private:
	std::string const m_filename ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps ;
	std::size_t m_number_of_snps_written ;
	db::Connection::UniquePtr m_connection ;

	void setup( db::Connection& connection ) ;
} ;

#endif

