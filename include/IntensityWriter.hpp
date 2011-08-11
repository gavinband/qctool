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

	std::map< std::string, db::Connection::RowId > m_entities ;

	std::vector< std::vector< genfile::VariantEntry > > m_data ;
	std::vector< char > m_buffer ;
	std::vector< char > m_compressed_buffer ;
	void setup() ;
	void get_or_create_entity( std::string const& name, std::string const& description ) ;
	void set_relationship( std::string const& left, std::string const& relation, std::string const& right ) const ;


} ;

#endif

