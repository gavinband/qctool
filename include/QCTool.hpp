#ifndef QCTOOL_PROCESSOR_HPP
#define QCTOOL_PROCESSOR_HPP

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <memory>
#include <numeric>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "genfile/SNPDataSourceProcessor.hpp"

#include "appcontext/UIContext.hpp"

#include "Timer.hpp"
#include "GenRow.hpp"
#include "SampleRow.hpp"
#include "AlleleProportions.hpp"
#include "GToolException.hpp"
#include "RowCondition.hpp"
#include "GenRowStatistics.hpp"
#include "SampleRowStatistics.hpp"
#include "ObjectSource.hpp"
#include "SimpleFileObjectSource.hpp"
#include "SimpleFileObjectSink.hpp"
#include "OstreamTee.hpp"
#include "QCToolContext.hpp"

struct QCToolException: public std::exception
{
	char const* what() const throw() {return "QCToolException" ; }
} ;

struct GenAndSampleFileMismatchException: public QCToolException
{
	GenAndSampleFileMismatchException( std::size_t actual_number_of_samples, std::size_t expected_number_of_samples )
		: m_expected_number_of_samples( expected_number_of_samples ),
			m_actual_number_of_samples( actual_number_of_samples )
	{}

	~GenAndSampleFileMismatchException() throw() {}
	
	char const* what() const throw() {return "GenAndSampleFileMismatchException" ; }

	std::size_t expected_number_of_samples() const { return m_expected_number_of_samples ; }
	std::size_t actual_number_of_samples() const { return m_actual_number_of_samples ; }
private:
	std::size_t m_expected_number_of_samples, m_actual_number_of_samples ;
} ;

struct QCTool: public genfile::SNPDataSourceProcessor::Callback
{
public:
	QCTool(
		QCToolContext & context,
		appcontext::UIContext& ui_context
	) ;
	
	void begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) ;
	void processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

private:
	QCToolContext& m_context ;
	appcontext::UIContext& m_ui_context ;
	std::size_t m_number_of_samples ;
	std::size_t m_number_of_snps ;
	std::size_t m_number_of_snps_processed ;
	std::size_t m_number_of_autosomal_snps_processed ;
	std::vector< GenotypeProportions > m_per_column_amounts ;
	Timer m_timer ;

	void unsafe_call_processed_snp(
		genfile::SNPIdentifyingData const& id_data,
		genfile::VariantDataReader& data_reader
	) ;
	void process_gen_row( GenRow const& row, std::size_t row_number ) ;
	void output_gen_row_stats( GenotypeAssayStatistics const& row_statistics ) ;
	void output_missing_gen_row_stats( GenRow const& row, GenotypeAssayStatistics const& row_statistics, std::size_t row_number ) ;
	void do_snp_filter_diagnostics( GenRowStatistics const& row_statistics, std::size_t const row_index ) ;
	void accumulate_per_column_amounts( GenRow& row, std::vector< GenotypeProportions >& per_column_amounts ) ;
	void process_sample_rows() ;
	void output_sample_stats( std::size_t index, GenotypeAssayStatistics const& stats ) ;
	void apply_sample_filter() ;
	bool sample_row_is_filtered_out( std::size_t const sample_row_index ) ;
	void do_sample_filter_diagnostics( SampleRow const& sample_row, std::size_t const sample_row_index ) ;
	void construct_plots() ;
} ;

#endif
