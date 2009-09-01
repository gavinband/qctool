#ifndef QCTOOL_CONTEXT_HPP
#define QCTOOL_CONTEXT_HPP

#include <vector>

#include "SNPDataSource.hpp"
#include "SNPDataSink.hpp"
#include "ObjectSource.hpp"
#include "Timer.hpp"
#include "GenRow.hpp"
#include "GenRowStatistics.hpp"
#include "SampleRow.hpp"
#include "SampleRowStatistics.hpp"
#include "GToolException.hpp"
#include "RowCondition.hpp"
#include "SimpleFileObjectSource.hpp"
#include "SimpleFileObjectSink.hpp"
#include "OstreamTee.hpp"
#include "statfile/StatSink.hpp"

struct QCToolContext
{
	typedef genfile::SNPDataSource SNPDataSource ;
	typedef genfile::SNPDataSink SNPDataSink ;
	
	virtual ~QCToolContext() {}
	
	virtual SNPDataSource& snp_data_source() const = 0 ;
	virtual SNPDataSink& fltrd_in_snp_data_sink() const = 0 ;
	virtual SNPDataSink& fltrd_out_snp_data_sink() const = 0 ;
	virtual ObjectSink< SampleRow >& fltrd_in_sample_sink() const = 0 ;
	virtual ObjectSink< SampleRow >& fltrd_out_sample_sink() const = 0 ;
	virtual statfile::BuiltInTypeStatSink& snp_stats_sink() const = 0 ;
	virtual statfile::BuiltInTypeStatSink& sample_stats_sink() const = 0 ;
	virtual AndRowCondition& snp_filter() const = 0 ;
	virtual AndRowCondition& sample_filter() const = 0 ;
	virtual std::vector< std::size_t >& snp_filter_failure_counts() = 0 ;
	virtual std::vector< std::size_t >& sample_filter_failure_counts() = 0 ;
	virtual GenRowStatistics& snp_statistics() = 0 ;
	virtual SampleRowStatistics& sample_statistics() = 0 ;
	virtual std::vector< SampleRow >& sample_rows() = 0 ;
	virtual std::vector< std::size_t > const& indices_of_filtered_out_samples() const = 0 ;
	virtual void print_progress_if_necessary() = 0 ;
	virtual OstreamTee& logger() = 0 ;
} ;

#endif
