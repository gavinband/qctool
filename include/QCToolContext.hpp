
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_CONTEXT_HPP
#define QCTOOL_CONTEXT_HPP

#include <vector>

#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSink.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "appcontext/OptionProcessor.hpp"

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
	virtual genfile::CohortIndividualSource const& samples() const = 0 ;
	virtual appcontext::OptionProcessor const& options() const = 0 ;
} ;

#endif
