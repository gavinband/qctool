
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_CLUSTER_PLOTTER_HPP
#define QCTOOL_CLUSTER_PLOTTER_HPP

#include <memory>
#include <string>
#include <boost/ptr_container/ptr_vector.hpp>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "worker/Worker.hpp"

struct ClusterPlotter: public genfile::SNPDataSourceProcessor::Callback
{
public:
	typedef std::auto_ptr< ClusterPlotter > UniquePtr ;
	static void declare_options( appcontext::OptionProcessor& options ) ;
	static UniquePtr create( appcontext::OptionProcessor const& options, worker::Worker* worker ) ;

public:
	ClusterPlotter(
		std::string const& filename_template,
		std::vector< std::string > const& call_fields,
		worker::Worker* worker
	) ;
	void begin_processing_snps( std::size_t number_of_samples ) ;
	void processed_snp( genfile::SNPIdentifyingData const&, genfile::VariantDataReader& data_reader ) ;
	void end_processing_snps() ;

private:
	std::string const m_filename_template ;
	std::vector< std::string > const m_call_fields ;
	std::string const m_intensity_field ;
	std::size_t m_number_of_samples ;
	double m_call_threshhold ;
	worker::Worker* m_worker ;
	boost::ptr_vector< worker::Task > m_tasks ;
	std::size_t const m_max_tasks ;
} ;

#endif

