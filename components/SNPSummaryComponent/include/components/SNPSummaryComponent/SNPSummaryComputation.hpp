
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SNP_SUMMARY_COMPUTATION_HPP
#define QCTOOL_SNP_SUMMARY_COMPUTATION_HPP

#include <string>
#include <memory>
#include <boost/noncopyable.hpp>
#include <boost/function.hpp>
#include <Eigen/Core>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantDataReader.hpp"
#include "appcontext/OptionProcessor.hpp"

struct SNPSummaryComputation: public boost::noncopyable {
	typedef std::auto_ptr< SNPSummaryComputation > UniquePtr ;
	virtual ~SNPSummaryComputation() {}
	static UniquePtr create( std::string const& name ) ;

	typedef genfile::VariantIdentifyingData VariantIdentifyingData ;
	typedef Eigen::MatrixXd Genotypes ;
	typedef boost::function< void ( std::string const& value_name, genfile::VariantEntry const& value ) > ResultCallback ;
	typedef boost::function< void ( std::size_t sample_i, std::string const& value_name, genfile::VariantEntry const& value ) > PerSampleResultCallback ;
	typedef std::vector< char > SampleSexes ;

	virtual std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const = 0 ;

	virtual void begin_processing_snps( std::size_t ) {}
	virtual void operator()( VariantIdentifyingData const&, Genotypes const&, SampleSexes const&, genfile::VariantDataReader&, ResultCallback ) = 0 ;
	virtual void end_processing_snps( PerSampleResultCallback ) {}
} ;

#endif
