
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_SAMPLESUMMARYCOMPONENT_FLAT_FILE_OUTPUTTER_HPP
#define QCTOOL_SAMPLESUMMARYCOMPONENT_FLAT_FILE_OUTPUTTER_HPP

#include <string>
#include <memory>
#include <map>
#include <vector>
#include <utility>
#include <boost/bimap.hpp>
#include "genfile/VariantEntry.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "components/SampleSummaryComponent/SampleStorage.hpp"

namespace sample_stats {
	struct FlatFileOutputter: public SampleStorage {
		typedef std::map< std::string, std::pair< std::vector< std::string >, std::string > > Metadata ;
		static UniquePtr create(
			genfile::CohortIndividualSource const& samples,
			std::string const& filename,
			std::string const& analysis_name,
			Metadata const& metadata 
		) ;
		static SharedPtr create_shared(
			genfile::CohortIndividualSource const& samples,
			std::string const& filename,
			std::string const& analysis_name,
			Metadata const& metadata
		) ;

		FlatFileOutputter(
			genfile::CohortIndividualSource const& samples,
			std::string const& filename,
			std::string const& analysis_name,
			Metadata const& metadata
		) ;
		~FlatFileOutputter() ;

		void add_variable( std::string const& ) ;

		void store_per_sample_data(
			std::string const& computation_name,
			std::size_t sample,
			std::string const& variable,
			std::string const& description,
			genfile::VariantEntry const& value
		) ;

		void finalise( long options = qcdb::eCreateIndices ) ;

		AnalysisId analysis_id() const ;
	private:
		std::vector< genfile::VariantEntry > m_samples ;
		std::string const m_filename ;
		std::string const m_analysis_name ;
		Metadata const m_metadata ;
		std::size_t const m_max_samples_per_block ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
		std::vector< std::size_t > m_sample_indices ;
		typedef boost::bimap< std::string, std::size_t > VariableMap ;
		VariableMap m_variables ;
		typedef std::map< std::pair< std::size_t, std::size_t >, genfile::VariantEntry > ValueMap ;
		ValueMap m_values ;
	
	private:
		void store_block() ;
		std::string format_metadata() const ;
	} ;
}

#endif
