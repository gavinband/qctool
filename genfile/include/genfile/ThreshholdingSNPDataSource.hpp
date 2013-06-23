
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_THRESHHOLDING_SNP_DATA_SOURCE_HPP
#define QCTOOL_THRESHHOLDING_SNP_DATA_SOURCE_HPP

#include <queue>
#include <memory>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/ptr_container/ptr_deque.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	// This SNPDataSource applies a threshhold to genotype calls.
	class ThreshholdingSNPDataSource: public SNPDataSource {
	public:
		ThreshholdingSNPDataSource( SNPDataSource::UniquePtr source, double const threshhold ) ;
		operator bool() const ;
		unsigned int number_of_samples() const ;
		OptionalSnpCount total_number_of_snps() const ;
		std::string get_source_spec() const ;
		SNPDataSource const& get_parent_source() const ;
		SNPDataSource const& get_base_source() const ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;
		void reset_to_start_impl() ;
		
	private:
		SNPDataSource::UniquePtr m_source ; // Accessed only by background thread.
		double const m_threshhold ;
	} ;
}

#endif
