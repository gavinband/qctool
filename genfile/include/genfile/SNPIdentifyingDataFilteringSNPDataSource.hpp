
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GENFILE_SNPIDENTIFYINGDATAFILTERINGSNPDATASOURCE_HPP
#define GENFILE_SNPIDENTIFYINGDATAFILTERINGSNPDATASOURCE_HPP

#include <vector>
#include <set>
#include <memory>
#include <boost/signals2/signal.hpp>
#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingDataTest.hpp"
#include "genfile/CommonSNPFilter.hpp"

namespace genfile {
	// class SNPIdentifyingDataFilteringSNPDataSource
	// This class represents a view of its parent SNPDataSource such that
	// only SNPs whose chromosome, position, id fields and alleles satisfy the
	// given test are visible.
	class SNPIdentifyingDataFilteringSNPDataSource: public SNPDataSource
	{
	public:
		typedef std::auto_ptr< SNPIdentifyingDataFilteringSNPDataSource > UniquePtr ;
		// Create a SNPIdentifyingDataFilteringSNPDataSource from the given source and the given test.
		// The test should return true for tests that are in and false for those excluded.
		static UniquePtr create( SNPDataSource::UniquePtr source, SNPIdentifyingDataTest::UniquePtr test ) ;
		static UniquePtr create( SNPDataSource::UniquePtr source, SNPIdentifyingDataTest const& test ) ;
		
	public:
		~SNPIdentifyingDataFilteringSNPDataSource() ;
		SNPIdentifyingDataFilteringSNPDataSource( SNPDataSource::UniquePtr source, SNPIdentifyingDataTest::UniquePtr test ) ;
		SNPIdentifyingDataFilteringSNPDataSource( SNPDataSource::UniquePtr source, SNPIdentifyingDataTest const& test ) ;

		operator bool() const ;
		Metadata get_metadata() const ;
		unsigned int number_of_samples() const ;
		void get_sample_ids( GetSampleIds ) const ;
		OptionalSnpCount total_number_of_snps() const ;
		std::string get_source_spec() const ;

		SNPDataSource const& get_parent_source() const ;
		SNPDataSource const& get_base_source() const ;
		
		typedef boost::signals2::signal< void ( SNPIdentifyingData const& ) > SNPSignal ;
		void send_filtered_out_SNPs_to( SNPSignal::slot_type ) ;

	private:
		void reset_to_start_impl() ;

		void get_snp_identifying_data_impl( 
VariantIdentifyingData* variant		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;

	protected:
		SNPDataSource& source() { return *m_source ; }
		
	private:
		
		SNPDataSource::UniquePtr m_source ;
		bool const m_manage_test ;
		SNPIdentifyingDataTest const* m_test ;
		SNPSignal m_filtered_out_snp_signal ;
	} ;
}

#endif
