
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef QCTOOL_ASYNCHRONOUS_SNP_DATA_SOURCE_HPP
#define QCTOOL_ASYNCHRONOUS_SNP_DATA_SOURCE_HPP

#include <queue>
#include <memory>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/ptr_container/ptr_deque.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	// This SNPDataSource asynchronously and eagerly reads data from the source given to it
	// in a background thread.  The intent is to keep data throughput at its maximum.
	class AsynchronousSNPDataSource: public SNPDataSource {
	public:
		AsynchronousSNPDataSource( SNPDataSource::UniquePtr source ) ;
		operator bool() const ;
		Metadata get_metadata() const ;
		unsigned int number_of_samples() const ;
		OptionalSnpCount total_number_of_snps() const ;
		std::string get_source_spec() const ;
		SNPDataSource const& get_parent_source() const ;
		SNPDataSource const& get_base_source() const ;
		std::string get_summary( std::string const& prefix = "", std::size_t column_width = 20 ) const ;
		void get_snp_identifying_data_impl( 
			VariantIdentifyingData* variant
		) ;

		VariantDataReader::UniquePtr read_variant_data_impl() ;

		void ignore_snp_probability_data_impl() ;
		void reset_to_start_impl() ;
		
	private:
		SNPDataSource::UniquePtr m_source ; // Accessed only by background thread.
		std::size_t const m_max_queue_size ;
		
		// This mutex and condition variables control access to the queues and the flag variables.
		typedef boost::mutex Mutex ;
		typedef boost::condition Condition ;
		typedef Mutex::scoped_lock ScopedLock ;
		Mutex m_mutex ;
		Condition m_queue_not_empty ;
		Condition m_queue_not_full ;

		// The following are accessed by both background and main threads and are protected by the above mutex.
		std::queue< VariantIdentifyingData > m_snp_queue ;
		boost::ptr_deque< VariantDataReader > m_data_queue ;
		bool m_source_empty ;
		bool m_reset ;

		// This is accessed by foreground thread only
		bool m_good ;
		
		// Background thread for data.
		std::auto_ptr< boost::thread > m_data_thread ;
		
	private:
		// This function runs in the background thread.
		void data_thread_read_data() ;
	} ;
}

#endif
