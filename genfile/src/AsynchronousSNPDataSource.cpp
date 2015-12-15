
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)


#include <boost/thread.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/AsynchronousSNPDataSource.hpp"

// #define DEBUG_ASYNCHRONOUS_SNP_DATA_SOURCE 1

namespace genfile {
	AsynchronousSNPDataSource::AsynchronousSNPDataSource( SNPDataSource::UniquePtr source ):
		m_source( source ),
		m_max_queue_size( 100 ),
		m_source_empty( false ),
		m_reset( false ),
		m_good( true )
	{
		assert( m_source.get() ) ;
		m_data_thread.reset(
			new boost::thread(
				boost::bind(
					&AsynchronousSNPDataSource::data_thread_read_data,
					this
				)
			)
		) ;
	}

	AsynchronousSNPDataSource::operator bool() const {
		return m_good ;
	}

	SNPDataSource::Metadata AsynchronousSNPDataSource::get_metadata() const {
		return m_source->get_metadata() ;
	}
	
	unsigned int AsynchronousSNPDataSource::number_of_samples() const {
		return m_source->number_of_samples() ;
	}
	SNPDataSource::OptionalSnpCount AsynchronousSNPDataSource::total_number_of_snps() const {
		return m_source->total_number_of_snps() ;
	}
	std::string AsynchronousSNPDataSource::get_source_spec() const {
		return "AsynchronousSNPDataSource( " + m_source->get_source_spec() + " )" ;
	}
	SNPDataSource const& AsynchronousSNPDataSource::get_parent_source() const {
		return *m_source ;
	}
	SNPDataSource const& AsynchronousSNPDataSource::get_base_source() const {
		return m_source->get_base_source() ;
	}
	
	std::string AsynchronousSNPDataSource::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + "AsynchronousSNPDataSource( " + m_source->get_summary() + ")" ;
	}
	
	namespace {
		
	}
	
	void AsynchronousSNPDataSource::get_snp_identifying_data_impl( 
		VariantIdentifyingData* variant
	) {
		VariantIdentifyingData snp ;
		{
			ScopedLock lock( m_mutex ) ;
			while( m_snp_queue.empty() && !m_source_empty ) {
#if DEBUG_ASYNCHRONOUS_SNP_DATA_SOURCE
				std::cerr << "r" ;
#endif
				m_queue_not_empty.wait( lock ) ;
			}

			if( m_source_empty ) {
				m_good = false ;
				return ;
			}

			snp = m_snp_queue.front() ;
		}
		*variant = snp ;
#if DEBUG_ASYNCHRONOUS_SNP_DATA_SOURCE
		std::cerr << "i" ;
#endif
	}

	VariantDataReader::UniquePtr AsynchronousSNPDataSource::read_variant_data_impl() {
		VariantDataReader::UniquePtr result ;
		{
			ScopedLock lock( m_mutex ) ;
			assert( !m_data_queue.empty() ) ;
			result.reset( m_data_queue.pop_front().release() ) ;
			m_snp_queue.pop() ;
			if( m_snp_queue.size() == m_max_queue_size - 1 ) {
				m_queue_not_full.notify_one() ;
			}
		}
#if DEBUG_ASYNCHRONOUS_SNP_DATA_SOURCE
		std::cerr << "-" ;
#endif		
		return result ;
	}

	void AsynchronousSNPDataSource::ignore_snp_probability_data_impl() {
		ScopedLock lock( m_mutex ) ;
		assert( !m_data_queue.empty() ) ;
		m_snp_queue.pop() ;
		m_data_queue.pop_front() ;
		if( m_snp_queue.size() == m_max_queue_size - 1 ) {
			m_queue_not_full.notify_one() ;
		}
#if DEBUG_ASYNCHRONOUS_SNP_DATA_SOURCE
		std::cerr << "-" ;
#endif
	}

	void AsynchronousSNPDataSource::reset_to_start_impl() {
		ScopedLock lock( m_mutex ) ;
		while( !m_snp_queue.empty() ) {
			m_snp_queue.pop() ;
		}
		m_data_queue.clear() ;
		m_reset = true ;
		m_good = true ;
		m_queue_not_full.notify_one() ;
	}
	
	void AsynchronousSNPDataSource::data_thread_read_data() {
		VariantIdentifyingData snp ;
		
		while( m_source->get_snp_identifying_data( &snp ) ) {
			VariantDataReader::UniquePtr variant_data = m_source->read_variant_data() ;
			
			ScopedLock lock( m_mutex ) ;
			
			if( m_reset ) {
				m_source->reset_to_start() ;
				m_reset = false ;
			}
			else {
				// wait for space in the queue.
				while( m_snp_queue.size() == m_max_queue_size ) {
#if DEBUG_ASYNCHRONOUS_SNP_DATA_SOURCE
					std::cerr << "w" ;
#endif
					m_queue_not_full.wait( lock ) ;
				}
				
#if DEBUG_ASYNCHRONOUS_SNP_DATA_SOURCE
				std::cerr << "+" ;
#endif			
				m_snp_queue.push( snp ) ;
				m_data_queue.push_back( variant_data ) ;
			
				// wake up the other thread if necessary
				// It will wait if there was nothing in the queue.
				if( m_snp_queue.size() == 1 ) {
					m_queue_not_empty.notify_one() ;
				}
			}
		}

		ScopedLock lock( m_mutex ) ;
		m_source_empty = true ;
		m_queue_not_empty.notify_one() ;
	}
}
