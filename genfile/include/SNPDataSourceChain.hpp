#ifndef SNPDATAsourceCHAIN_HPP
#define SNPDATAsourceCHAIN_HPP

#include <iostream>
#include <string>
#if HAVE_BOOST_FUNCTION
#include <boost/function.hpp>
#endif
#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"

namespace genfile {
	// class SNPDataSourceChain represnets a SNPDataSource
	// which gets it data sequentially from a collection of other SNPDataSources
	class SNPDataSourceChain: public SNPDataSource
	{
	public:
		SNPDataSourceChain(): m_current_source(0), m_number_of_samples(0), m_moved_to_next_source_callback(0) {}

		~SNPDataSourceChain() {
			for( std::size_t i = 0; i < m_sources.size(); ++i ) {
				delete m_sources[i] ;
			}
		} ;

		void add_source( std::auto_ptr< SNPDataSource > source ) {
			if( m_sources.empty() ) {
				m_number_of_samples = source->number_of_samples() ;
				m_current_source = 0 ;
			}
			else if( source->number_of_samples() != m_number_of_samples ) {
				throw FileContainsSNPsOfDifferentSizes() ;
			}
			m_sources.push_back( source.release() ) ;
			m_number_of_snps_read.push_back( 0u ) ;
		}

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const {
			unsigned int total_number_of_snps = 0 ;
			for( std::size_t i = 0; i < m_sources.size(); ++i ) {
				total_number_of_snps += m_sources[i]->total_number_of_snps() ;
			}
			return total_number_of_snps ;
		}

		FormatType format() const {
			assert( m_current_source < m_sources.size() ) ;
			return m_sources[ m_current_source ]->format() ;
		}

		operator bool() const {
			if( m_current_source < m_sources.size() ) {
				bool result = static_cast< bool >( *m_sources[ m_current_source ] ) ;
				return result ;
			}
			else {
				return false ;
			}
		}

		std::istream const& stream() const {
			assert( m_current_source < m_sources.size() ) ;
			return m_sources[ m_current_source ]->stream() ;
		}

		std::istream& stream() {
			assert( m_current_source < m_sources.size() ) ;
			return m_sources[ m_current_source ]->stream() ;
		}

	#if HAVE_BOOST_FUNCTION
		typedef boost::function< void( int index ) > moved_to_next_source_callback_t ;
	#else
		typedef void( *moved_to_next_source_callback_t )( std::size_t ) ;
	#endif

		void set_moved_to_next_source_callback( moved_to_next_source_callback_t callback ) { m_moved_to_next_source_callback = callback ; }

	private:

		void pre_read_snp() {
			// Make sure we switch to the next source when necessary.
			while(( m_current_source < m_sources.size())
				&& (m_number_of_snps_read[ m_current_source ] >= m_sources[m_current_source]->total_number_of_snps())) {
				move_to_next_source() ;
			}
		}
		
		void move_to_next_source() {
			++m_current_source ;
			if( m_moved_to_next_source_callback ) {
				m_moved_to_next_source_callback( m_current_source ) ;
			}
		}

		void post_read_snp() {
			++m_number_of_snps_read[ m_current_source ] ;
		}
	
		std::vector< SNPDataSource* > m_sources ;
		std::vector< std::size_t > m_number_of_snps_read ;
		std::size_t m_current_source ;
		unsigned int m_number_of_samples ;
		
		moved_to_next_source_callback_t m_moved_to_next_source_callback ;
	} ;
}

#endif