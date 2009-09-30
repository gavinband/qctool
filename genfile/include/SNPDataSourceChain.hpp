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
		}

		unsigned int number_of_samples() const { return m_number_of_samples ; }
		unsigned int total_number_of_snps() const {
			unsigned int total_number_of_snps = 0 ;
			for( std::size_t i = 0; i < m_sources.size(); ++i ) {
				total_number_of_snps += m_sources[i]->total_number_of_snps() ;
			}
			return total_number_of_snps ;
		}

		unsigned int number_of_snps_in_source( std::size_t source_index ) {
			assert( source_index < m_sources.size() ) ;
			return m_sources[ source_index ]->total_number_of_snps() ;
		}

		operator bool() const {
			if( m_current_source < m_sources.size() ) {
				return *m_sources[ m_current_source ] ;
			}
			else {
				return false ;
			}
		}

	protected:

		void get_snp_identifying_data_impl( 
			IntegerSetter const& set_number_of_samples,
			StringSetter const& set_SNPID,
			StringSetter const& set_RSID,
			ChromosomeSetter const& set_chromosome,
			SNPPositionSetter const& set_SNP_position,
			AlleleSetter const& set_allele1,
			AlleleSetter const& set_allele2
		) {
			move_to_next_nonempty_source_if_necessary() ;
			if( m_current_source < m_sources.size() ) {
				m_sources[m_current_source]->get_snp_identifying_data( set_number_of_samples, set_SNPID, set_RSID, set_chromosome, set_SNP_position, set_allele1, set_allele2 ) ;
			}
		}

		void read_snp_probability_data_impl(
			uint32_t* number_of_samples,
			GenotypeProbabilitySetter const& set_genotype_probabilities
		) {
			assert( m_current_source < m_sources.size() ) ;
			m_sources[m_current_source]->read_snp_probability_data( number_of_samples, set_genotype_probabilities ) ;
		}

		void ignore_snp_probability_data_impl(
		) {
			assert( m_current_source < m_sources.size() ) ;
			m_sources[m_current_source]->ignore_snp_probability_data() ;
		}

	public:

	#if HAVE_BOOST_FUNCTION
		typedef boost::function< void( int index ) > moved_to_next_source_callback_t ;
	#else
		typedef void( *moved_to_next_source_callback_t )( std::size_t ) ;
	#endif

		void set_moved_to_next_source_callback( moved_to_next_source_callback_t callback ) { m_moved_to_next_source_callback = callback ; }

	private:

		void move_to_next_source() {
			++m_current_source ;
			if( m_moved_to_next_source_callback ) {
				m_moved_to_next_source_callback( m_current_source ) ;
			}
		}

		void move_to_next_nonempty_source_if_necessary() {
			while(( m_current_source < m_sources.size())
				&& (m_sources[ m_current_source ]->number_of_snps_read() >= m_sources[m_current_source]->total_number_of_snps())) {
				move_to_next_source() ;
			}
		}

		std::vector< SNPDataSource* > m_sources ;
		std::size_t m_current_source ;
		unsigned int m_number_of_samples ;
		
		moved_to_next_source_callback_t m_moved_to_next_source_callback ;
	} ;
}

#endif