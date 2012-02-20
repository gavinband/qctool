//#include <boost/signal.hpp>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	SNPDataSourceProcessor::Callback::~Callback() {
	}
	
	SNPDataSourceProcessor::~SNPDataSourceProcessor() {
		for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
			if( m_should_manage[i] ) {
				delete m_callbacks[i] ;
			}
		}
	}

	void SNPDataSourceProcessor::add_callback( Callback::UniquePtr callback ) {
		m_should_manage.push_back( true ) ;
		m_callbacks.push_back( callback.release() ) ;
	}

	void SNPDataSourceProcessor::add_callback( Callback& callback ) {
		m_should_manage.push_back( false ) ;
		m_callbacks.push_back( &callback ) ;
	}

	void SNPDataSourceProcessor::call_begin_processing_snps( std::size_t const& number_of_samples ) const {
		for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
			m_callbacks[i]->begin_processing_snps( number_of_samples ) ;
		}
	}
	
	void SNPDataSourceProcessor::call_processed_snp(  SNPIdentifyingData const& id_data, VariantDataReader& data_reader ) const {
		for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
			m_callbacks[i]->processed_snp( id_data, data_reader ) ;
		}
	}

	void SNPDataSourceProcessor::call_end_processing_snps() const {
		for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
			m_callbacks[i]->end_processing_snps() ;
		}
	}

	void SimpleSNPDataSourceProcessor::process( genfile::SNPDataSource& source, ProgressCallback progress_callback ) {
		SNPIdentifyingData id_data ;
		SingleSNPGenotypeProbabilities genotypes( source.number_of_samples() ) ;
		
		call_begin_processing_snps( source.number_of_samples() ) ;

		while( source.get_snp_identifying_data(
				ignore(),
				set_value( id_data.SNPID() ),
				set_value( id_data.rsid() ),
				set_value( id_data.position().chromosome() ),
				set_value( id_data.position().position() ),
				set_value( id_data.first_allele() ),
				set_value( id_data.second_allele() )
			)
		) {
			VariantDataReader::UniquePtr data_reader = source.read_variant_data() ;
			call_processed_snp( id_data, *data_reader ) ;
			if( progress_callback ) {
				progress_callback( source.number_of_snps_read(), source.total_number_of_snps() ) ;
			}
		}

		call_end_processing_snps() ;
	}
}

