#include <set>
#include <boost/bind.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SampleFilteringSNPDataSource.hpp"

namespace genfile {
	// Create a chain of SNPDataSources taking data from the specified files.
	std::auto_ptr< SampleFilteringSNPDataSource > SampleFilteringSNPDataSource::create(
		SNPDataSource& source,
		std::set< std::size_t > const& indices_of_samples_to_filter_out
	) {
		return std::auto_ptr< SampleFilteringSNPDataSource >( new SampleFilteringSNPDataSource( source, indices_of_samples_to_filter_out )) ;
	}
	
	SampleFilteringSNPDataSource::SampleFilteringSNPDataSource(
		SNPDataSource& source,
		std::set< std::size_t > const& indices_of_samples_to_filter_out
	)
		: m_source( source ),
		  m_indices_of_samples_to_filter_out( indices_of_samples_to_filter_out ),
		  m_genotype_data( m_source.number_of_samples() * 3 )
	{
		verify_indices( m_indices_of_samples_to_filter_out ) ;
	}
	
	void SampleFilteringSNPDataSource::verify_indices( std::set< std::size_t > const& indices ) const {
		std::set< std::size_t >::const_iterator
			i = indices.begin(),
			end_i = indices.end() ;
		for( ; i != end_i; ++i ) {
			if( *i >= m_source.number_of_samples() ) {
				throw SampleIndexOutOfRangeError( *i, m_source.number_of_samples() ) ;
			}
		}
	}
	
	SampleFilteringSNPDataSource::~SampleFilteringSNPDataSource() {}

	unsigned int SampleFilteringSNPDataSource::number_of_samples() const {
		return m_source.number_of_samples() - m_indices_of_samples_to_filter_out.size() ;
	}
	unsigned int SampleFilteringSNPDataSource::total_number_of_snps() const {
		return m_source.total_number_of_snps() ;
	}

	SampleFilteringSNPDataSource::operator bool() const {
		return m_source ;
	}

	void SampleFilteringSNPDataSource::reset_to_start_impl() {
		m_source.reset_to_start() ;
	}

	void SampleFilteringSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		uint32_t number_of_samples ;
		m_source.get_snp_identifying_data( 
			set_value( number_of_samples ),
			set_SNPID,
			set_RSID,
			set_chromosome,
			set_SNP_position,
			set_allele1,
			set_allele2
		) ;
		
		set_number_of_samples( m_source.number_of_samples() - m_indices_of_samples_to_filter_out.size() ) ;
	}

	void SampleFilteringSNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		read_source_probability_data() ;
		return_filtered_genotype_probabilities( set_genotype_probabilities ) ;
	}
	
	void SampleFilteringSNPDataSource::read_source_probability_data() {
		m_source.read_snp_probability_data(
		 	boost::bind( &SampleFilteringSNPDataSource::set_unfiltered_genotype_probabilities, this, _1, _2, _3, _4 )
		) ;
	}

	void SampleFilteringSNPDataSource::return_filtered_genotype_probabilities(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		std::set< std::size_t >::const_iterator
			i( m_indices_of_samples_to_filter_out.begin() ),
			end_i( m_indices_of_samples_to_filter_out.end() ) ;

		std::size_t index_of_last_unfiltered_sample = 0 ;
		std::size_t index_of_filtered_sample = 0 ;

		for( ; i != end_i; ++i ) {
			for( std::size_t sample_i = index_of_last_unfiltered_sample; sample_i < *i; ++sample_i ) {
				std::size_t genotype_data_index = sample_i * 3 ;
				set_genotype_probabilities(
					index_of_filtered_sample++,
					m_genotype_data[ genotype_data_index ],
					m_genotype_data[ genotype_data_index + 1 ],
					m_genotype_data[ genotype_data_index + 2 ]
				) ;
			}
			index_of_last_unfiltered_sample = (*i) + 1 ;
		}

		for( std::size_t sample_i = index_of_last_unfiltered_sample; sample_i < m_source.number_of_samples() ; ++sample_i ) {
			std::size_t genotype_data_index = sample_i * 3 ;
			set_genotype_probabilities(
				index_of_filtered_sample++,
				m_genotype_data[ genotype_data_index ],
				m_genotype_data[ genotype_data_index + 1 ],
				m_genotype_data[ genotype_data_index + 2 ]
			) ;
		}
	}
	
	void SampleFilteringSNPDataSource::set_unfiltered_genotype_probabilities( std::size_t i, double aa, double ab, double bb ) {
		i *= 3 ;
		assert( i + 2 < m_genotype_data.size() ) ;
		m_genotype_data[i++] = aa ;
		m_genotype_data[i++] = ab ;
		m_genotype_data[i] = bb ;
	}

	void SampleFilteringSNPDataSource::ignore_snp_probability_data_impl() {
		m_source.ignore_snp_probability_data() ;
	}

}

