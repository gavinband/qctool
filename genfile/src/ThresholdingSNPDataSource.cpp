#include "genfile/SNPDataSource.hpp"
#include "genfile/ThresholdingSNPDataSource.hpp"

namespace genfile {
	ThresholdingSNPDataSource::ThresholdingSNPDataSource( SNPDataSource::UniquePtr source, double threshhold ) :
		m_source( source ),
		m_threshhold( threshhold )
	{
		assert( m_threshhold >= 0.5 ) ;
	}
	
	unsigned int ThresholdingSNPDataSource::number_of_samples() const { return m_source->number_of_samples() ; }
	unsigned int ThresholdingSNPDataSource::total_number_of_snps() const { return m_source->total_number_of_snps() ; }
	ThresholdingSNPDataSource::operator bool() const { return *m_source ; }
	void ThresholdingSNPDataSource::reset_to_start_impl() { m_source->reset_to_start() ; }
	void ThresholdingSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		m_source->get_snp_identifying_data(
			set_number_of_samples,
			set_SNPID,
			set_RSID,
			set_chromosome,
			set_SNP_position,
			set_allele1,
			set_allele2
		) ;
	}

	namespace impl {
		template< typename Setter >
		struct ThreshholdingSetter
		{
			ThreshholdingSetter( Setter const& setter, double threshhold ):
				m_setter( setter ),
				m_threshhold( threshhold )
			{}
			
			void operator()( std::size_t i, double AA, double AB, double BB ) {
				AA = ( AA >= m_threshhold ) ? 1.0 : 0.0 ;
				AB = ( AB >= m_threshhold ) ? 1.0 : 0.0 ;
				BB = ( BB >= m_threshhold ) ? 1.0 : 0.0 ;
				m_setter( i, AA, AB, BB ) ;
			}
		private:
			Setter const& m_setter ;
			double const m_threshhold ;
		} ;

		template< typename Setter >
		ThreshholdingSetter< Setter > make_threshholding_setter( Setter const& setter, double threshhold ) {
			return ThreshholdingSetter< Setter >( setter, threshhold ) ;
		}
	}

	void ThresholdingSNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		m_source->read_snp_probability_data(
			impl::make_threshholding_setter( set_genotype_probabilities, m_threshhold )
		) ;
	}
	
	void ThresholdingSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}
}
