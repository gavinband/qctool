#include "compute_maximum_likelihood_allele_frequency.hpp"
#include "KinshipCoefficientComputation.hpp"

KinshipCoefficientComputation::KinshipCoefficientComputation(
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
):
	m_threshhold( 0.9 )
{}

void KinshipCoefficientComputation::prepare(
	std::vector< genfile::SNPIdentifyingData > const& snps,
	std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
) {
	m_allele_frequencies.resize( snps.size() ) ;
	for( std::size_t i = 0; i < snps.size(); ++i ) {
		double allele_sum = 0.0 ;
		double N = 0.0 ;
		for( std::size_t sample_i = 0; sample_i < genotypes[i].get_number_of_samples(); ++sample_i ) {
			for( std::size_t g = 0; g < 3; ++g ) {
				if( genotypes[ i ]( sample_i, g ) > m_threshhold ) {
					++N ;
					allele_sum += g ;
					break ;
				}
			}
		}
		m_allele_frequencies[i] = allele_sum / ( 2 * N ) ;
	}
}

double KinshipCoefficientComputation::operator()(
	std::size_t const sample1,
	std::size_t const sample2,
	std::vector< genfile::SingleSNPGenotypeProbabilities > const& genotypes
) {
	double non_missing_count = 0.0 ;
	double result = 0.0 ;
	for( std::size_t snp_i = 0; snp_i < genotypes.size(); ++snp_i ) {
		if( m_allele_frequencies[ snp_i] > 0.0 && m_allele_frequencies[ snp_i] < 1.0 ) {
			std::size_t
				genotype1 = std::numeric_limits< std::size_t >::max(),
				genotype2 = std::numeric_limits< std::size_t >::max()
			;

			for( std::size_t g = 0; g < 3; ++g ) {
				if( genotypes[snp_i]( sample1, g ) >= m_threshhold ) {
					genotype1 = g ;
					break ;
				}
			}

			for( std::size_t g = 0; g < 3; ++g ) {
				if( genotypes[snp_i]( sample2, g ) >= m_threshhold ) {
					genotype2 = g ;
					break ;
				}
			}

			if( genotype1 < 3 && genotype2 < 3 ) {
				++non_missing_count ;
				if( sample1 != sample2 ) {
					// Relationship between distinct individuals.
					double twice_allele_frequency = 2.0 * m_allele_frequencies[ snp_i ] ;
					result += ( genotype1 - twice_allele_frequency ) * ( genotype2 - twice_allele_frequency )
						/ ( 2.0 * twice_allele_frequency * ( 1.0 - m_allele_frequencies[ snp_i ] ) ) ;
				}
				else {
					// Relationship of individual with itself.
					double Fhat ;
					switch( genotype1 ) {
						case 0:
							Fhat = m_allele_frequencies[ snp_i ] / ( 1.0 - m_allele_frequencies[ snp_i ] ) ;
						break ;
						case 1:
							Fhat = -1 ;
						break ;
						case 2:
							Fhat = (1.0 - m_allele_frequencies[ snp_i ] ) / ( m_allele_frequencies[ snp_i ] ) ;
						break ;
						default:
							assert(0) ;
					}
					result += 1.0 + Fhat ;
				}
			}
		}
	}
	return result / non_missing_count ;
}
	
std::string KinshipCoefficientComputation::get_summary() const {
	return "Kinship coefficients (\\hat{A} as in Powell, Visscher and Goddard, NRG 2010) calculated using SNPs well-called in both samples." ;
}
