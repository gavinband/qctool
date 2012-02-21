#include "genfile/vcf/get_set.hpp"
#include "components/CallComparerComponent/ConsensusCaller.hpp"

ConsensusCaller::ConsensusCaller():
	m_call_threshhold( 0.9 )
{}

void ConsensusCaller::begin_processing_snps( std::size_t number_of_samples ) {}

void ConsensusCaller::begin_comparisons( genfile::SNPIdentifyingData const& snp ) {
	m_call_names.clear() ;
}

void ConsensusCaller::set_result(
	std::string const& comparison,
	std::string const& comparison_value,
	genfile::VariantEntry const& value
) {
	if( comparison_value == "concordant_calls" ) {
		m_call_names = genfile::string_utils::split( value.as< string >(), "," ) ;
	}
}

void ConsensusCaller::end_comparisons() {}

void ConsensusCaller::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	m_calls.resize( m_call_names.size() ) ;
	std::vector< double > missingness( m_call_names.size() ) ;
	for( std::size_t i = 0; i < m_call_names.size(); ++i ) {
		{
			genfile::vcf::GenotypeSetter< Eigen::MatrixXd > setter( m_calls[i] ) ;
			data_reader.get( m_call_names[i], setter ) ;
		}
		double missing_calls = 0 ;
		for( int i = 0; i < m_calls[i].rows(); ++i ) {
			if( genotypes.row( i ).maxCoeff() < m_call_threshhold ) {
				++missing_calls ;
			}
		}
		missingness[i] = missing_calls ;
	}
	
	std::size_t which = ( std::min_element( missingness.begin(), missingness.end() ) - missingness.begin() ) ;
	m_sink->write_snp( snp, ... )
}

