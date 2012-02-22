#include "genfile/vcf/get_set.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/CallComparerComponent/ConsensusCaller.hpp"

ConsensusCaller::ConsensusCaller( genfile::SNPDataSink::UniquePtr sink ):
	m_call_threshhold( 0.9 ),
	m_sink( sink )
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
		m_call_names = genfile::string_utils::split( value.as< std::string >(), "," ) ;
	}
}

void ConsensusCaller::end_comparisons() {}

void ConsensusCaller::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	m_genotypes.resize( m_call_names.size() ) ;
	std::vector< double > missingness( m_call_names.size() ) ;
	for( std::size_t i = 0; i < m_call_names.size(); ++i ) {
		{
			genfile::vcf::GenotypeSetter< Eigen::MatrixBase< Eigen::MatrixXd > > setter( m_genotypes[i] ) ;
			data_reader.get( m_call_names[i], setter ) ;
		}
		double missing_calls = 0 ;
		for( int j = 0; j < m_genotypes[i].rows(); ++i ) {
			if( m_genotypes[i].row( j ).maxCoeff() < m_call_threshhold ) {
				++missing_calls ;
			}
		}
		missingness[i] = missing_calls ;
	}
	
/*
	std::size_t which = ( std::min_element( missingness.begin(), missingness.end() ) - missingness.begin() ) ;
	genfile::GenotypeGetter< Eigen::MatrixXd > getterAA( m_genotypes[ which ], 0 );
	m_sink->write_snp(
		snp,
		genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes[ which ], 0 ),
		genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes[ which ], 1 ),
		genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes[ which ], 2 )
	) ;
	*/
}

