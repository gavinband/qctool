#include "genfile/vcf/get_set.hpp"
#include "genfile/get_set_eigen.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/CallComparerComponent/ConsensusCaller.hpp"

ConsensusCaller::ConsensusCaller( genfile::SNPDataSink::UniquePtr sink ):
	m_call_threshhold( 0.9 ),
	m_sink( sink )
{}

void ConsensusCaller::begin_processing_snps( std::size_t number_of_samples ) {
	m_number_of_samples = number_of_samples ;
}

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

namespace impl {
	
}

void ConsensusCaller::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	m_genotypes.resize( m_call_names.size() ) ;

	std::map< std::string, std::vector< genfile::VariantEntry > > info ;
	std::size_t chosen_call = 0 ;
	double chosen_call_missingness = std::numeric_limits< double >::max() ;
	if( m_call_names.size() > 0 ) {
		for( std::size_t i = 0; i < m_call_names.size(); ++i ) {
			{
				genfile::vcf::GenotypeSetter< Eigen::MatrixBase< Eigen::MatrixXd > > setter( m_genotypes[i] ) ;
				data_reader.get( m_call_names[i], setter ) ;
			}
			double missing_calls = 0 ;
			for( int j = 0; j < m_genotypes[i].rows(); ++j ) {
				if( m_genotypes[i].row( j ).maxCoeff() < m_call_threshhold ) {
					++missing_calls ;
				}
			}
			info[ m_call_names[i] ].push_back( missing_calls / m_genotypes[i].rows() ) ;
			if( missing_calls < chosen_call_missingness ) {
				chosen_call = i ;
				chosen_call_missingness = missing_calls ;
			}
		}
	
		info[ "call" ].push_back( m_call_names[ chosen_call ] ) ;

		m_sink->write_snp(
			m_genotypes[ chosen_call ].rows(),
			snp,
			genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes[ chosen_call ], 0 ),
			genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes[ chosen_call ], 1 ),
			genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes[ chosen_call ], 2 ),
			info
		) ;
	} else {
		m_genotypes.resize( 1 ) ;
		m_genotypes[0].setZero( m_number_of_samples, 3 ) ;
		info[ "call" ].push_back( genfile::MissingValue() ) ;
		m_sink->write_snp(
			m_genotypes[ 0 ].rows(),
			snp,
			genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes[ 0 ], 0 ),
			genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes[ 0 ], 1 ),
			genfile::GenotypeGetter< Eigen::MatrixXd >( m_genotypes[ 0 ], 2 ),
			info
		) ;
	}
}

