#include <memory>
#include <vector>
#include "Eigen/Core"
#include "Eigen/Cholesky"
#include "genfile/VariantEntry.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/Error.hpp"
#include "genfile/vcf/get_set.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "ClusterFitter.hpp"
#include "NormalClusterFitter.hpp"

namespace {
	std::vector< std::pair< std::string, std::string > > parse_spec( std::string spec ) {
		std::vector< std::pair< std::string, std::string > > result ;
		std::vector< std::string > elts = genfile::string_utils::split_and_strip_discarding_empty_entries( spec, "," ) ;
		for( std::size_t i = 0; i < elts.size(); ++i ) {
			std::vector< std::string > these_elts = genfile::string_utils::split_and_strip( elts[i], "/" ) ;
			if( these_elts.size() != 2 ) {
				throw genfile::BadArgumentError( "(ClusterFitter.cpp) impl::parse_spec()", "spec=\"" + elts[i] + "\"" ) ;
			}
			result.push_back( std::make_pair( these_elts[0], these_elts[1] )) ;
		}
		return result ;
	}
}

NormalClusterFitter::NormalClusterFitter( appcontext::OptionProcessor const& options ):
	m_options( options ),
	m_spec( parse_spec( options.get_value< std::string >( "-fit-clusters" ))),
	m_call_threshhold( options.get_value< double >( "-call-threshhold" ))
{}

void NormalClusterFitter::begin_processing_snps( std::size_t number_of_samples ) {
	m_number_of_samples = number_of_samples ;
}

void NormalClusterFitter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
	Base::processed_snp( snp, data_reader ) ;

	genfile::SingleSNPGenotypeProbabilities genotypes ;
	Eigen::MatrixXd intensities ;

	Eigen::MatrixXd fit( 2, 9 ) ;
	std::vector< std::size_t > non_missing_counts( 3 ) ;

	for( std::size_t spec_i = 0; spec_i < m_spec.size(); ++spec_i ) {
		genfile::vcf::MatrixSetter< Eigen::MatrixXd > intensity_setter( intensities ) ;
		data_reader
			.get( m_spec[spec_i].first, genotypes )
			.get( m_spec[spec_i].second, intensity_setter )
		;
		IntensityModel::SharedPtr model( IntensityModel::estimate( intensities, genotypes ) ) ;
		send_results( snp, "NormalClusterFit:" + m_spec[ spec_i ].first + "/" + m_spec[ spec_i ].second, intensities, genotypes, model ) ;
	}
}

void NormalClusterFitter::end_processing_snps() {
	// nothing to do.
}
