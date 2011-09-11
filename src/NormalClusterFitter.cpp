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
	m_filename( options.get_value< std::string >( "-fit-cluster-file" )),
	m_file( genfile::open_text_file_for_output( m_filename )),
	m_spec( parse_spec( options.get_value< std::string >( "-fit-clusters" ))),
	m_call_threshhold( options.get_value< double >( "-call-threshhold" ))
{}

void NormalClusterFitter::begin_processing_snps( std::size_t number_of_samples, std::size_t number_of_snps ) {
	m_number_of_samples = number_of_samples ;
	m_number_of_snps = number_of_snps ;
	write_output_columns( m_spec ) ;
}

void NormalClusterFitter::processed_snp( genfile::SNPIdentifyingData const& id_data, genfile::VariantDataReader& data_reader ) {
	std::vector< genfile::VariantEntry > genotypes ;
	Eigen::MatrixXd intensities ;

	std::vector< Eigen::Vector2d > means( 3, Eigen::Vector2d( 2 ) ) ;
	std::vector< Eigen::Matrix2d > variances( 3, Eigen::Matrix2d( 2, 2 ) ) ;
	std::vector< std::size_t > non_missing_counts( 3 ) ;

	write_id_data( id_data ) ;

	for( std::size_t spec_i = 0; spec_i < m_spec.size(); ++spec_i ) {
		for( std::size_t g = 0; g < 3; ++g ) {
			means[g].setZero() ;
			variances[g].setZero() ;
			non_missing_counts[g] = 0 ;
		}
		genfile::vcf::GenotypeSetter< std::vector< genfile::VariantEntry > > genotype_setter( genotypes, m_call_threshhold ) ;
		genfile::vcf::MatrixSetter< Eigen::MatrixXd > intensity_setter( intensities ) ;
		data_reader
			.get( m_spec[spec_i].first, genotype_setter )
			.get( m_spec[spec_i].second, intensity_setter )
		;
		get_cluster_fit( genotypes, intensities, means, variances, non_missing_counts ) ;
		write_cluster_fit( means, variances, non_missing_counts ) ;
	}
	(*m_file) << "\n" ;
}

void NormalClusterFitter::get_cluster_fit(
	std::vector< genfile::VariantEntry > const& genotypes,
	Eigen::MatrixXd const& intensities,
	std::vector< Eigen::Vector2d >& means,
	std::vector< Eigen::Matrix2d >& variances,
	std::vector< std::size_t >& non_missing_counts
) {
	assert( genotypes.size() == m_number_of_samples ) ;
	assert( std::size_t( intensities.cols() ) == m_number_of_samples ) ;
	assert( intensities.rows() == 2 ) ;
	assert( means.size() == 3 ) ;
	assert( variances.size() == 3 ) ;
	assert( non_missing_counts.size() == 3 ) ;
	for( std::size_t i = 0; i < m_number_of_samples; ++i ) {
		if( genotypes[i].is_int() && intensities( 0, i ) == intensities( 0, i ) && intensities( 1, i ) == intensities( 1, i ) ) {
			std::size_t g = genotypes[i].as< int >() ;
			means[g].noalias() += intensities.col( i ) ;
			variances[g].noalias() += intensities.col( i ) * intensities.col( i ).transpose() ;
			++non_missing_counts[g] ;
		} else {
			// something missing; ignore this.
		}
	}
	for( std::size_t g = 0; g < 3; ++g ) {
		means[g] /= non_missing_counts[g] ;
		variances[g].noalias() -= ( means[g] * means[g].transpose() ) ;
	}
}

void NormalClusterFitter::write_output_columns(
	std::vector< std::pair< std::string, std::string > > const& spec
) const {
	(*m_file) << "# Written by qctool:NormalClusterFitter, " << appcontext::get_current_time_as_string() << ".\n" ;
	(*m_file) << "SNPID rsid chromosome position" ;
	for( std::size_t i = 0 ; i < spec.size(); ++i ) {
		std::string const prefix = spec[i].first + "/" + spec[i].second ;
		for( std::size_t g = 0; g < 3; ++g ) {
			(*m_file)
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":non_missing_count"
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":mean_x"
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":mean_y"
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":variance_xx"
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":variance_xy"
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":variance_yx"
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":variance_yy"
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":T_xx"
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":T_xy"
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":T_yx"
				<< " " << prefix << ":" << genfile::string_utils::to_string( g ) << ":T_yy"
			;
		}
	}
	(*m_file) << "\n" ;
}

void NormalClusterFitter::write_id_data( genfile::SNPIdentifyingData const& id_data ) const {
	(*m_file)
		<< id_data.get_SNPID()
		<< " " << id_data.get_rsid()
		<< " " << id_data.get_position().chromosome()
		<< " " << id_data.get_position().position() ;
}

void NormalClusterFitter::write_cluster_fit(
	std::vector< Eigen::Vector2d > const& means,
	std::vector< Eigen::Matrix2d > const& variances,
	std::vector< std::size_t > const& non_missing_counts
) const {
	assert( means.size() == 3 ) ;
	assert( variances.size() == 3 ) ;
	assert( non_missing_counts.size() == 3 ) ;
	for( std::size_t g = 0; g < 3; ++g ) {
		(*m_file)
			<< " " << non_missing_counts[g] ;
		if( non_missing_counts[g] == 0 ) {
			(*m_file) << " NA NA NA NA NA NA NA NA NA NA" ;
		}
		else {
			Eigen::Matrix2d square_root = variances[g].llt().matrixL() ;
			(*m_file)
				<< " " << means[g](0)
				<< " " << means[g](1)
				<< " " << variances[g](0,0)
				<< " " << variances[g](0,1)
				<< " " << variances[g](1,0)
				<< " " << variances[g](1,1)
				<< " " << square_root(0,0)
				<< " " << square_root(0,1)
				<< " " << square_root(1,0)
				<< " " << square_root(1,1)
			;
		}
	}
}

void NormalClusterFitter::end_processing_snps() {
	
}
