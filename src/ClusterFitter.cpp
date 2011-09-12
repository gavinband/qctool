#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "ClusterFitter.hpp"
#include "NormalClusterFitter.hpp"

void ClusterFitter::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Cluster options" ) ;
	options[ "-fit-clusters" ]
		.set_description( "Using intensity and genotype data, re-fit clusters for each SNP and output them. "
		 	"The argument specifies the key of the genotype and intensity fields to work with. "
			"This should be a comma-separated list of values of the form \"genotypes/intensities\"." )
		.set_takes_single_value()
		.set_default_value( "genotypes/intensities" ) ;
		
	options[ "-fit-cluster-file" ]
		.set_description( "Override the default output file for cluster fits.")
		.set_takes_single_value()
		.set_default_value( "qctool.clusters" ) ;
	options[ "-call-threshhold" ]
		.set_description( "Set the threshhold for making genotype calls from call probabilities, where applicable." )
		.set_takes_single_value()
		.set_default_value( 0.9 ) ;
}

namespace impl {
	struct ClusterFitFileOutputter: public ClusterFitter::ResultCallback {
		static UniquePtr create( std::string const& filename ) { return UniquePtr( new ClusterFitFileOutputter( filename ) ) ; }
		ClusterFitFileOutputter( std::string const& filename ):
			m_filename( filename ),
			m_sink( statfile::BuiltInTypeStatSink::open( filename ))
		{
			(*m_sink) | "SNPID" | "rsid" | "chromosome" | "position" | "alleleA" | "alleleB" | "name" ;
			using genfile::string_utils::to_string ;
			for( std::size_t g = 0; g < 3; ++g ) {
				std::string prefix = "class" + to_string( g ) ;
				(*m_sink)
					| ( prefix + ":non_missing_count" )
					| ( prefix + ":mean_x" )
					| ( prefix + ":mean_y" )
					| ( prefix + ":variance_xx" )
					| ( prefix + ":variance_xy" )
					| ( prefix + ":variance_yx" )
					| ( prefix + ":variance_yy" )
					| ( prefix + ":T_xx" )
					| ( prefix + ":T_xy" )
					| ( prefix + ":T_yx" )
					| ( prefix + ":T_yy" ) ;
			}
		}

		virtual void set_SNP( genfile::SNPIdentifyingData const& snp ) {
			m_snp = snp ;
		}

		void write_cluster_fit(
			std::string const& name,
			std::vector< std::size_t > const& non_missing_counts,
			Eigen::MatrixXd const& fit
		) {
			(*m_sink)
				<< m_snp.get_SNPID()
				<< m_snp.get_rsid()
				<< m_snp.get_position().chromosome()
				<< m_snp.get_position().position()
				<< m_snp.get_first_allele()
				<< m_snp.get_second_allele()
				<< name ;

			for( std::size_t g = 0; g < 3; ++g ) {
				if( non_missing_counts[g] == 0 ) {
					for( std::size_t i = 0; i < 11; ++i ) {
						(*m_sink) << genfile::MissingValue() ;
					}
				}
				else {
					Eigen::Matrix2d square_root = fit.block( 0, 3*g+1, 2, 2 ).llt().matrixL() ;
					(*m_sink)
						<< non_missing_counts[g]
						<< fit( 0, 3*g )
						<< fit( 1, 3*g )
						<< fit( 0, 3*g+1 )
						<< fit( 0, 3*g+2 )
						<< fit( 1, 3*g+1 )
						<< fit( 1, 3*g+2 )
						<< square_root( 0, 0 )
						<< square_root( 0, 1 )
						<< square_root( 1, 0 )
						<< square_root( 1, 1 ) ;
				}
			}
			(*m_sink) << statfile::end_row() ;
		}

	private:
		std::string const m_filename ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
		genfile::SNPIdentifyingData m_snp ;
	} ;
}

ClusterFitter::UniquePtr ClusterFitter::create( appcontext::OptionProcessor const& options ) {
	ClusterFitter::UniquePtr result ;
	result.reset( new NormalClusterFitter( options ) ) ;
	result->connect(
		impl::ClusterFitFileOutputter::create( options.get< std::string >( "-fit-cluster-file" ))
	) ;
	return result ;
}

void ClusterFitter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& ) {
	for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
		m_callbacks[i].set_SNP( snp ) ;
	}
}

void ClusterFitter::send_results(
	std::string const& name,
	std::vector< std::size_t > const& non_missing_counts,
	Eigen::MatrixXd const& fit
) {
	for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
		m_callbacks[i].write_cluster_fit( name, non_missing_counts, fit ) ;
	}
}



