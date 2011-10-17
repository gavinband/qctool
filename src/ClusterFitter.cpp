#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/Error.hpp"
#include "genfile/endianness_utils.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "ClusterFitter.hpp"
#include "NormalClusterFitter.hpp"
#include "DataStore.hpp"

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
		}

		virtual void set_SNP( genfile::SNPIdentifyingData const& snp ) {
			m_snp = snp ;
		}

		void write_intensity_model(
			genfile::SNPIdentifyingData const& snp,
			std::string const& name,
			Eigen::MatrixXd const& intensities,
			genfile::SingleSNPGenotypeProbabilities const& genotypes,
			IntensityModel::SharedPtr model
		) {
			assert(0) ;
			(*m_sink)
				<< m_snp.get_SNPID()
				<< m_snp.get_rsid()
				<< m_snp.get_position().chromosome()
				<< m_snp.get_position().position()
				<< m_snp.get_first_allele()
				<< m_snp.get_second_allele()
				<< name ;
			(*m_sink) << statfile::end_row() ;
		}

	private:
		std::string const m_filename ;
		statfile::BuiltInTypeStatSink::UniquePtr m_sink ;
		genfile::SNPIdentifyingData m_snp ;
	} ;
	
	struct ClusterFitDataStoreOutputter: public ClusterFitter::ResultCallback {
		static UniquePtr create( std::string const& spec ) {
			UniquePtr result ;
			result.reset(
				new ClusterFitDataStoreOutputter(
					DataStore::create( spec )
				)
			) ;
			return result ;
		}

		ClusterFitDataStoreOutputter( DataStore::UniquePtr store ):
			m_store( store ),
			m_storage_id( m_store->get_or_create_entity( "double_matrix", "A matrix of doubles" ) ),
			m_cohort_id( m_store->get_or_create_entity( "unnamed cohort", "An unnamed cohort" ) )
		{
		}
		
		virtual void set_SNP( genfile::SNPIdentifyingData const& snp ) {
			m_snp_id = m_store->get_or_create_SNP( snp ) ;
		}

		void write_intensity_model(
			genfile::SNPIdentifyingData const& snp,
			std::string const& name,
			Eigen::MatrixXd const& intensities,
			genfile::SingleSNPGenotypeProbabilities const& genotypes,
			IntensityModel::SharedPtr model
		) {
			assert(0) ;
		}

	private:
		DataStore::UniquePtr m_store ;
		db::Connection::RowId m_snp_id ;
		db::Connection::RowId m_storage_id ;
		db::Connection::RowId m_cohort_id ;
	} ;
}

ClusterFitter::UniquePtr ClusterFitter::create( appcontext::OptionProcessor const& options ) {
	ClusterFitter::UniquePtr result ;
	result.reset( new NormalClusterFitter( options ) ) ;
	std::string filename = options.get< std::string >( "-fit-cluster-file" ) ;
	try {
		db::SQLite3Connection connection( filename, false ) ;
		result->connect(
			impl::ClusterFitDataStoreOutputter::create( "sqlite3://" + filename )
		) ;
	}
	catch( db::Error const& ) {
		// not a sqlite3 file.  Use a file instead.
		result->connect(
			impl::ClusterFitFileOutputter::create( filename )
		) ;
	}
	return result ;
}

void ClusterFitter::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& ) {
	for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
		m_callbacks[i].set_SNP( snp ) ;
	}
}

void ClusterFitter::send_results(
	genfile::SNPIdentifyingData const& snp,
	std::string const& name,
	Eigen::MatrixXd const& intensities,
	genfile::SingleSNPGenotypeProbabilities const& genotypes,
	IntensityModel::SharedPtr model
) {
	for( std::size_t i = 0; i < m_callbacks.size(); ++i ) {
		m_callbacks[i].write_intensity_model( snp, name, intensities, genotypes, model ) ;
	}
}



