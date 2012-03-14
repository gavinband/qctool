#include <boost/function.hpp>
#include "../config.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "worker/Task.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "components/RelatednessComponent/RelatednessComponent.hpp"
#include "components/RelatednessComponent/PCAComputer.hpp"
#include "components/RelatednessComponent/KinshipCoefficientComputer.hpp"
#include "components/RelatednessComponent/PCALoadingComputer.hpp"
#include "components/RelatednessComponent/UDUTDecompositionLoader.hpp"
#include "components/RelatednessComponent/PCALoadingLoader.hpp"
#include "components/RelatednessComponent/PCAProjector.hpp"
#include "components/RelatednessComponent/names.hpp"
#include "components/RelatednessComponent/write_matrix.hpp"

void RelatednessComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "Kinship options" ) ;
	options[ "-kinship" ]
		.set_description( "Perform kinship computation using threshholded genotype calls and cblas or Eigen libraries." )
		.set_takes_single_value() ;
	options[ "-load-kinship" ]
		.set_description( "Load a previously-computed kinship matrix from the specified file." )
		.set_takes_single_value() ;
	options[ "-load-UDUT" ]
		.set_description( "Load a previously-computed eigenvalue decomposition of a relatedness matrix." )
		.set_takes_single_value() ;
	options[ "-UDUT" ]
		.set_description( "Compute the UDUT decomposition of the matrix passed to -load-kinship." ) ;
	options[ "-PCAs" ]
		.set_description( "Compute PCA components of kinship matrix. "
		 	"The argument should be the number of PCA components to compute.")
		.set_takes_single_value()
		.set_default_value( 10 ) ;
	options[ "-PCA-prefix" ]
		.set_description( "Set the prefix for filenames used to store eigendecomposition and PCA-related files. "
		 	"By default, this is formed by removing the extension from the argument to -load-kinship." )
		.set_takes_single_value() ;
	options[ "-loadings" ]
		.set_description( "Compute SNP loadings for each PCA." )
		.set_takes_single_value() ;
	options[ "-project-onto" ]
		.set_description( "Project samples onto a PCA based on the loadings in the file given as the first argument. "
		 	"Put the output in the file specified in the second argument." )
		.set_takes_values(2)
		.set_maximum_multiplicity(1) ;
	options[ "-PCA-exclusions" ]
		.set_description( "Output a list of exclusions based on outliers in the first N PCA components." )
		.set_takes_single_value() ;
	options[ "-no-lapack" ]
		.set_description( "Don't use lapack to perform computations. "
			"Usually this means to use the Eigen library (http://eigen.tuxfamily.org) instead." ) ;

	options.option_implies_option( "-kinship", "-s" ) ;
	options.option_implies_option( "-load-kinship", "-s" ) ;
	options.option_implies_option( "-UDUT", "-load-kinship" ) ;
	options.option_implies_option( "-PCAs", "-UDUT" ) ;
	options.option_implies_option( "-UDUT", "-PCAs" ) ;
	options.option_implies_option( "-PCA-prefix", "-PCAs" ) ;
	options.option_excludes_option( "-load-kinship", "-kinship" ) ;
	options.option_excludes_option( "-project-onto", "-loadings" ) ;
	options.option_excludes_option( "-load-UDUT", "-UDUT" ) ;
}

RelatednessComponent::UniquePtr RelatednessComponent::create(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	worker::Worker* worker,
	appcontext::UIContext& ui_context
) {
	RelatednessComponent::UniquePtr result(
		new RelatednessComponent( options, samples, worker, ui_context )
	) ;
	return result ;
}

RelatednessComponent::RelatednessComponent(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	worker::Worker* worker,
	appcontext::UIContext& ui_context
):
	m_options( options ),
	m_samples( samples ),
	m_worker( worker ),
	m_ui_context( ui_context )
{}


void RelatednessComponent::setup( genfile::SNPDataSourceProcessor& processor ) const {
	PCAComputer::UniquePtr pca_computer ;
	if( m_options.check( "-kinship" )) {
		KinshipCoefficientComputer::UniquePtr result( new KinshipCoefficientComputer( m_options, m_samples, m_worker, m_ui_context ) ) ;
		result->send_results_to(
			boost::bind(
				&pca::write_matrix,
				m_options.get< std::string >( "-kinship" ),
				_2, _3, _4, _5, _6
			)
		) ;
		processor.add_callback(
			genfile::SNPDataSourceProcessor::Callback::UniquePtr( result.release() )
		) ;
	}
	else if( m_options.check( "-PCAs" ) ) {
		pca_computer.reset( new PCAComputer( m_options, m_samples, m_worker, m_ui_context ) ) ;
		pca_computer->send_UDUT_to(
			boost::bind(
				&pca::write_matrix,
				get_PCA_filename_prefix() + ".UDUT.csv",
				_3, "qctool:PCAComputer" ,_1, _4, _5
			)
		) ;
		pca_computer->send_PCAs_to(
			boost::bind(
				&pca::write_matrix,
				get_PCA_filename_prefix() + ".PCAs.csv",
				_3, "qctool:PCAComputer", _1, _4, _5
			)
		) ;
	}
	if( m_options.check( "-loadings" )) {
		PCALoadingComputer::UniquePtr loading_computer( new PCALoadingComputer( m_options.get< int >( "-loadings" ) ) ) ;
		if( pca_computer.get() ) {
			pca_computer->send_UDUT_to( boost::bind( &PCALoadingComputer::set_UDUT, loading_computer.get(), _2, _3 ) ) ;
			pca_computer->compute_PCA() ;
		} else if( m_options.check( "-load-UDUT" )) {
			relatedness::UDUTDecompositionLoader udut( m_samples ) ;
			udut.send_UDUT_to( boost::bind( &PCALoadingComputer::set_UDUT, loading_computer.get(), _3, _1 )) ;
			udut.load_matrix( m_options.get< std::string >( "-load-UDUT" )) ;
		} else {
			throw appcontext::OptionProcessorImpliedOptionNotSuppliedException( "-loadings", "-PCAs or -load-UDUT" ) ;
		}

		// Set up an output location for the loadings.
		boost::shared_ptr< statfile::BuiltInTypeStatSink > loadings_file(
			statfile::BuiltInTypeStatSink::open( get_PCA_filename_prefix() + ".loadings.csv" ).release()
		) ;
		
		loadings_file->write_metadata( 
			pca::get_metadata( "PCALoadingComputer", loading_computer->get_metadata() )
		) ;

		loading_computer->send_results_to(
			boost::bind(
				&pca::write_snp_and_vector_to_sink,
				loadings_file,
				_1, _2, _3
			)
		) ;

		processor.add_callback(
			genfile::SNPDataSourceProcessor::Callback::UniquePtr( loading_computer.release() )
		) ;
	} else if( pca_computer.get() ) {
		pca_computer->compute_PCA() ;
	}
	
	if( m_options.check( "-project-onto" )) {
		std::vector< std::string > elts = m_options.get_values< std::string >( "-project-onto" ) ;
		assert( elts.size() == 2 ) ;
		pca::PCAProjector::UniquePtr projector = pca::PCAProjector::create( m_samples, m_ui_context ) ;
		projector->send_results_to(
			boost::bind(
				&pca::write_matrix,
				elts[1],
				_2, "qctool:PCAProjector", _1, _3, _4
			)
		) ;
		pca::PCALoadingLoader loader( m_samples ) ;
		loader.send_loadings_to( boost::bind( &pca::PCAProjector::set_loadings, projector.get(), _1, _2, _5 )) ;
		loader.load_loadings( elts[0] ) ;
		
		processor.add_callback(
			genfile::SNPDataSourceProcessor::Callback::UniquePtr( projector.release() )
		) ;
	}
	
	
}

std::string RelatednessComponent::get_PCA_filename_prefix() const {
	std::string result ;
	if( m_options.check( "-PCA-prefix" ) ) {
		result = m_options.get< std::string >( "-PCA-prefix" ) ;
	} else if( m_options.check( "-load-kinship" )) {
		result = genfile::replace_or_add_extension( m_options.get< std::string >( "-load-kinship" ), "" ) ;
	} else if( m_options.check( "-load-UDUT" )) {
		result = genfile::replace_or_add_extension( m_options.get< std::string >( "-load-UDUT" ), "" ) ;
		if( result.size() > 5 && result.substr( result.size() - 5, result.size() ) == ".UDUT" ) {
			result = result.substr( 0, result.size() - 5 ) ;
		}
	}
	return result ;
}



