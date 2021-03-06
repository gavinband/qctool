
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <boost/function.hpp>
#include "config/config.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "components/SampleSummaryComponent/SampleStorage.hpp"
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
		.set_description( "Perform kinship computation using threshholded genotype calls." )
		.set_takes_single_value() ;
	options[ "-load-kinship" ]
		.set_description( "Load a previously-computed kinship matrix from the specified file." )
		.set_takes_single_value() ;
	options[ "-load-UDUT" ]
		.set_description( "Load a previously-computed eigenvalue decomposition of a relatedness matrix." )
		.set_takes_single_value() ;
	options[ "-UDUT" ]
		.set_description( "Compute the UDUT decomposition of the matrix passed to -load-kinship, and save it in the specified file." )
		.set_takes_single_value() ;
	options[ "-PCs" ]
		.set_description(
			"Compute the specified number of principal components, using a kinship matrix "
			" that is either computed using -kinship or loaded using -load-kinship. "
		 	"The argument should be the number of principal components to compute."
		)
		.set_takes_single_value()
		.set_default_value( 20 ) ;
	options[ "-loadings" ]
		.set_description( "Compute SNP loadings for each principal component, and store them in the specified file." )
		.set_takes_single_value() ;
	options[ "-project-onto" ]
		.set_description( "Project samples onto a principal components based on the loadings in the file given as the first argument. "
		 	"Put the output in the file specified in the second argument." )
		.set_takes_values_until_next_option()
		.set_maximum_multiplicity(1) ;
	options[ "-kinship-method" ]
		.set_description( "Method to use for relatedness matrix computation."
			"The only permitted value is \"lookup-table\", which uses a lookup table to compute kinship values across"
			" several SNPs at a time.  This is usually fastest."
		)
		.set_takes_single_value()
		.set_default_value( "lookup-table" )
	;

	options.option_implies_option( "-kinship", "-s" ) ;
	options.option_implies_option( "-kinship-method", "-kinship" ) ;
	options.option_implies_option( "-load-kinship", "-s" ) ;
	options.option_implies_option( "-UDUT", "-PCs" ) ;
	options.option_implies_option( "-PCs", "-UDUT" ) ;
	options.option_implies_option( "-PCs", "-osample" ) ;
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


namespace impl {
	std::string get_ids( genfile::CohortIndividualSource const* samples, std::size_t i ) {
		return samples->get_entry( i, "id_1" ).as< std::string >() ;
	}

	std::string get_tagged_ids( genfile::CohortIndividualSource const* samples, std::size_t i, std::vector< std::string > const& tags ) {
		std::size_t const sample_i = i / tags.size() ;
		std::size_t const tag_i = i % tags.size() ;
		return samples->get_entry( sample_i, "id_1" ).as< std::string >()
			+ ":"
			+ tags[tag_i]
		;
	}
}

void RelatednessComponent::setup(
	genfile::SNPDataSourceProcessor& processor,
	sample_stats::SampleStorage::SharedPtr storage
) const {
	PCAComputer::SharedPtr pca_computer ;
	KinshipCoefficientManager::GetNames get_ids = boost::bind(
		&impl::get_ids,
		&m_samples,
		_1
	) ;

	if( m_options.check( "-PCs" ) ) {
		if( !storage ) {
			throw genfile::BadArgumentError(
				"RelatednessComponent::setup()",
				"storage",
				"A valid storage destination for PCs must be supplied."
			) ;
		}
		if( !m_options.check( "-kinship" ) && !m_options.check( "-load-kinship" )) {
			throw genfile::BadArgumentError(
				"RelatednessComponent::setup()",
				"-PCs",
				"Expected -kinship or -load-kinship to be supplied."
			) ;
		}
		pca_computer = PCAComputer::SharedPtr( new PCAComputer( m_options, m_samples, m_ui_context ) ) ;
		pca_computer->send_UDUT_to(
			boost::bind(
				&pca::write_matrix,
				m_options.get< std::string >( "-UDUT" ),
				_3, "qctool:PCAComputer" ,_1, _4, _5
			)
		) ;
		pca_computer->send_PCs_to( storage ) ;
	}

	if( m_options.check( "-loadings" )) {
		PCALoadingComputer::UniquePtr loading_computer( new PCALoadingComputer( m_options.get< int >( "-PCs" ) ) ) ;
		if( pca_computer.get() ) {
			pca_computer->send_UDUT_to( boost::bind( &PCALoadingComputer::set_UDUT, loading_computer.get(), _2, _3 ) ) ;
		} else if( m_options.check( "-load-UDUT" )) {
			relatedness::UDUTDecompositionLoader udut( m_samples ) ;
			udut.send_UDUT_to( boost::bind( &PCALoadingComputer::set_UDUT, loading_computer.get(), _3, _1 )) ;
			udut.load_matrix( m_options.get< std::string >( "-load-UDUT" )) ;
		} else {
			throw appcontext::OptionProcessorImpliedOptionNotSuppliedException( "-loadings", "-PCs or -load-UDUT" ) ;
		}

		// Set up an output location for the loadings.
		boost::shared_ptr< statfile::BuiltInTypeStatSink > loadings_file(
			statfile::BuiltInTypeStatSink::open( m_options.get< std::string >( "-loadings" ) ).release()
		) ;
	
		loadings_file->write_metadata( 
			pca::get_metadata( "PCALoadingComputer", loading_computer->get_metadata() )
		) ;

		loading_computer->send_results_to(
			boost::bind(
				&pca::write_loadings_to_sink,
				loadings_file,
				_1, _2, _3 ,_4, _5
			)
		) ;

		processor.add_callback(
			genfile::SNPDataSourceProcessor::Callback::UniquePtr( loading_computer.release() )
		) ;
	}
	
	if( m_options.check( "-kinship" )) {
		KinshipCoefficientComputer::Computation::UniquePtr computation ;
		std::string const& method = m_options.get< std::string >( "-kinship-method" ) ;
		if( method == "lookup-table" ) {
			computation.reset(
				impl::NormaliseGenotypesAndComputeXXtFast::create(
					m_worker, 4
				).release()
			) ;
		} else {
			throw genfile::BadArgumentError(
				"RelatednessComponent::setup()",
				"-method " + method,
				"Method must be \"lookup-table\""
			) ;
		}
		
		m_ui_context.logger() << "RelatednessComponent::setup(): using computation: " << computation->get_summary() << ".\n" ;

		KinshipCoefficientComputer::UniquePtr result(
			new KinshipCoefficientComputer(
				m_options,
				m_samples,
				m_ui_context,
				computation
			)
		) ;
		result->send_results_to(
			boost::bind(
				&pca::write_matrix_lower_diagonals_in_long_form,
				m_options.get< std::string >( "-kinship" ),
				_2, _3, _4, _5,
				get_ids, get_ids
			)
		) ;
		if( pca_computer ) {
			result->send_results_to(
				boost::bind(
					&PCAComputer::compute,
					pca_computer,
					_3, _1, m_options.get< std::string >( "-kinship" )
				)
			) ;
		}
		processor.add_callback(
			genfile::SNPDataSourceProcessor::Callback::UniquePtr( result.release() )
		) ;
	} else if( m_options.check( "-load-kinship" )) {
		assert( pca_computer.get() ) ;
		assert( !m_options.check( "-kinship" ) ) ;
		Eigen::MatrixXd relatednessMatrix ;
		std::size_t number_of_snps ;
		PCAComputer::load_long_form_matrix(
			m_samples,
			m_options.get< std::string >( "-load-kinship" ),
			&relatednessMatrix,
			&number_of_snps,
			m_ui_context
		) ;
		pca_computer->compute( relatednessMatrix, number_of_snps, m_options.get< std::string >( "-load-kinship" ) ) ;
	}


	if( m_options.check( "-project-onto" )) {
		std::vector< std::string > elts = m_options.get_values< std::string >( "-project-onto" ) ;
		assert( elts.size() >= 2 && elts.size() <= 3 ) ;
		std::string const match_fields = ( elts.size() == 3 ) ? elts[2] : m_options.get_value< std::string >( "-compare-variants-by" ) ;
		pca::PCAProjector::UniquePtr projector = pca::PCAProjector::create(
			m_samples,
			m_ui_context,
			genfile::VariantIdentifyingData::CompareFields( match_fields )
		) ;
		projector->send_results_to(
			boost::bind(
				&pca::write_sample_file,
				elts[1],
				boost::cref( m_samples ),
				_2, "qctool:PCAProjector", _1, _3, _4
			)
		) ;
		pca::PCALoadingLoader loader( m_samples ) ;
		loader.send_loadings_to( boost::bind( &pca::PCAProjector::set_loadings, projector.get(), _1, _2, _3, _4, _7 )) ;
		loader.load_loadings( elts[0] ) ;
		
		processor.add_callback(
			genfile::SNPDataSourceProcessor::Callback::UniquePtr( projector.release() )
		) ;
	}
	
	
}
