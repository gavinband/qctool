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
	options[ "-PCAs" ]
		.set_description( "Compute PCA components of kinship matrix. "
		 	"The argument should be the number of PCA components to compute.")
		.set_takes_single_value()
		.set_default_value( 10 ) ;
	options[ "-PCA-prefix" ]
		.set_description( "Set the prefix for filenames used to store eigendecomposition and PCA-related files. "
		 	"By default, this is formed by removing the extension from the argument to -load-kinship." )
		.set_takes_single_value() ;
		
	options[ "-PCA-exclusions" ]
		.set_description( "Output a list of exclusions based on outliers in the first N PCA components." )
		.set_takes_single_value() ;
	options[ "-loadings" ]
		.set_description( "Compute SNP loadings for each PCA." )
		.set_takes_single_value() ;
	options[ "-no-lapack" ]
		.set_description( "Don't use lapack to perform computations. "
			"Usually this means to use the Eigen library (http://eigen.tuxfamily.org) instead." ) ;

	options.option_implies_option( "-kinship", "-s" ) ;
	options.option_implies_option( "-load-kinship", "-s" ) ;
	options.option_implies_option( "-PCAs", "-load-kinship" ) ;
	options.option_implies_option( "-PCA-prefix", "-PCAs" ) ;
	options.option_excludes_option( "-load-kinship", "-kinship" ) ;
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
	template< typename Matrix >
	void write_matrix_as_csv(
		std::string const& filename,
		Matrix const& matrix,
		std::string const& source,
		std::string const& description,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_row_names = 0,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_column_names = 0
	) {
		std::cerr << "write_matrix_as_csv: description is \"" << description << "\".\n" ;
		appcontext::OUTPUT_FILE_PTR file = appcontext::open_file_for_output( filename ) ;
		(*file) << "# Created by " << source << ", " << appcontext::get_current_time_as_string() << "\n" ;
		if( description.size() > 0 ) {
			bool newLine = true ;
			for( std::size_t i = 0; i < description.size(); ++i ) {
				if( newLine ) {
					(*file) << "# " ;
				}
				(*file) << description[i] ;
				newLine = ( description[i] == '\n' ) ;
			}
		}
		(*file) << "\n" ;
		if( get_column_names ) {
			if( get_row_names ) {
				(*file) << "id" ;
			}
			for( int j = 0; j < matrix.cols(); ++j ) {
				if( get_row_names || ( j > 0 )) {
					(*file) << "," ;
				}
				(*file) << get_column_names( j ) ;
			}
			(*file) << "\n" ;
		}

		for( int i = 0; i < matrix.rows(); ++i ) {
			if( (!get_row_names.empty()) ) {
				(*file) << get_row_names(i) ;
			}
			for( int j = 0; j < matrix.cols(); ++j ) {
				if( (!get_row_names.empty()) || ( j > 0 ) ) {
					(*file) << "," ;
				}
				double const& value = matrix( i, j ) ;
				if( value == value ) {
					(*file) << value ;
				}
				else {
					(*file) << "NA" ;
				}
			}
			(*file) << "\n" ;
		}
	}
	
	template< typename Vector >
	void write_snp_and_vector(
		boost::shared_ptr< statfile::BuiltInTypeStatSink > sink,
		genfile::SNPIdentifyingData snp,
		Vector const& vector,
		boost::function< genfile::VariantEntry ( std::size_t ) > get_names
	) {
		if( sink->number_of_rows_written() == 0 && sink->current_column() == 0 ) {
			(*sink) | "SNPID" | "rsid" | "chromosome" | "position" | "allele_A" | "allele_B" ;
			if( get_names ) {
				for( int i = 0; i < vector.size(); ++i ) {
					(*sink) | get_names( i ).as< std::string >() ;
				}
			}
			else {
				for( int i = 0; i < vector.size(); ++i ) {
					(*sink) | ( "v" + genfile::string_utils::to_string( i )) ;
				}
			}
		}
		(*sink)
			<< snp.get_SNPID()
			<< snp.get_rsid()
			<< std::string( snp.get_position().chromosome() )
			<< snp.get_position().position()
			<< snp.get_first_allele()
			<< snp.get_second_allele() ;
		for( int i = 0; i < vector.size(); ++i ) {
			(*sink) << vector(i) ;
		}
		(*sink) << statfile::end_row() ;
	}
}

void RelatednessComponent::setup( genfile::SNPDataSourceProcessor& processor ) const {
	PCAComputer::UniquePtr pca_computer ;
	if( m_options.check( "-kinship" )) {
		KinshipCoefficientComputer::UniquePtr result( new KinshipCoefficientComputer( m_options, m_samples, m_worker, m_ui_context ) ) ;
		result->send_results_to( &impl::write_matrix_as_csv< Eigen::MatrixXd > ) ;
		processor.add_callback(
			genfile::SNPDataSourceProcessor::Callback::UniquePtr( result.release() )
		) ;
	}
	else if( m_options.check( "-PCAs" ) ) {
		pca_computer.reset( new PCAComputer( m_options, m_samples, m_worker, m_ui_context ) ) ;
		pca_computer->send_UDUT_to(
			boost::bind(
				&impl::write_matrix_as_csv< Eigen::MatrixXd >,
				get_PCA_filename_prefix() + ".UDUT.csv",
				_3, "qctool:PCAComputer" ,_1, _4, _5
			)
		) ;
		pca_computer->send_PCAs_to(
			boost::bind(
				&impl::write_matrix_as_csv< Eigen::MatrixXd >,
				get_PCA_filename_prefix() + ".PCAs.csv",
				_3, "qctool:PCAComputer", _1, _4, _5
			)
		) ;
	}
	if( m_options.check( "-loadings" )) {
		PCALoadingComputer::UniquePtr loading_computer( new PCALoadingComputer( m_options.get< int >( "-loadings" ) ) ) ;
		if( pca_computer.get() ) {
			pca_computer->send_UDUT_to( boost::bind( &PCALoadingComputer::set_UDUT, loading_computer.get(), _2, _3 ) ) ;
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
			"Created by qctool::PCALoadingComputer, " + appcontext::get_current_time_as_string() + "\n"
			+ loading_computer->get_metadata()
		) ;

		loading_computer->send_results_to(
			boost::bind(
				&impl::write_snp_and_vector< Eigen::VectorXd >,
				loadings_file,
				_1, _2, _3
			)
		) ;

		processor.add_callback(
			genfile::SNPDataSourceProcessor::Callback::UniquePtr( loading_computer.release() )
		) ;
	}
	
	if( m_options.check( "-project-onto" )) {
		pca::PCALoadingLoader loader( m_samples ) ;
	}
	
	
	if( pca_computer.get() ) {
		pca_computer->compute_PCA() ;
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



