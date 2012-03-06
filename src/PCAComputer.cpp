#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/function.hpp>
#include "../config.hpp"
#include "Eigen/Core"
#include "Eigen/Eigenvalues"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SingleSNPGenotypeProbabilities.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "worker/Worker.hpp"
#include "appcontext/FileUtil.hpp"
#include "appcontext/get_current_time_as_string.hpp"
#include "appcontext/OptionProcessor.hpp"
#include "appcontext/UIContext.hpp"
#include "PCAComputer.hpp"

PCAComputer::PCAComputer(
	appcontext::OptionProcessor const& options,
	genfile::CohortIndividualSource const& samples,
	worker::Worker* worker,
	appcontext::UIContext& ui_context
):
	m_options( options ),
	m_ui_context( ui_context ),
	m_samples( samples ),
	m_number_of_samples( samples.get_number_of_individuals() ),
	m_number_of_snps_processed( 0 ),
	m_number_of_PCAs_to_compute( 0 ),
	m_threshhold( 0.9 )
{
	assert( m_options.check( "-load-kinship" )) ;
	m_filename = m_options.get< std::string >( "-load-kinship" ) ;
	// remove .csv extension if present.
	std::string const extension = ".csv" ;
	if(
		m_filename.size() >= extension.size()
		&& ( m_filename.compare(
			m_filename.size() - extension.size(),
			extension.size(),
			extension
		) == 0 )
	) {
		m_filename.resize( m_filename.size() - extension.size() ) ;
	}
	load_matrix( m_filename + ".csv", &m_kinship_matrix ) ;

	if( m_options.check( "-PCAs" ) ) {	
		m_number_of_PCAs_to_compute = std::min( m_options.get< std::size_t >( "-PCAs" ), samples.get_number_of_individuals() ) ;
	}
}

void PCAComputer::end_processing_snps() {
	// nothing to do.
}
