#include <iostream>
#include <iomanip>
#include <string>
#include <Eigen/Core>
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/Error.hpp"
#include "statfile/BuiltInTypeStatSource.hpp"
#include "components/RelatednessComponent/PCALoadingLoader.hpp"

namespace pca {
	PCALoadingLoader::PCALoadingLoader( genfile::CohortIndividualSource const& samples ):
	 	m_samples( samples )
	{}

	void PCALoadingLoader::send_loadings_to( ResultSignal::slot_type callback ) {
		m_result_signal.connect( callback ) ;
	}

	void PCALoadingLoader::load_loadings( std::string const& filename ) const {
		Eigen::MatrixXd loadings ;
		std::vector< genfile::SNPIdentifyingData > snps ;
		std::vector< std::string > column_names ;
		load_loadings_impl( filename, &snps, &loadings, &column_names ) ;
		m_result_signal( snps, loadings, snps.size(), m_samples.get_number_of_individuals(), column_names ) ;
	}

	void PCALoadingLoader::load_loadings_impl( std::string const& filename, std::vector< genfile::SNPIdentifyingData >* snps, Eigen::MatrixXd* loadings, std::vector< std::string >* column_names ) const {
		statfile::BuiltInTypeStatSource::UniquePtr source = statfile::BuiltInTypeStatSource::open( filename ) ;
		assert( loadings ) ;
		assert( snps ) ;
		
		using namespace genfile::string_utils ;
		// Read the metadata
		std::size_t number_of_samples = 0 ;
		std::size_t number_of_snps = 0 ;
		std::vector< std::string > metadata = split( source->get_descriptive_text(), "\n" ) ;
		for( std::size_t i = 0; i < metadata.size(); ++i ) {
			std::vector< std::string > elts = split_and_strip( metadata[i], ":", " " ) ;
			if( elts.size() == 2 && elts[0] == "Number of SNPs" ) {
				number_of_snps = to_repr< std::size_t >( elts[1] ) ;
			}
			else if( elts.size() == 2 && elts[0] == "Number of samples" ) {
				number_of_samples = to_repr< std::size_t >( elts[1] ) ;
			}
		}
		if( source->number_of_columns() < 7 ) {
			throw genfile::MalformedInputError( source->get_source_spec(), 0, source->number_of_columns() ) ;
		}
		if( source->number_of_rows() != number_of_snps ) {
			throw genfile::MalformedInputError( source->get_source_spec(), 0, std::min( source->number_of_rows(), number_of_samples )) ;
		}
		if( ( source->number_of_columns() - 6 ) % 2 != 0 ) {
			throw genfile::MalformedInputError( source->get_source_spec(), 0, source->number_of_columns() ) ;
		}
		std::size_t const number_of_loadings = ( source->number_of_columns() - 6 ) / 2;

		snps->resize( number_of_snps ) ;
		loadings->resize( number_of_snps, number_of_loadings ) ;
		for( std::size_t i = 0; i < number_of_snps; ++i ) {
			(*source) >> (*snps)[i].SNPID() >> (*snps)[i].rsid()
				>> (*snps)[i].position().chromosome() >> (*snps)[i].position().position()
				>> (*snps)[i].first_allele() >> (*snps)[i].second_allele() ;
			for( std::size_t j = 0; j < number_of_loadings; ++j ) {
				(*source) >> (*loadings)(i,j) ;
			}
			(*source) >> statfile::ignore_all() ;
		}
		
		assert( column_names ) ;
		assert( column_names->empty() ) ;
		for( std::size_t i = 0; i < number_of_loadings; ++i ) {
			column_names->push_back( source->name_of_column( 6 + i )) ;
		}
	}
}
