
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "components/RelatednessComponent/PCAProjector.hpp"
#include "components/RelatednessComponent/PCAComputer.hpp"
#include "components/RelatednessComponent/mean_centre_genotypes.hpp"
#include "components/RelatednessComponent/names.hpp"
#include "components/RelatednessComponent/write_matrix_to_stream.hpp"

// #define DEBUG_PCA_PROJECTOR 1

namespace pca {
	PCAProjector::UniquePtr PCAProjector::create( genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ) {
		return PCAProjector::UniquePtr( new PCAProjector( samples, ui_context )) ;
	}

	PCAProjector::PCAProjector( genfile::CohortIndividualSource const& samples, appcontext::UIContext& ui_context ):
		m_samples( samples ),
	 	m_ui_context( ui_context )
	{}

	void PCAProjector::set_loadings( std::vector< genfile::SNPIdentifyingData > const& snps, Matrix const& loadings, std::vector< std::string > const& names ) {
		using genfile::string_utils::to_string ;
		if( snps.size() != std::size_t( loadings.rows() ) ) {
			throw genfile::MismatchError(
				"relatedness::PCAProjector::set_loadings()",
				"(function arguments)",
				"snp.size()=" + to_string( snps.size() ),
				"loadings.rows()=" + to_string( loadings.rows() )
			) ;
		}
		if( names.size() != std::size_t( loadings.cols() ) ) {
			throw genfile::MismatchError(
				"relatedness::PCAProjector::set_loadings()",
				"(function arguments)",
				"names.size()=" + to_string( names.size() ),
				"loadings.cols()=" + to_string( loadings.cols() )
			) ;
		}
		for( std::size_t i = 0; i < snps.size(); ++i ) {
			std::pair< SnpMap::const_iterator, bool > where = m_snps.insert( std::make_pair( snps[i].get_position(), int( i ) ) ) ;
			if( !where.second ) {
				throw genfile::DuplicateKeyError( "snps argument to pca::PCAProjector::set_loadings()", to_string( snps[i] )) ;
			}
		}
		m_loadings = loadings ;
		m_names = names ;
		
#if DEBUG_PCA_PROJECTOR
		std::cerr << "PCAProjector: loadings are:\n" ;
		pca::write_matrix_to_stream( std::cerr, m_loadings.block( 0, 0, 10, m_loadings.cols() )) ;
#endif
	}

	void PCAProjector::send_results_to( ResultSignal::slot_type callback ) {
		m_result_signal.connect( callback ) ;
	}

	void PCAProjector::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) {
		m_projections.setZero( number_of_samples, m_loadings.cols() ) ;
	}

	void PCAProjector::processed_snp( genfile::SNPIdentifyingData const& snp, genfile::VariantDataReader& data_reader ) {
		SnpMap::const_iterator where = m_snps.find( snp.get_position() ) ;
		if( where != m_snps.end() && ( m_loadings.row( where->second ).sum() == m_loadings.row( where->second ).sum() ) ) {
			data_reader.get(
				":genotypes:",
				genfile::vcf::get_threshholded_calls( m_genotype_calls, m_non_missingness, 0.9, 0, 0, 1, 2 )
			) ;
			assert( m_genotype_calls.size() == m_projections.rows() ) ;
			assert( m_non_missingness.size() == m_genotype_calls.size() ) ;
			double const allele_frequency = m_genotype_calls.sum() / ( 2.0 * m_non_missingness.sum() ) ;
			double const maf = std::min( allele_frequency, 1 - allele_frequency ) ;
			if( maf > 0.001 ) {
#if DEBUG_PCA_PROJECTOR
				std::cerr << std::resetiosflags( std::ios::floatfield ) << std::setprecision( 5 ) ;
				std::cerr << "SNP: " << snp << ": freq = " << allele_frequency << ", uncentred genotypes are: " << m_genotype_calls.transpose().head( 20 ) << "...\n" ;
#endif
				pca::mean_centre_genotypes( &m_genotype_calls, m_non_missingness, allele_frequency ) ;
				m_genotype_calls /= std::sqrt( 2.0 * allele_frequency * ( 1.0 - allele_frequency ) ) ;
#if DEBUG_PCA_PROJECTOR
				std::cerr << "mean-centred genotypes are: " << m_genotype_calls.transpose().head( 20 ) << "...\n" ;
				std::cerr << "non-missingness is: " << m_non_missingness.transpose().head( 20 ) << "...\n" ;
#endif
				
				//
				// We have loadings for C components as columns of m_loadings.
				// We have calls for N samples as entries of m_genotype_calls
				// We wish to get a matrix with samples as rows and projections as columns.
				// We therefore compute a Kronecker product of the column vectors of genotype calls
				// with the row vector of loadings for that SNP.
				//
				m_projections += m_genotype_calls * m_loadings.row( where->second ) ;
				assert( m_projections.sum() == m_projections.sum() ) ;
			}
			m_visited[ snp.get_position() ] = true ;
		}
	}

	void PCAProjector::end_processing_snps() {
		diagnose_projection() ;
		using genfile::string_utils::to_string ;
		m_result_signal(
			"Number of SNPs: " + to_string( m_snps.size() ) + "\n" +
			"Number of samples: " + to_string( m_projections.size() ) + "\n" +
			"Note: these PCAs are 1/sqrt(L) times the projection of samples onto unit eigenvectors of the variance-covariance matrix\n"
			"    1/(L-1) X X^t,\n"
			"where X is the L x N matrix of genotypes, L is the number of SNPs, and N the number of samples.\n"
			"The constant 1/sqrt(L) ensures that the PCAs do not grow with the number of SNPs.",
			m_projections / std::sqrt( m_snps.size() ),
			boost::bind(
				&pca::get_concatenated_sample_ids,
				&m_samples,
				_1
			),
			boost::bind(
				&PCAComputer::get_pca_name,
				_1
			)
		) ;
	}

	void PCAProjector::diagnose_projection() const {
		{
			std::size_t missing_count = 0 ;
			for( SnpMap::const_iterator i = m_snps.begin(); i != m_snps.end(); ++i ) {
				VisitedSnpMap::const_iterator where = m_visited.find( i->first ) ;
				if( where == m_visited.end() ) {
					if( missing_count == 0 ) {
						m_ui_context.logger() << "!! ( pca::PCAProjector::diagnose_projection() ): "
							<< "The following SNPs with loadings were not used in the projection computation:\n" ;
					}
					m_ui_context.logger() << "    " << i->first << "\n" ;
					++missing_count ;
				}
			}
			m_ui_context.logger() << "++ ( pca::PCAProjector::diagnose_projection() ): a total of " << m_snps.size() - missing_count
				<< " of " << m_snps.size() << " loading SNPs were used in the projection computation.\n" ;
		}
	}
	
	std::string PCAProjector::get_metadata() const {
		using namespace genfile::string_utils ;
		return "Number of SNPs: " + to_string( m_snps.size() ) ;
	}
	
}
