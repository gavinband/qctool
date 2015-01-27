
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <Eigen/Core>
#include <strstream>
#include "metro/mean_and_variance.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "components/SampleSummaryComponent/SampleSummaryComputation.hpp"
#include "components/SampleSummaryComponent/RiskScoreComputation.hpp"

//#define DEBUG_RISK_SCORE_COMPUTATION 1

namespace sample_stats {
	RiskScoreComputation::RiskScoreComputation( genfile::CohortIndividualSource const& samples, genfile::SNPIdentifyingData::CompareFields comparator ):
		m_samples( samples ),
		m_number_of_samples( samples.get_number_of_individuals() ),
		m_map( comparator )
	{}

	void RiskScoreComputation::accumulate( genfile::SNPIdentifyingData const& snp, Genotypes const& genotypes, genfile::VariantDataReader& ) {
		// do nothing
		RiskScoreMap::const_iterator const where = m_map.find( snp ) ;
		if( where != m_map.end() ) {
			RiskScoreIdBetaMap::const_iterator i = where->second.begin() ;
			RiskScoreIdBetaMap::const_iterator const end_i = where->second.end() ;
			for( ; i != end_i; ++i ) {
				std::string const& risk_score_identifier = i->first ;
#if	DEBUG_RISK_SCORE_COMPUTATION
				std::cerr << "Accumulating risk scores for " << snp << ", " << risk_score_identifier << "...\n" ;
#endif		
				accumulate( genotypes, i->second, m_scores[ risk_score_identifier ], m_counts[ risk_score_identifier ] ) ;
			}
		}
	}

	void RiskScoreComputation::accumulate( Genotypes const& genotypes, Betas const& betas, Eigen::VectorXd& scores, Eigen::VectorXd& counts ) const {
		assert( genotypes.cols() == 3 ) ;
		Eigen::VectorXd non_missingness = ( genotypes.rowwise().sum().array() > 0.0 ).cast< double >() ;
		Eigen::VectorXd contribution = ( genotypes * betas ).rowwise().sum() ;
		
#if DEBUG_RISK_SCORE_COMPUTATION > 1
		std::cerr << "betas = \n"
			<< betas << ".\n" ;
		std::cerr << "non_missingness = " << non_missingness.head(10).transpose() << "...\n" ;
		std::cerr << "contribution = " << contribution.head(10).transpose() << "...\n" ;
#endif

		scores += contribution ;
		counts += non_missingness ;
	}

	void RiskScoreComputation::compute( int sample, ResultCallback callback ) {
		for( Identifiers::const_iterator i = m_identifiers.begin(); i != m_identifiers.end(); ++i ) {
			compute_impl( sample, *i, callback ) ;
		}
	}

	void RiskScoreComputation::compute_impl( int sample, std::string const& risk_score_identifier, ResultCallback callback ) {
		Eigen::VectorXd const& scores = m_scores.at( risk_score_identifier ) ;
		Eigen::VectorXd const& counts = m_counts.at( risk_score_identifier ) ;
		assert( std::size_t( scores.size() ) == m_number_of_samples ) ;
		assert( std::size_t( counts.size() ) == m_number_of_samples ) ;
		callback( sample, risk_score_identifier + "_risk_score", scores( sample ) ) ;
		callback( sample, risk_score_identifier + "_risk_score_count", counts( sample ) ) ;
	}

	std::string RiskScoreComputation::get_summary( std::string const& prefix, std::size_t column_width ) const {
		std::ostringstream ostr ;
		ostr << prefix << "RiskScoreComputation\n"
			<< prefix << " - computing these risk scores:" ;
		std::size_t count = 0 ;
		for( std::set< std::string >::const_iterator i = m_identifiers.begin(); i != m_identifiers.end(); ++i, ++count ) {
			if( count > 10 ) {
				ostr << prefix << "     ..." ;
				break ;
			}
			ostr << prefix << "     " << *i << "(" << m_snp_counts.at(*i) << " SNPs)\n" ;
		}
		return ostr.str() ;
	}

	void RiskScoreComputation::add_effects( statfile::BuiltInTypeStatSource& source, ProgressCallback progress_callback ) {
		using namespace genfile::string_utils ;
		if(
			( source.number_of_columns() < 8)
			|| ( to_lower( source.name_of_column( 0 ) ) != "snpid" )
			|| ( to_lower( source.name_of_column( 1 ) ) != "rsid" )
			|| ( to_lower( source.name_of_column( 2 ) ) != "chromosome" )
			|| ( to_lower( source.name_of_column( 3 ) ) != "position" )
			|| ( to_lower( source.name_of_column( 4 ) ) != "allelea" )
			|| ( to_lower( source.name_of_column( 5 ) ) != "alleleb" )
		) {
			throw genfile::MalformedInputError(
				source.get_source_spec(),
				"Expected first six column names SNPID rsid chromosome position alleleA alleleB, but found: "
					+ source.name_of_column( 0 ) + " "
					+ source.name_of_column( 1 ) + " "
					+ source.name_of_column( 2 ) + " "
					+ source.name_of_column( 3 ) + " "
					+ source.name_of_column( 4 ) + " "
					+ source.name_of_column( 5 ),
				1
			) ;
		}
			
		if( !source.has_column( "risk_score_identifier" ) || !source.has_column( "additive_beta" ) || !source.has_column( "heterozygote_beta" ) ) {
			throw genfile::MalformedInputError( source.get_source_spec(), "Expected to find risk_score_identifier, additive_beta and heterozygote_beta columns.", 1 ) ;
		}
		
		std::size_t const identifier_column = source.index_of_column( "risk_score_identifier" ) ;
		std::size_t const additive_beta_column = source.index_of_column( "additive_beta" ) ;
		std::size_t const heterozygote_beta_column = source.index_of_column( "heterozygote_beta" ) ;
		std::size_t const min_beta_column = std::min( additive_beta_column, heterozygote_beta_column ) ;
		std::size_t const max_beta_column = std::max( additive_beta_column, heterozygote_beta_column ) ;
		bool additive_first = ( additive_beta_column < heterozygote_beta_column ) ;
		
		if( identifier_column > min_beta_column ) {
			throw genfile::MalformedInputError( source.get_source_spec(), "Expected to find the risk_score_identifier before the additive_beta and heterozygote_beta columns.", 1 ) ;
		}
		
		if( progress_callback ) {
			progress_callback( 0, source.number_of_rows() ) ;
		}
		
		std::set< std::string > identifiers ;
		
		genfile::SNPIdentifyingData snp ;
		std::string risk_score_identifier ;
		Eigen::MatrixXd betas( 3, 2 ) ;
		double beta1, beta2 ;
		while( source >> snp.SNPID() >> snp.rsid() >> snp.position().chromosome() >> snp.position().position() >> snp.first_allele() >> snp.second_allele() ) {
			source >> statfile::ignore( identifier_column - 6 ) >> risk_score_identifier ;
			source >> statfile::ignore( min_beta_column - identifier_column - 1 ) >> beta1 ;
			source >> statfile::ignore( max_beta_column - min_beta_column - 1 ) >> beta2 ;
			
			betas.setZero() ;
			if( additive_first ) {
				betas(1,0) = beta1 ;
				betas(2,0) = 2 * beta1 ;
				betas(1,1) = beta2 ;
			} else {
				betas(1,0) = beta2 ;
				betas(2,0) = 2 * beta2 ;
				betas(1,1) = beta1 ;
			}

			RiskScoreIdBetaMap::const_iterator where = m_map[ snp ].find( risk_score_identifier ) ;
			if( where != m_map[ snp ].end() ) {
				throw genfile::MalformedInputError( source.get_source_spec(), "SNP " + snp.get_rsid() + " has duplicated entries for risk score " + risk_score_identifier, source.number_of_rows_read() ) ;
			}

			m_map[ snp ][ risk_score_identifier ] = betas ;
			++m_snp_counts[ risk_score_identifier ] ;
			identifiers.insert( risk_score_identifier ) ;
			source >> statfile::end_row() ;
			if( progress_callback ) {
				progress_callback( source.number_of_rows_read(), source.number_of_rows() ) ;
			}
		}

		for( std::set< std::string >::const_iterator i = identifiers.begin(); i != identifiers.end(); ++i ) {
			m_scores[ *i ] = Eigen::VectorXd::Zero( m_number_of_samples ) ;
			m_counts[ *i ] = Eigen::VectorXd::Zero( m_number_of_samples ) ;
		}
		
		m_identifiers = identifiers ;

#if DEBUG_RISK_SCORE_COMPUTATION
		std::cerr << "I have these SNPs:\n" ;
		for( RiskScoreMap::const_iterator i = m_map.begin(); i != m_map.end(); ++i ) {
			std::cerr << i->first << ".\n" ;
		}
#endif
	}
}

