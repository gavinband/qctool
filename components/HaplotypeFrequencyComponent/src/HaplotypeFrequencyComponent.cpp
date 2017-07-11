
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread.hpp>
#include <Eigen/Core>
#include "genfile/SNPDataSourceProcessor.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/SampleFilteringSNPDataSource.hpp"
#include "genfile/SampleMappingSNPDataSource.hpp"
#include "genfile/wildcard.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "db/Error.hpp"
#include "components/HaplotypeFrequencyComponent/FlatTableDBOutputter.hpp"
#include "components/HaplotypeFrequencyComponent/HaplotypeFrequencyComponent.hpp"

// #define DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT 1
struct HaplotypeFrequencyLogLikelihood ;

void HaplotypeFrequencyComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "LD computation options" ) ;
	options[ "-compute-ld-with" ]
		.set_description( "Compute LD pairwise metrics between the main dataset and SNPs." )
		.set_takes_values(2) ;
	options[ "-max-ld-distance" ]
		.set_description( "Maximum physical distance between SNPs, above which LD will not be computed. "
			"A value of zero indicates LD between all SNPs will be computed. "
			"A plain number indicates distance in base pairs, or you add a Mb or kb suffix to specify the "
			"distance in megabases or kilobases if desired." )
		.set_takes_single_value()
		.set_default_value( "0" ) ;
	options[ "-min-r2" ]
		.set_description( "Minimum squared correlation between variants.  LD results for pairs of variants with "
			"lower than this squared correlation will not be output." )
		.set_takes_single_value()
		.set_default_value( "0" ) ;
	options.option_implies_option( "-compute-ld-with", "-osnp" ) ;
}

namespace {
	uint64_t parse_physical_distance( std::string distance ) {
		std::size_t number_part_length = std::string::npos ;
		uint64_t multiplier = 1 ;
		// allowable suffixes are mb, kb, bp.
		
		if( distance.size() == 0 ) {
			throw genfile::BadArgumentError( "parse_physical_distance()", "distance=\"" + distance + "\"" ) ;
		}
		
		if( distance[ distance.size() - 1 ] != 'b' && distance[ distance.size() - 1 ] != 'b' ) {
			// no suffix.
			number_part_length = distance.size() ;
			multiplier = 1 ;
		}
		else {
			if( distance.size() < 3 ) {
				throw genfile::BadArgumentError( "parse_physical_distance()", "distance=\"" + distance + "\"" ) ;
			}
			number_part_length = distance.size() - 2 ;
			std::string suffix = genfile::string_utils::to_lower( distance.substr( distance.size() - 2, 2 ) ) ;
			if( suffix == "bp" ) {
				multiplier = 1 ;
			} else if( suffix == "kb" ) {
				multiplier = 1000 ;
			} else if( suffix == "mb" ) {
				multiplier = 1000000 ;
			} else {
				throw genfile::BadArgumentError( "parse_physical_distance()", "distance=\"" + distance + "\"" ) ;
			}
		}
		
		uint64_t const result = genfile::string_utils::to_repr< double >( distance.substr( 0, number_part_length ) ) * multiplier ;
		std::cerr << "Parsed distance \"" + distance + "\" as " + genfile::string_utils::to_string( result ) + " base pairs.\n" ;
		return result ;
	}
}

HaplotypeFrequencyComponent::UniquePtr HaplotypeFrequencyComponent::create(
	genfile::CohortIndividualSource const& source_samples,
	std::string const& source_sample_id_column,
	genfile::CohortIndividualSource::UniquePtr samples,
	std::string const& ld_sample_id_column,
	genfile::SNPDataSource::UniquePtr source,
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
) {
	source.reset(
		new genfile::SampleMappingSNPDataSource(
			source_samples,
			source_sample_id_column,
			*samples,
			ld_sample_id_column,
			source
		)
	) ;
	
	HaplotypeFrequencyComponent::UniquePtr result ;

#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
	std::cerr << "LD SNP samples: " << source->number_of_samples() << " (matched from " << source->get_parent_source().number_of_samples() << " in raw data).\n" ;
#endif
	result.reset(
		new HaplotypeFrequencyComponent(
			source,
			ui_context
		)
	) ;
	
	result->set_max_distance( parse_physical_distance( options.get< std::string >( "-max-ld-distance" ))) ;
	// result->set_min_r2( parse_physical_distance( options.get< double >( "-min-r2" ))) ;

	haplotype_frequency_component::FlatTableDBOutputter::UniquePtr outputter
		= haplotype_frequency_component::FlatTableDBOutputter::create(
			options.get_value< std::string >( "-osnp" ),
			options.get_value< std::string >( "-analysis-name" ),
			"HaplotypeFrequencyComponent",
			options.get_values_as_map()
	) ;

	result->send_results_to( outputter ) ;
	return result ;
}

HaplotypeFrequencyComponent::HaplotypeFrequencyComponent(
	genfile::SNPDataSource::UniquePtr source,
	appcontext::UIContext& ui_context
):
	m_source( source ),
	m_ui_context( ui_context ),
	m_threshhold( 0.9 ),
	m_max_distance( 200000 )
{}

void HaplotypeFrequencyComponent::set_max_distance( uint64_t distance ) {
	m_max_distance = distance ;
}

void HaplotypeFrequencyComponent::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) {
	std::cerr << m_source->number_of_samples() << " : " <<  number_of_samples << ".\n" ;
	assert( m_source->number_of_samples() == number_of_samples ) ;
}

void HaplotypeFrequencyComponent::processed_snp(
	genfile::VariantIdentifyingData const& target_snp,
	genfile::VariantDataReader& target_data_reader
) {
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
	std::cerr << "Processing " << target_snp << "...\n" ;
#endif
	genfile::VariantIdentifyingData source_snp ;
	m_source->reset_to_start() ;
	while( m_source->get_snp_identifying_data( &source_snp )) {
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
		std::cerr << "Comparing to " << source_snp << "...\n" ;
#endif
		if(
			( m_max_distance == 0 )
			||
			(
				( source_snp.get_position().chromosome() == target_snp.get_position().chromosome() )
				&&
				( std::abs( int64_t( source_snp.get_position().position() ) - int64_t( target_snp.get_position().position() ) ) <= m_max_distance )
			)
		) {
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
			std::cerr << "Computing LD measures for " << source_snp << " : " << target_snp << "...\n" ;
#endif
			genfile::VariantDataReader::UniquePtr source_data_reader = m_source->read_variant_data() ;
			compute_ld_measures( source_snp, *source_data_reader, target_snp, target_data_reader ) ;
		}
		else {
			m_source->ignore_snp_probability_data() ;
		}
	}
}

#if 0
namespace {
	struct PairwiseGetter {
		PairwiseGetter() {}
		
		struct Getter {
			Getter( std::vector< int >* result ):
				m_result( result )
			{
				// missing
				m_coding.push_back( -1 ) ; // missing

				// haploid
				m_coding.push_back( 0 ) ; // 00
				m_coding.push_back( 1 ) ; // 01
				m_coding.push_back( 2 ) ; // 10
				m_coding.push_back( 3 ) ; // 11
				
				// diploid
				m_coding.push_back( 2 ) ;
				m_coding.push_back( 3 ) ;
				m_coding.push_back( 4 ) ;
			}
				
			void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
				assert( number_of_alleles == 2 ) ;
				m_result->resize( number_of_samples ) ;
				m_sample_i = 0 ;
			}
			
			void set_min_max_ploidy(
				std::size_t min_ploidy,
				std::size_t max_ploidy,
				std::size_t min_entries,
				std::size_t max_entries
			) {
				assert( min_ploidy > 0 ) ;
				assert( max_ploidy < 3 ) ;
			}
			
			bool set_sample( std::size_t i ) {
				m_sample_i = i ;
				return true ;
			}
			
			void set_number_of_entries(
				std::size_t ploidy,
				std::size_t entries,
				genfile::OrderType order_type,
				genfile::ValueType value_type
			) {
				assert( m_ploidy == ploidy )
				assert( value_type == eAlleleIndex ) ;
				assert( order_type == ePerOrderedHaplotype || order_type == ePerUnorderedGenotype ) ;
				m_order_type = order_type ;
				m_ploidy = ploidy ;

				// store genotypes as bit representation.
				// bit 4 means 'missing'.
				// bit 3 means 'unphased'.
				// bit 0 means haploid genotype
				// bits 0 means haploid genotype
				// bits 1-2 mean diploid genotype
				(*m_result)[ m_sample_i ] = ((m_ploidy > 1) && (order_type == ePerUnorderedGenotype)) ? (1 << 3) : 0;
			}
			
			void set_value( std::size_t entry_i, int value ) {
				assert( value == 0 || value == 1 ) ;
				if( m_order_type == ePerOrderedHaplotype ) {
					(*m_result)[ m_sample_i ] |= ( value << ( entry_i + (ploidy-1))) ;
				} else {
					(*m_result)[ m_sample_i ] += ( value << (ploidy-1)) ;
				}
			}

			void set_value( std::size_t entry_i, genfile::MissingValue value ) {
				(*m_result)[ m_sample_i ] = 0x8 ;
			}
			
			void finalise() {
				// convert wrong-way-round diploid genotypes, e.g. 1/0 to right-way-round, e.g. 0/1
				for( std::size_t i = 0; i < m_result->size(); ++i ) {
					
				}
			}

		private:
			std::vector< int > m_coding ;
			std::size_t m_sample_i ;
			std::vector< int >* m_result ;
		} ;
			
			
		void build_tables( Eigen::MatrixXd* haploid_table, Eigen::MatrixXd* diploid_table ) {
			haploid_table->setZero( 2, 2 ) ;
			diploid_table->setZero( 3, 3 ) ;
			for( std::size_t i = 0; i < m_genotypes[0].size(); ++i ) {
				int const& genotype1 = m_genotypes[0][i] ;
				int const& genotype2 = genotype2 ;
				if( !is_missing( genotype1 ) & !is_missing( genotype2 )) {
					assert( ploidy( genotype1 ) == ploidy( genotype2 )) ;
					if( is_phased( m_genotypes[0][i] ) & is_phased( genotype2 )) {
						(*haploid_table)( genotype1 & 1, genotype1 & 1 ) += 1 ;
						if( ploidy( genotype1 ) == 2 ) {
							(*haploid_table)( genotype1 & 2, genotype1 & 2 ) += 1 ;
						}
					} else {
						(*diploid_table)( genotype1 & 0x3, genotype2 & 0x3 ) += 1 ;
					}
				}
			}
		}

	private:
			std::vector< std::vector< int > > m_genotypes ;

			bool haploid( int genotype ) ;
	} ;
}
#endif

void HaplotypeFrequencyComponent::compute_ld_measures(
	genfile::VariantIdentifyingData const& source_snp,
	genfile::VariantDataReader& source_data_reader,
	genfile::VariantIdentifyingData const& target_snp,
	genfile::VariantDataReader& target_data_reader
) {
	std::vector< Eigen::VectorXd > genotypes( 2 ) ;
	source_data_reader.get( ":genotypes:", genfile::vcf::get_threshholded_calls< Eigen::VectorXd >( genotypes[0], m_threshhold ) ) ;
	target_data_reader.get( ":genotypes:", genfile::vcf::get_threshholded_calls< Eigen::VectorXd >( genotypes[1], m_threshhold ) ) ;
	assert( genotypes[0].size() == m_source->number_of_samples() ) ;
	assert( genotypes[0].size() == genotypes[1].size() ) ;
	try {
		compute_ld_measures(
			source_snp,
			target_snp,
			genotypes
		) ;
	} catch( genfile::OperationFailedError const& e ) {
		m_ui_context.logger() << "!! HaplotypeFrequencyComponent::compute_ld_measures(): could not compute LD measures between SNPs "
			<< source_snp << " and " << target_snp << ".\n"
			<< "!! reason: " << e.get_message() << "\n" ;
	}
}

void HaplotypeFrequencyComponent::compute_ld_measures(
	genfile::VariantIdentifyingData const& source_snp,
	genfile::VariantIdentifyingData const& target_snp,
	std::vector< Eigen::VectorXd > const& genotypes
) {
	// Construct table of genotypes at each SNP.
	Eigen::Matrix3d table = Eigen::Matrix3d::Zero() ;
	for( std::size_t i = 0; i < m_source->number_of_samples(); ++i ) {
		if( genotypes[0](i) != -1 && genotypes[1](i) != -1 ) {
			++table( genotypes[0](i), genotypes[1](i) ) ;
		}
	}
	
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
	std::cerr << "SNP1: " << source_snp << "\n"
		<< "SNP2: " << target_snp << "\n"
		<< "table: " << table << ".\n" ;
#endif
	genfile::VariantEntry const missing = genfile::MissingValue() ;
	if( table.array().maxCoeff() > 0 ) {
		try {
			HaplotypeFrequencyLogLikelihood ll( table ) ;
			HaplotypeFrequencyLogLikelihood::Vector pi( 4 ) ;
			pi.tail( 3 ) = ll.get_MLE_by_EM() ;
			pi(0) = 1.0 - pi.tail( 3 ).sum() ;
			double D = pi(0) * pi(3) - pi(1) * pi(2) ;
			double max_D ;
			if( D < 0 ) {
				max_D = std::min( (pi(0) + pi(1)) * (pi(0)+pi(2)), (pi(1)+pi(3))*(pi(2)+pi(3)) ) ;
			}
			else {
				max_D = std::min( (pi(0) + pi(1)) * (pi(1)+pi(3)), (pi(0)+pi(2))*(pi(2)+pi(3)) ) ;
			}
			double Dprime = D / max_D ;
			double r = D / std::sqrt( (pi(0)+pi(1)) * (pi(2)+pi(3)) * (pi(0)+pi(2)) * (pi(1)+pi(3))) ;

			m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi00", pi(0) ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi01", pi(1) ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi10", pi(2) ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi11", pi(3) ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "D", D ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "Dprime", Dprime ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "r", r ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "r_squared", r * r ) ;
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
		std::cerr << "Output compltee.\n" ;
#endif
		}
		catch( genfile::OperationFailedError const& ) {
			m_ui_context.logger() << "!! Could not compute haplotype frequencies for " << source_snp << ", " << target_snp << ".\n"
				<< "!! table is:\n" << table << ".\n" ;
			
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi00", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi01", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi10", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi11", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "D", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "Dprime", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "r", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "r_squared", missing ) ;
		}
	
		{
			// compute allele dosage r squared too.
			Eigen::Vector2d means = Eigen::Vector2d::Zero() ;
			for( int i = 0; i < 3; ++i ) {
				means(0) += table.row( i ).sum() * i ;
				means(1) += table.col( i ).sum() * i ;
			}
		
			means /= table.sum() ;
	#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT	
			std::cerr << "table: " << table << "\n" ;
			std::cerr << "means:\n" << std::fixed << std::setprecision( 5 ) <<  means << ".\n" ;
	#endif
			double dosage_cov = 0.0 ;
			Eigen::Vector2d variances = Eigen::Vector2d::Zero() ;
			for( int i = 0; i < 3; ++i ) {
				for( int j = 0; j < 3; ++j ) {
					dosage_cov += table(i,j)  * ( i - means(0) ) * ( j - means(1) ) ;
				}
				variances(0) += table.row( i ).sum() * ( i - means(0) ) * ( i - means(0) ) ;
				variances(1) += table.col( i ).sum() * ( i - means(1) ) * ( i - means(1) ) ;
			}
		
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT	
			std::cerr << "dosage_cov:\n" << dosage_cov << ".\n" ;
			std::cerr << "variances:\n" << variances << ".\n" ;
#endif
		
			double dosage_r = dosage_cov / std::sqrt( variances(0) * variances( 1 ) ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "dosage_r", dosage_r ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, "dosage_r_squared", dosage_r * dosage_r ) ;
		}
	} else {
		m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi00", missing ) ;
		m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi01", missing ) ;
		m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi10", missing ) ;
		m_sink->store_per_variant_pair_data( source_snp, target_snp, "pi11", missing ) ;
		m_sink->store_per_variant_pair_data( source_snp, target_snp, "D", missing ) ;
		m_sink->store_per_variant_pair_data( source_snp, target_snp, "Dprime", missing ) ;
		m_sink->store_per_variant_pair_data( source_snp, target_snp, "r", missing ) ;
		m_sink->store_per_variant_pair_data( source_snp, target_snp, "r_squared", missing ) ;
		m_sink->store_per_variant_pair_data( source_snp, target_snp, "dosage_r", missing ) ;
		m_sink->store_per_variant_pair_data( source_snp, target_snp, "dosage_r_squared", missing ) ;
	}
}

void HaplotypeFrequencyComponent::end_processing_snps() {
	if( m_sink.get() ) {
		m_sink->finalise() ;
	}
}

void HaplotypeFrequencyComponent::send_results_to( haplotype_frequency_component::FlatTableDBOutputter::UniquePtr sink ) {
	m_sink = sink ;
}
