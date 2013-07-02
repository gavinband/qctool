
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
#include "genfile/RandomAccessSNPDataSource.hpp"
#include "genfile/SNPDataSourceChain.hpp"
#include "genfile/SampleFilteringSNPDataSource.hpp"
#include "genfile/wildcard.hpp"
#include "genfile/vcf/get_set.hpp"
#include "statfile/BuiltInTypeStatSink.hpp"
#include "integration/NewtonRaphson.hpp"
#include "integration/Derivative.hpp"
#include "db/Connection.hpp"
#include "db/SQLStatement.hpp"
#include "db/Error.hpp"
#include "components/HaplotypeFrequencyComponent/DBOutputter.hpp"
#include "components/HaplotypeFrequencyComponent/HaplotypeFrequencyComponent.hpp"

// #define DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT 1
struct HaplotypeFrequencyLogLikelihood ;

void HaplotypeFrequencyComponent::declare_options( appcontext::OptionProcessor& options ) {
	options.declare_group( "LD computation options" ) ;
	options[ "-compute-ld-with" ]
		.set_description( "Compute LD pairwise metrics between the main dataset and SNPs." )
		.set_takes_single_value() ;
	options[ "-max-ld-distance" ]
		.set_description( "Maximum physical distance between SNPs, above which LD will not be computed. "
			"A value of zero indicates LD between all SNPs will be computed. "
			"A plain number indicates distance in base pairs, or you add a Mb or kb suffix to specify the "
			"distance in megabases or kilobases if desired." )
		.set_takes_single_value()
		.set_default_value( "0" ) ;
	options.option_implies_option( "-compute-ld-with", "-o" ) ;
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
	genfile::SNPDataSource::UniquePtr source,
	appcontext::OptionProcessor const& options,
	appcontext::UIContext& ui_context
) {
	HaplotypeFrequencyComponent::UniquePtr result ;

#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
	std::cerr << "LD SNP samples: " << source->number_of_samples() << " (filtered from " << source->get_parent_source().number_of_samples() << ").\n" ;
#endif
	result.reset(
		new HaplotypeFrequencyComponent(
			source,
			ui_context
		)
	) ;
	
	result->set_max_distance( parse_physical_distance( options.get< std::string >( "-max-ld-distance" ))) ;

	haplotype_frequency_component::DBOutputter::SharedPtr outputter = haplotype_frequency_component::DBOutputter::create_shared(
		options.get_value< std::string >( "-o" ),
		options.get_value< std::string >( "-analysis-name" ),
		"HaplotypeFrequencyComponent",
		options.get_values_as_map()
	) ;

	result->send_results_to(
		boost::bind(
			&haplotype_frequency_component::DBOutputter::operator(),
			outputter,
			options.get< std::string >( "-analysis-name" ),
			_1,
			_2,
			_3,
			_4
		)
	) ;
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

void HaplotypeFrequencyComponent::begin_processing_snps( std::size_t number_of_samples ) {
	std::cerr << m_source->number_of_samples() << " : " <<  number_of_samples << ".\n" ;
	assert( m_source->number_of_samples() == number_of_samples ) ;
}

void HaplotypeFrequencyComponent::processed_snp( genfile::SNPIdentifyingData const& target_snp, genfile::VariantDataReader& target_data_reader ) {
	genfile::SNPIdentifyingData source_snp ;
	genfile::SingleSNPGenotypeProbabilities source_probs, target_probs ;
	target_data_reader.get( "genotypes", target_probs ) ;
	m_source->reset_to_start() ;
	while( m_source->get_snp_identifying_data( source_snp )) {
		if(
			( m_max_distance == 0 )
			||
			(
				( source_snp.get_position().chromosome() == target_snp.get_position().chromosome() )
				&&
				( std::abs( int64_t( source_snp.get_position().position() ) - int64_t( target_snp.get_position().position() ) ) <= m_max_distance )
			)
		) {
			genfile::VariantDataReader::UniquePtr source_data_reader = m_source->read_variant_data() ;
			compute_ld_measures( source_snp, *source_data_reader, target_snp, target_data_reader ) ;
		}
		else {
			m_source->ignore_snp_probability_data() ;
		}
	}
}

void HaplotypeFrequencyComponent::compute_ld_measures(
	genfile::SNPIdentifyingData const& source_snp,
	genfile::VariantDataReader& source_data_reader,
	genfile::SNPIdentifyingData const& target_snp,
	genfile::VariantDataReader& target_data_reader
) {
	std::vector< std::vector< int > > genotypes( 2 ) ;
	genfile::vcf::ThreshholdingGenotypeSetter< std::vector< int > > source_getter( genotypes[0], m_threshhold ) ;
	genfile::vcf::ThreshholdingGenotypeSetter< std::vector< int > > target_getter( genotypes[1], m_threshhold ) ;
	source_data_reader.get( "genotypes", source_getter ) ;
	target_data_reader.get( "genotypes", target_getter ) ;
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
	genfile::SNPIdentifyingData const& source_snp,
	genfile::SNPIdentifyingData const& target_snp,
	std::vector< std::vector< int > > const& genotypes
) {
	// Construct table of genotypes at each SNP.
	Eigen::Matrix3d table = Eigen::Matrix3d::Zero() ;
	for( std::size_t i = 0; i < m_source->number_of_samples(); ++i ) {
		if( genotypes[0][i] != -1 && genotypes[1][i] != -1 ) {
			++table( genotypes[0][i], genotypes[1][i] ) ;
		}
	}
	
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
	std::cerr << "SNP1: " << source_snp << "\n"
		<< "SNP2: " << target_snp << "\n"
		<< "table: " << table << ".\n" ;
#endif
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

		m_result_signal( source_snp, target_snp, "pi00", pi(0) ) ;
		m_result_signal( source_snp, target_snp, "pi01", pi(1) ) ;
		m_result_signal( source_snp, target_snp, "pi10", pi(2) ) ;
		m_result_signal( source_snp, target_snp, "pi11", pi(3) ) ;
		m_result_signal( source_snp, target_snp, "D", D ) ;
		m_result_signal( source_snp, target_snp, "Dprime", Dprime ) ;
		m_result_signal( source_snp, target_snp, "r", r ) ;
		m_result_signal( source_snp, target_snp, "r_squared", r * r ) ;
	}
	catch( genfile::OperationFailedError const& ) {
		m_ui_context.logger() << "!! Could not compute haplotype frequencies for " << source_snp << ", " << target_snp << ".\n"
			<< "!! table is:\n" << table << ".\n" ;
			
		genfile::VariantEntry const missing = genfile::MissingValue() ;
		m_result_signal( source_snp, target_snp, "pi00", missing ) ;
		m_result_signal( source_snp, target_snp, "pi01", missing ) ;
		m_result_signal( source_snp, target_snp, "pi10", missing ) ;
		m_result_signal( source_snp, target_snp, "pi11", missing ) ;
		m_result_signal( source_snp, target_snp, "D", missing ) ;
		m_result_signal( source_snp, target_snp, "Dprime", missing ) ;
		m_result_signal( source_snp, target_snp, "r", missing ) ;
		m_result_signal( source_snp, target_snp, "r_squared", missing ) ;
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
		m_result_signal( source_snp, target_snp, "dosage_r", dosage_r ) ;
		m_result_signal( source_snp, target_snp, "dosage_r_squared", dosage_r * dosage_r ) ;
	}
}

void HaplotypeFrequencyComponent::end_processing_snps() {}

void HaplotypeFrequencyComponent::send_results_to( ResultCallback callback ) {
	m_result_signal.connect( callback ) ;
}
