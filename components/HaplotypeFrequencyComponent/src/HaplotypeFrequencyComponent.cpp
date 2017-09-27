
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
#include "genfile/ToGP.hpp"
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
	options[ "-old" ]
		.set_description( "Specify file to write LD metrics to" )
		.set_takes_single_value() ;
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
	options.option_implies_option( "-compute-ld-with", "-old" ) ;
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

namespace {
	genfile::SampleStratification compute_stratification( genfile::CohortIndividualSource const& samples, std::string const& variable ) {
		genfile::SampleStratification result ;
		genfile::CohortIndividualSource::ColumnSpec const spec = samples.get_column_spec() ;
		if( !spec[ variable ].is_discrete() ) {
			throw genfile::BadArgumentError( "void impl::compute_strata()", "variable=\"" + variable + "\"" ) ;
		}

		for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
			genfile::VariantEntry const& entry = samples.get_entry( i, variable ) ;
			if( !entry.is_missing() ) {
				result.add_sample( variable + "=" + entry.as< std::string >(), i ) ;
			}
		}

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
	result->set_min_r2( options.get< double >( "-min-r2" )) ;

	using namespace genfile::string_utils ;
	std::string filename = options.get_value< std::string >( "-old" ) ;
	if( filename.size() >= 9 && filename.substr( 0, 9 ) == "sqlite://" ) {
		filename = filename.substr( 9, filename.size() ) ;
	}
	std::vector< std::string > outputFilenameElts = genfile::string_utils::split( filename, ":" ) ;
	assert( outputFilenameElts.size() > 0 ) ;
	if( outputFilenameElts.size() > 2 ) {
		throw genfile::BadArgumentError(
			"HaplotypeFrequencyComponent::create()",
			"-old \"" + options.get_value< std::string >( "-old" ) + "\"",
			"Value should be of the form [sqlite://]<filename>[:<tablename>]"
		) ;
	}
	
	haplotype_frequency_component::FlatTableDBOutputter::UniquePtr outputter
		= haplotype_frequency_component::FlatTableDBOutputter::create(
			outputFilenameElts[0],
			options.get_value< std::string >( "-analysis-name" ),
			"HaplotypeFrequencyComponent",
			options.get_values_as_map()
	) ;
	if( outputFilenameElts.size() == 2 ) {
		outputter->set_table_name( outputFilenameElts[1] ) ;
	}
	
	if( options.check( "-stratify" )) {
		result->set_stratification( compute_stratification( source_samples, options.get< std::string >( "-stratify" )) ) ;
	}

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
	m_max_distance( 200000 ),
	m_min_r2( 0.0 )
{}

void HaplotypeFrequencyComponent::set_max_distance( uint64_t distance ) {
	m_max_distance = distance ;
}

void HaplotypeFrequencyComponent::set_min_r2( double min_r2 ) {
	assert( min_r2 >= 0.0 ) ;
	m_min_r2 = min_r2 ;
}

void HaplotypeFrequencyComponent::set_stratification( genfile::SampleStratification stratification ) {
	m_stratification = stratification ;
}

void HaplotypeFrequencyComponent::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) {
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

namespace {
	struct ThreshholdCalls: public genfile::VariantDataReader::PerSampleSetter {
		ThreshholdCalls( Eigen::VectorXd* data, double const threshhold ):
			m_data( data ),
			m_missing_value( -1.0 ),
			m_threshhold( threshhold )
		{}

		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
			assert( nAlleles == 2 ) ;
			m_data->resize( nSamples ) ;
			m_data->setConstant( nSamples, m_missing_value ) ;
		}

		bool set_sample( std::size_t n ) {
			m_sample_i = n ;
			return true ;
		}

		void set_number_of_entries( uint32_t ploidy, std::size_t n, genfile::OrderType const, genfile::ValueType const ) {
			assert( ploidy == 2 ) ;
			assert( n == 3 ) ;
		}

		void set_value( std::size_t entry_i, genfile::MissingValue const value ) {
			(*m_data)[m_sample_i] = m_missing_value ;
		}

		void set_value( std::size_t entry_i, double const value ) {
			assert( m_value >= 0.0 ) ;
			//std::cerr << m_sample_i << ": " << entry_i << ": " << value << "\n" ;
			bool passesThreshhold = ( value >= m_threshhold ) ;
			if( passesThreshhold ) {
				// should only get one value passing threshhold.
				//assert( (*m_data)[m_sample_i] == m_missing_value ) ;
				(*m_data)[m_sample_i] = entry_i ;
			}
		}

		void finalise() {
			// nothing to do
		}
		
		
	private:
		Eigen::VectorXd* const m_data ;
		double const m_missing_value ;
		double const m_threshhold ;
		std::size_t m_sample_i ;
		double m_value ;
	} ;
}


void HaplotypeFrequencyComponent::compute_ld_measures(
	genfile::VariantIdentifyingData const& source_snp,
	genfile::VariantDataReader& source_data_reader,
	genfile::VariantIdentifyingData const& target_snp,
	genfile::VariantDataReader& target_data_reader
) {
	if( source_snp.number_of_alleles() == 2 && target_snp.number_of_alleles() == 2 ) {
		std::vector< Eigen::VectorXd > genotypes( 2 ) ;
		ThreshholdCalls setter1( &(genotypes[0] ), m_threshhold ) ;
		ThreshholdCalls setter2( &(genotypes[1] ), m_threshhold ) ;
		source_data_reader.get( ":genotypes:", genfile::to_GP_unphased( setter1 ) ) ;
		target_data_reader.get( ":genotypes:", genfile::to_GP_unphased( setter2 ) ) ;
	//	source_data_reader.get( ":genotypes:", setter1 ) ;
	//	target_data_reader.get( ":genotypes:", setter2 ) ;
		assert( genotypes[0].size() == m_source->number_of_samples() ) ;
		assert( genotypes[0].size() == genotypes[1].size() ) ;

		compute_ld_measures(
			source_snp,
			genotypes[0],
			target_snp,
			genotypes[1]
		) ;
	} else {
		// can't compute LD for multiallelics right now
	}
}

void HaplotypeFrequencyComponent::compute_ld_measures(
	genfile::VariantIdentifyingData const& source_snp,
	Eigen::VectorXd const& source_genotypes,
	genfile::VariantIdentifyingData const& target_snp,
	Eigen::VectorXd const& target_genotypes
) {
	if( m_stratification ) {
		for( std::size_t i = 0; i < m_stratification->number_of_strata(); ++i ) {
			try {
				compute_ld_measures(
					source_snp,
					source_genotypes,
					target_snp,
					target_genotypes,
					m_stratification->stratum_name(i) + ":",
					m_stratification->stratum(i)
				) ;
			} catch( genfile::OperationFailedError const& e ) {
				m_ui_context.logger() << "!! HaplotypeFrequencyComponent::compute_ld_measures(): "
					<< "stratum " << m_stratification->stratum_name(i) << ": "
					<< "could not compute LD measures between SNPs "
					<< source_snp << " and " << target_snp << ".\n"
					<< "!! reason: " << e.get_message() << "\n" ;
			}
		}
	} else {
		try {
			compute_ld_measures(
				source_snp,
				source_genotypes,
				target_snp,
				target_genotypes,
				"",
				std::vector< genfile::SampleRange  >( 1, genfile::SampleRange( 0, source_genotypes.size() ) )
			) ;
		} catch( genfile::OperationFailedError const& e ) {
			m_ui_context.logger() << "!! HaplotypeFrequencyComponent::compute_ld_measures(): "
				<< "could not compute LD measures between SNPs "
				<< source_snp << " and " << target_snp << ".\n"
				<< "!! reason: " << e.get_message() << "\n" ;
		}
	}
}

void HaplotypeFrequencyComponent::compute_ld_measures(
	genfile::VariantIdentifyingData const& source_snp,
	Eigen::VectorXd const& source_genotypes,
	genfile::VariantIdentifyingData const& target_snp,
	Eigen::VectorXd const& target_genotypes,
	std::string const& variable_name_stub,
	std::vector< genfile::SampleRange > const& sample_set
) {
	// Construct table of genotypes at each SNP.
	Eigen::Matrix3d table = Eigen::Matrix3d::Zero() ;
	for( std::size_t i = 0; i < sample_set.size(); ++i ) {
		for( std::size_t j = sample_set[i].begin(); j < sample_set[i].end(); ++j ) {
			assert( source_genotypes(j) > -2 && source_genotypes(j) < 3 ) ;
			assert( target_genotypes(j) > -2 && target_genotypes(j) < 3 ) ;
			if( source_genotypes(j) != -1 && target_genotypes(j) != -1 ) {
				++table( source_genotypes(j), target_genotypes(j) ) ;
			}
		}
	}
	
	// We estimate LD under a prior, represented as data augmentation 
	// that assumes we have previously seen at least one of each haplotype.
	// That corresponds to the assumptions:
	// 1: that both variants are polymorphic, and
	// 2: there is at least some recombination (or difference) between them.
	Eigen::Matrix2d prior ;
	prior
		<< 0, 0,
		   0, 0
	;
	
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
	std::cerr << "SNP1: " << source_snp << "\n"
		<< "SNP2: " << target_snp << "\n"
		<< "table: " << table << ".\n" ;
#endif
	
	genfile::VariantEntry::Integer const N = table.sum() ;
	bool includeInOutput = (m_min_r2 == 0.0 ) ;
	genfile::VariantEntry const missing = genfile::MissingValue() ;

	if( includeInOutput ) {
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "number_of_genotypes", N ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "number_of_haplotypes", missing ) ;
	}

	if( N == 0 ) {
		if( includeInOutput ) {
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi00", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi01", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi10", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi11", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "D", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "Dprime", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "r", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "r2", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "dosage_r", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "dosage_r2", missing ) ;
		}
	} else {
		try {
			HaplotypeFrequencyLogLikelihood ll( table, prior ) ;
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
			double const Dprime = D / max_D ;
			double const r = D / std::sqrt( (pi(0)+pi(1)) * (pi(2)+pi(3)) * (pi(0)+pi(2)) * (pi(1)+pi(3))) ;
			double const r2 = r*r ;
			includeInOutput = (r2 >= m_min_r2 ) ;
			if( includeInOutput ) {
				m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi00", pi(0) ) ;
				m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi01", pi(1) ) ;
				m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi10", pi(2) ) ;
				m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi11", pi(3) ) ;
				m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "D", D ) ;
				m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "Dprime", Dprime ) ;
				m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "r", r ) ;
				m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "r2", r2 ) ;
			}
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
		std::cerr << "Output complete.\n" ;
#endif
		}
		catch( genfile::OperationFailedError const& ) {
			m_ui_context.logger() << "!! Could not compute haplotype frequencies for " << source_snp << ", " << target_snp << ".\n"
				<< "!! table is:\n" << table << ".\n" ;

			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi00", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi01", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi10", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi11", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "D", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "Dprime", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "r", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "r2", missing ) ;
		}
	
		if( includeInOutput ) {
			// compute allele dosage r squared too.
			Eigen::Vector2d means = Eigen::Vector2d::Zero() ;
			for( int i = 0; i < 3; ++i ) {
				means(0) += table.row( i ).sum() * i ;
				means(1) += table.col( i ).sum() * i ;
			}
			means /= table.sum() ;
	#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT	
			std::cerr << "table: " << table << "\n" ;
			std::cerr << "prior: " << prior << "\n" ;
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
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "dosage_r", dosage_r ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "dosage_r2", dosage_r * dosage_r ) ;
		}
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
