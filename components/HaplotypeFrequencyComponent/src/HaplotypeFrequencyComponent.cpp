
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
	options[ "-prior-ld-weight" ]
		.set_description( "Weight w to place on shrinkage prior in computation of pairwise LD. "
			"This is interpreted as adding dummy observations of w/4 for each of the four possible haplotypes, "
			"or equivalently, to placing a Dirichlet(w/4, w/4, w/4, w/4) prior on the vector of haplotype frequencies." )
		.set_takes_single_value()
		.set_default_value( 1.0 ) ;
	options.option_implies_option( "-compute-ld-with", "-old" ) ;
	options.option_implies_option( "-prior-ld-weight", "-compute-ld-with" ) ;
	options.option_implies_option( "-max-ld-distance", "-compute-ld-with" ) ;
	options.option_implies_option( "-min-r2", "-compute-ld-with" ) ;
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
	{
		double const w = options.get< double >( "-prior-ld-weight" ) ;
		Eigen::MatrixXd prior(2,2) ;
		prior
			<< w/4, w/4,
			   w/4, w/4
		;
		result->set_prior( prior ) ;
	}

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
	m_min_r2( 0.0 ),
	m_prior( Eigen::Matrix2d::Zero() )
{
}

void HaplotypeFrequencyComponent::set_max_distance( uint64_t distance ) {
	m_max_distance = distance ;
}

void HaplotypeFrequencyComponent::set_min_r2( double min_r2 ) {
	assert( min_r2 >= 0.0 ) ;
	m_min_r2 = min_r2 ;
}

void HaplotypeFrequencyComponent::set_prior( Eigen::Matrix2d const& matrix ) {
	assert( matrix.array().minCoeff() >= 0.0 ) ;
	m_prior = matrix ;
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
	struct CallSetter: public genfile::VariantDataReader::PerSampleSetter {
		enum { ePhased = 0x80000000 } ;

		CallSetter(
			std::vector< int >* result,
			std::vector< uint32_t >* ploidy
		):
			m_result( result ),
			m_ploidy( ploidy )
		{
			assert( result ) ;
		}

		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
			if( nAlleles != 2 ) {
				throw genfile::BadArgumentError(
					"CallSetter::initialise()",
					"nAlleles=" + genfile::string_utils::to_string( nAlleles ),
					"I only support biallelic variants"
				) ;
			}
			m_result->clear() ; // ploidy 2
			m_result->resize( nSamples, -1 ) ;
			m_ploidy->clear() ;
			m_ploidy->resize( nSamples, 0 ) ;
		}

		bool set_sample( std::size_t n ) {
			m_sample_i = n ;
			return true ;
		}

		void set_number_of_entries(
			uint32_t ploidy, std::size_t n,
			genfile::OrderType const order_type,
			genfile::ValueType const value_type
		) {
			if( ploidy > 2 ) {
				throw genfile::BadArgumentError(
					"Callsetter::set_number_of_entries()",
					"ploidy=" + genfile::string_utils::to_string( ploidy ),
					"Only haploid or diploid samples are currently supported."
				) ;
			}
			if( value_type != genfile::eAlleleIndex ) {
				throw genfile::BadArgumentError(
					"Callsetter::set_number_of_entries()",
					"value_type=" + genfile::string_utils::to_string( value_type ),
					"Expected a hard-called genotype, consider using threshholded calls."
				) ;
			}
			if( order_type != genfile::ePerOrderedHaplotype && order_type != genfile::ePerUnorderedHaplotype ) {
				throw genfile::BadArgumentError(
					"Callsetter::set_number_of_entries()",
					"order_type=" + genfile::string_utils::to_string( order_type ),
					"Expected a hard-called genotype, consider using threshholded calls."
				) ;
			}
			assert( value_type == genfile::eAlleleIndex ) ;
			m_order_type = order_type ;
			uint32_t storedPloidy = ploidy ;
			if( m_order_type == genfile::ePerOrderedHaplotype ) {
				storedPloidy |= ePhased ;
			}
			(*m_ploidy)[m_sample_i] = storedPloidy ;
			(*m_result)[m_sample_i] = 0 ;
		}

		void set_value( std::size_t entry_i, genfile::MissingValue const value ) {
			(*m_result)[m_sample_i] = -1; 
		}

		void set_value( std::size_t entry_i, Integer const value ) {
			// Only accumulate if not missing
			int& stored = (*m_result)[m_sample_i] ;
			if( stored != -1 ) {
				if( m_order_type == genfile::ePerUnorderedHaplotype ) {
					// Compute count of 2nd allele
					(*m_result)[m_sample_i] += value ;
				}
				else if( m_order_type == genfile::ePerOrderedHaplotype ) {
					// Put haplotype allele counts in seperate bits.
					(*m_result)[m_sample_i] += (value << entry_i) ;
				}
				else {
					assert(0) ;
				}
			}
		}

		void set_value( std::size_t entry_i, double const value ) {
			assert(0) ; // expecting GT field
		}

		void finalise() {
			// nothing to do
		}
		
	private:
		std::vector< int >* m_result ;
		std::vector< uint32_t >* m_ploidy ;
		std::size_t m_sample_i ;
		genfile::OrderType m_order_type ;
	} ;
	
	void tabulate_calls(
		std::vector< genfile::SampleRange > const& sample_set,
		std::vector< int > left_calls,
		std::vector< uint32_t > left_ploidy,
		std::vector< int > right_calls,
		std::vector< uint32_t > right_ploidy,
		Eigen::Matrix3d* diploid_table,
		Eigen::Matrix2d* haploid_table
	) {
		for( std::size_t range_i = 0; range_i < sample_set.size(); ++range_i ) {
			for( std::size_t j = sample_set[range_i].begin(); j < sample_set[range_i].end(); ++j ) {
				if( right_ploidy[j] == left_ploidy[j] && left_calls[j] != -1 && right_calls[j] != -1 ) {
					uint32_t const ploidy = left_ploidy[j] & 0xF;
					uint32_t const phased = left_ploidy[j] & CallSetter::ePhased ;
					switch( ploidy ) {
						case 0: break; // zeroploid, nothing to do
						case 1:
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT > 1
							std::cerr << "Setting haploid: " << left_calls[j] << ", " << right_calls[j] << ".\n" ;
#endif
							++((*haploid_table)( left_calls[j], right_calls[j] )) ;
							break ;
						case 2:
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT > 1
							std::cerr << "Setting diploid: " << left_calls[j] << ", " << right_calls[j] << ".\n" ;
#endif
							if( phased != 0 ) {
								++(*haploid_table)(
									(left_calls[j] & 0x1),
									(right_calls[j] & 0x1)
								) ;
								++(*haploid_table)(
									(left_calls[j] & 0x2) >> 1,
									(right_calls[j] & 0x2) >> 1
								) ;
							} else {
								++(*diploid_table)( left_calls[j], right_calls[j] ) ;
							}
							break ;
						default:
							assert(0) ;
					}
				}
			}
		}
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
		std::cerr << "Haploid table:\n" << *haploid_table << "\n" ;
		std::cerr << "Diploid table:\n" << *diploid_table << "\n" ;
#endif
	} ;

	struct DosageSetter: public genfile::VariantDataReader::PerSampleSetter {
		DosageSetter(
			Eigen::MatrixXd* result,
			Eigen::MatrixXd* nonmissingness,
			int column
		):
			m_result( result ),
			m_nonmissingness( nonmissingness ),
			m_order_type( genfile::eUnknownOrderType ),
			m_result_column( column )
		{
			assert( result ) ;
			assert( column >= 0 && column < result->cols() ) ;
		}

		void initialise( std::size_t nSamples, std::size_t nAlleles ) {
			assert( nAlleles == 2 ) ;
			assert( m_result->rows() == nSamples ) ;
			assert( m_nonmissingness->rows() == nSamples ) ;
		}

		bool set_sample( std::size_t n ) {
			m_sample_i = n ;
			m_order_type = genfile::eUnknownOrderType ;
			return true ;
		}

		void set_number_of_entries(
			uint32_t ploidy, std::size_t n,
			genfile::OrderType const order_type,
			genfile::ValueType const value_type
		) {
			assert( ploidy <= 2 ) ;
			assert( (order_type == genfile::ePerOrderedHaplotype || genfile::ePerUnorderedGenotype) && value_type == genfile::eProbability ) ;
			m_order_type = order_type ;
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT > 1
			std::cerr << m_order_type << ", " << m_sample_i << ", " << m_result_column << "\n" << std::flush ;
			std::cerr << m_result->rows() << " x " << m_result->cols() << ".\n" << std::flush ;
			std::cerr << m_nonmissingness->rows() << " x " << m_nonmissingness->cols() << ".\n" << std::flush ;
#endif
			(*m_result)(m_sample_i,m_result_column) = 0 ;
			(*m_nonmissingness)(m_sample_i,m_result_column) = 0 ;
		}

		void set_value( std::size_t entry_i, genfile::MissingValue const value ) {
			(*m_result)(m_sample_i, m_result_column) = 0 ;
			(*m_nonmissingness)(m_sample_i,m_result_column) = 0 ;
		}

		void set_value( std::size_t entry_i, double const value ) {
			// For biallelic variants, genotypes come in the order of
			// the number of B alleles, so this computes dosage:
			if( m_order_type == genfile::ePerOrderedHaplotype ) {
				// genotypes come in the order A, B (1st hap); A, B (2nd hap)
				// assumption is variant is biallelic.
				(*m_result)(m_sample_i, m_result_column) += (entry_i % 2) * value ;
			} else {
				// order type = genfile::ePerUnorderedGenotype
				// genotypes come in the order AA, AB, BB, 
				(*m_result)(m_sample_i, m_result_column) += entry_i * value ;
			}
			(*m_nonmissingness)(m_sample_i,m_result_column) = 1 ;
		}

		void set_value( std::size_t entry_i, Integer const value ) {
			assert(0) ; // expecting GT field
		}

		void finalise() {
			// nothing to do
		}
		
	private:
		Eigen::MatrixXd* m_result ;
		Eigen::MatrixXd* m_nonmissingness ;
		genfile::OrderType m_order_type ;
		int m_result_column ;
		std::size_t m_sample_i ;
	} ;
}

void HaplotypeFrequencyComponent::compute_ld_measures(
	genfile::VariantIdentifyingData const& source_snp,
	genfile::VariantDataReader& source_data_reader,
	genfile::VariantIdentifyingData const& target_snp,
	genfile::VariantDataReader& target_data_reader
) {
	if( source_snp.number_of_alleles() == 2 && target_snp.number_of_alleles() == 2 ) {
		std::vector< int > source_calls ;
		std::vector< uint32_t > source_ploidy ;
		std::vector< int > target_calls ;
		std::vector< uint32_t > target_ploidy ;
		Eigen::MatrixXd dosages( source_data_reader.get_number_of_samples(), 2 ) ;
		Eigen::MatrixXd nonmissingness( source_data_reader.get_number_of_samples(), 2 ) ;

		bool haveCalls = false ;
		try {
			source_data_reader.get( ":genotypes:", CallSetter( &source_calls, &source_ploidy ) ) ;
			target_data_reader.get( ":genotypes:", CallSetter( &target_calls, &target_ploidy ) ) ;
			haveCalls = true ;
		} catch( ... ) {
		}

		{
			DosageSetter source_setter( &dosages, &nonmissingness, 0 ) ;
			DosageSetter target_setter( &dosages, &nonmissingness, 1 ) ;
			source_data_reader.get( ":genotypes:", genfile::to_GP_unphased( source_setter ) ) ;
			target_data_reader.get( ":genotypes:", genfile::to_GP_unphased( target_setter ) ) ;
		}
		
		// we treat any data point that is missing in one sample as missing in both
		for( int i = 0; i < nonmissingness.rows(); ++i ) {
			if( nonmissingness.row(i).sum() < 2 ) {
				nonmissingness.row(i).setZero() ;
			}
		}
		
		compute_ld_measures(
			source_snp,
			source_calls,
			source_ploidy,
			target_snp,
			target_calls,
			target_ploidy,
			dosages,
			nonmissingness,
			haveCalls
		) ;
	}
}

void HaplotypeFrequencyComponent::compute_ld_measures(
	genfile::VariantIdentifyingData const& source_snp,
	std::vector< int > const& source_calls,
	std::vector< uint32_t > const& source_ploidy,
	genfile::VariantIdentifyingData const& target_snp,
	std::vector< int > const& target_calls,
	std::vector< uint32_t > const& target_ploidy,
	Eigen::MatrixXd const& dosages,
	Eigen::MatrixXd const& nonmissingness,
	bool const runEM
) {
	if( m_stratification ) {
		bool output = (m_min_r2 == 0.0 ) ;
		for( std::size_t i = 0; i < m_stratification->number_of_strata(); ++i ) {
			try {
				if( runEM ) {
					output = compute_em_ld_measures(
						source_snp,
						source_calls,
						source_ploidy,
						target_snp,
						target_calls,
						target_ploidy,
						m_stratification->stratum_name(i) + ":",
						m_stratification->stratum(i),
						output
					) || output ;
				}
				output = compute_dosage_ld_measures(
					source_snp,
					target_snp,
					dosages,
					nonmissingness,
					m_stratification->stratum_name(i) + ":",
					m_stratification->stratum(i),
					output
				) || output ;
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
			bool output = (m_min_r2 == 0.0 ) ;
			if( runEM ) {
				output = compute_em_ld_measures(
					source_snp,
					source_calls, source_ploidy,
					target_snp,
					target_calls, target_ploidy,
					"",
					std::vector< genfile::SampleRange  >( 1, genfile::SampleRange( 0, dosages.rows() ) ),
					output
				) ;
			}
			compute_dosage_ld_measures(
				source_snp,
				target_snp,
				dosages,
				nonmissingness,
				"",
				std::vector< genfile::SampleRange  >( 1, genfile::SampleRange( 0, dosages.rows() ) ),
				output
			) ;
		} catch( genfile::OperationFailedError const& e ) {
			m_ui_context.logger() << "!! HaplotypeFrequencyComponent::compute_ld_measures(): "
				<< "could not compute LD measures between SNPs "
				<< source_snp << " and " << target_snp << ".\n"
				<< "!! reason: " << e.get_message() << "\n" ;
		}
	}
}

bool HaplotypeFrequencyComponent::compute_dosage_ld_measures(
	genfile::VariantIdentifyingData const& source_snp,
	genfile::VariantIdentifyingData const& target_snp,
	Eigen::MatrixXd const& dosages,
	Eigen::MatrixXd const& nonmissingness,
	std::string const& variable_name_stub,
	std::vector< genfile::SampleRange > const& sample_set,
	bool alwaysOutput
) {
	bool includeInOutput = alwaysOutput ;
	Eigen::RowVectorXd means = Eigen::RowVectorXd::Zero(2) ;
	Eigen::RowVectorXd totals = Eigen::RowVectorXd::Zero(2) ;
	for( std::size_t i = 0; i < sample_set.size(); ++i ) {

#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
		std::cerr << sample_set[i].begin() << "-" << sample_set[i].end() << ":\n" ;
		std::cerr << dosages.rows() << "x" << dosages.cols()
			<< ", " << nonmissingness.rows() << "x" << nonmissingness.cols() << ".\n" ;
		std::cerr << dosages.block( 0, 0, std::min( int( dosages.rows() ), 10 ), dosages.cols() ) << ".\n" ;
#endif
		means.array() += (
			dosages.block( sample_set[i].begin(), 0, sample_set[i].end() - sample_set[i].begin(), 2 ).array()
			* nonmissingness.block( sample_set[i].begin(), 0, sample_set[i].end() - sample_set[i].begin(), 2 ).array()
		).colwise().sum() ;
		
		totals += nonmissingness.block( sample_set[i].begin(), 0, sample_set[i].end() - sample_set[i].begin(), 2 ).colwise().sum() ;
	}

	means.array() /= totals.array() ;

	Eigen::Matrix2d covariance = Eigen::Matrix2d::Zero() ;
	Eigen::Matrix2d nonmissing = Eigen::Matrix2d::Zero() ;
	for( std::size_t i = 0; i < sample_set.size(); ++i ) {
		//Eigen::MatrixBase< Eigen::MatrixXd > block = (
		Eigen::MatrixXd block = (
			(dosages.block( sample_set[i].begin(), 0,  sample_set[i].end() - sample_set[i].begin(), 2 ).rowwise() - means ).array()
				* nonmissingness.block( sample_set[i].begin(), 0,  sample_set[i].end() - sample_set[i].begin(), 2 ).array()
		) ;
		Eigen::Block< Eigen::MatrixXd const > nonmissing_block = nonmissingness.block(
			sample_set[i].begin(), 0,  sample_set[i].end() - sample_set[i].begin(), 2
		) ;

		covariance += block.transpose() * block ;
		nonmissing += nonmissing_block.transpose() * nonmissing_block ;
	}
	// Compute correlation as:
	// 
	double r = covariance(0,1) / std::sqrt( covariance(0,0) * covariance(1,1) ) ;

	includeInOutput = includeInOutput || (r*r >= m_min_r2 ) ;

#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
	std::cerr << "means = " << means << ".\n" ;
	std::cerr << "cov =\n" << covariance << ".\n" ;
	std::cerr << "r =\n" << r << ".\n" ;
#endif

	if( includeInOutput ) {
		m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "dosage_r", r  ) ;
		m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "dosage_r2", r*r ) ;
	}
	return includeInOutput ;
}

bool HaplotypeFrequencyComponent::compute_em_ld_measures(
	genfile::VariantIdentifyingData const& source_snp,
	std::vector< int > const& source_calls,
	std::vector< uint32_t > const& source_ploidy,
	genfile::VariantIdentifyingData const& target_snp,
	std::vector< int > const& target_calls,
	std::vector< uint32_t > const& target_ploidy,
	std::string const& variable_name_stub,
	std::vector< genfile::SampleRange > const& sample_set,
	bool alwaysOutput
) {
	// Construct table of genotypes at each SNP.
	Eigen::Matrix3d diploid_table = Eigen::Matrix3d::Zero() ;
	Eigen::Matrix2d haploid_table = Eigen::Matrix2d::Zero() ;
	
	tabulate_calls(
		sample_set,
		source_calls, source_ploidy,
		target_calls, target_ploidy,
		&diploid_table,
		&haploid_table
	) ;
	
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
	std::cerr << "SNP1: " << source_snp << "\n"
		<< "SNP2: " << target_snp << "\n"
		<< "diploids: " << diploid_table << ".\n"
		<< "haploids: " << haploid_table << ".\n" ;
	std::cerr << "alwaysOutput is " << ( alwaysOutput ? "true" : "false" ) << ".\n" ;
#endif

	genfile::VariantEntry::Integer const N = diploid_table.sum() + haploid_table.sum() ;
	bool includeInOutput = alwaysOutput ;
	genfile::VariantEntry const missing = genfile::MissingValue() ;

	if( includeInOutput ) {
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "number_of_genotypes", diploid_table.sum() ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "number_of_haplotypes", haploid_table.sum() ) ;
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
			HaplotypeFrequencyLogLikelihood ll( diploid_table, haploid_table + m_prior ) ;
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
#if DEBUG_HAPLOTYPE_FREQUENCY_COMPONENT
			std::cerr << "r = \n" << std::setprecision(3) << r << ".\n" ;
			std::cerr << "r^2 = " << r2 << ".\n" ;
			std::cerr << "D =\n" << D << ".\n" ;
			std::cerr << "Dprime =\n" << Dprime << ".\n" ;
			std::cerr << "pi = " << pi.transpose() << "\n" ;
#endif
			
			includeInOutput = includeInOutput || (r2 >= m_min_r2 ) ;
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
				<< "!! table is:\n" << diploid_table << ".\n" ;

			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi00", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi01", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi10", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "pi11", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "D", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "Dprime", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "r", missing ) ;
			m_sink->store_per_variant_pair_data( source_snp, target_snp, variable_name_stub + "r2", missing ) ;
		}
	

	}
	return includeInOutput ;
}

void HaplotypeFrequencyComponent::end_processing_snps() {
	if( m_sink.get() ) {
		m_sink->finalise() ;
	}
}

void HaplotypeFrequencyComponent::send_results_to( haplotype_frequency_component::FlatTableDBOutputter::UniquePtr sink ) {
	m_sink = sink ;
}
