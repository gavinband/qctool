
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/thread.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/VariantEntry.hpp"
#include "genfile/vcf/get_set_eigen.hpp"
#include "genfile/Error.hpp"
#include "genfile/ToGP.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputation.hpp"
#include "components/SNPSummaryComponent/SNPSummaryComputationManager.hpp"
#include "components/SNPSummaryComponent/StratifyingSNPSummaryComputation.hpp"

// #define DEBUG_SNP_SUMMARY_COMPUTATION_MANAGER 1

namespace stats {
	SNPSummaryComputationManager::SNPSummaryComputationManager(
		genfile::CohortIndividualSource const& samples,
		std::string const& sex_column_name
	):
		m_samples( samples ),
		m_haploid_coding_column( -1 )
	{}

	std::vector< char > SNPSummaryComputationManager::get_sexes( genfile::CohortIndividualSource const& samples, std::string const& sex_column_name ) const {
		std::vector< char > result( samples.get_number_of_individuals(), '.' ) ;
		genfile::CohortIndividualSource::ColumnSpec column_spec = samples.get_column_spec() ;
		if( column_spec.check_for_column( sex_column_name ) ) {
			if( column_spec[ sex_column_name ].is_discrete() ) {
				for( std::size_t i = 0; i < samples.get_number_of_individuals(); ++i ) {
					genfile::CohortIndividualSource::Entry const entry = samples.get_entry( i, sex_column_name ) ;
					if( !entry.is_missing() ) {
						std::string const sex = genfile::string_utils::to_lower( entry.as< std::string >() ) ;
						if( sex == "m" || sex == "male" || sex == "1" ) {
							result[i] = 'm' ;
						} else if( sex == "f" || sex == "female" || sex == "2" ) {
							result[i] = 'f' ;
						}
						else {
							throw genfile::MalformedInputError( samples.get_source_spec(), i+2, column_spec.find_column( sex_column_name )) ;
						}
					}
				}
			}
			else {
				std::cerr << "!! (SNPSummaryComputationManager::get_sexes): sex column found but it has the wrong type!\n" ;
				throw genfile::MalformedInputError( samples.get_source_spec(), 1, column_spec.find_column( sex_column_name )) ;
			}
		}
		return result ;
	}

	std::map< char, std::vector< int > > SNPSummaryComputationManager::get_samples_by_sex( std::vector< char > const& sexes ) const {
		std::map< char, std::vector< int > > result ;
		result[ 'm' ] ;
		result[ 'f' ] ;
		result[ '.' ] ;

		for( std::size_t i = 0; i < sexes.size(); ++i ) {
			result[ sexes[i] ].push_back( i ) ;
		}

		return result ;
	}

	void SNPSummaryComputationManager::add_computation( std::string const& name, SNPSummaryComputation::UniquePtr computation ) {
		m_computations.insert( name, computation ) ;
	}

	void SNPSummaryComputationManager::add_result_callback( ResultCallback callback ) {
		m_result_signal.connect( callback ) ;
	}

	void SNPSummaryComputationManager::add_per_sample_result_callback( PerSampleResultCallback callback ) {
		m_per_sample_result_signal.connect( callback ) ;
	}

	void SNPSummaryComputationManager::set_haploid_genotype_coding( int coding ) {
		assert( coding == 1 || coding == 2 ) ;
		m_haploid_coding_column = coding ;	
	}

	void SNPSummaryComputationManager::begin_processing_snps( std::size_t number_of_samples, genfile::SNPDataSource::Metadata const& ) {
		m_snp_index = 0 ;
		m_genotypes.resize( number_of_samples, 3 ) ;
		Computations::iterator i = m_computations.begin(), end_i = m_computations.end() ;
		for( ; i != end_i; ++i ) {
			i->second->begin_processing_snps( number_of_samples ) ;
		}
	}

	namespace {
		struct GPSetter: public genfile::VariantDataReader::PerSampleSetter {
			~GPSetter() throw() {}
			GPSetter(
				Eigen::MatrixXd* genotypes,
				SNPSummaryComputation::Ploidy* ploidy
			):
				m_genotypes( genotypes ),
				m_ploidy( ploidy ),
				m_sample_i( 0 )
			{
				assert( genotypes != 0 ) ;
				assert( ploidy != 0 ) ;
			}

			void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
	#if DEBUG_PER_VARIANT_COMPUTATION_MANAGER
				std::cerr << "GPSetter::initialise( " << number_of_samples << ", " << number_of_alleles << " ).\n" ;
	#endif
				// assume at most diploidy, max entries is number of alleles choose 2.
				m_genotypes->resize( number_of_samples, number_of_alleles * (number_of_alleles+1)/2 ) ;
				m_genotypes->setZero() ;
				m_ploidy->resize( number_of_samples ) ;
				m_ploidy->setConstant( -1 ) ;
				m_sample_i = 0 ;
				m_max_ploidy = 0 ;
				m_min_ploidy = std::numeric_limits< uint32_t >::max() ;
			}

			bool set_sample( std::size_t i ) {
				m_sample_i = i ;
				return true ;
			}
		
			void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
				assert( order_type == genfile::ePerUnorderedGenotype ) ;
				assert( value_type == genfile::eProbability ) ;
				(*m_ploidy)( m_sample_i ) = ploidy ;
				m_min_ploidy = std::min( ploidy, m_min_ploidy ) ;
				m_max_ploidy = std::max( ploidy, m_max_ploidy ) ;
			}
		
			void set_value( std::size_t value_i, genfile::MissingValue const value ) {
				// GEnotypes already set to 0, nothing to do.
			}

			void set_value( std::size_t value_i, double const value ) {
	#if DEBUG_PER_VARIANT_COMPUTATION_MANAGER
				std::cerr << "GPSetter::set_value( " << value_i << ", " << value << " ).\n" ;
				std::cerr << "GPSetter: m_ploidy[" << m_sample_i << "] = " << m_ploidy[ m_sample_i ] << ".\n" ;
	#endif
				assert( value_i < m_genotypes->cols() ) ;
				(*m_genotypes)( m_sample_i, value_i ) = value ;
			}

			void finalise() {}

			uint32_t max_ploidy() const { return m_max_ploidy ; }
			uint32_t min_ploidy() const { return m_min_ploidy ; }

		private:
			Eigen::MatrixXd* m_genotypes ;
			SNPSummaryComputation::Ploidy* m_ploidy ;
			int m_sample_i ;
			uint32_t m_max_ploidy ;
			uint32_t m_min_ploidy ;
		} ;
	}

	void SNPSummaryComputationManager::processed_snp(
		genfile::VariantIdentifyingData const& snp,
		genfile::VariantDataReader& data_reader
	) {
		try {
			GPSetter setter( &m_genotypes, &m_ploidy ) ;
			data_reader.get( ":genotypes:", genfile::to_GP_unphased( setter )) ;

	#if DEBUG_SNP_SUMMARY_COMPUTATION_MANAGER
			std::cerr << "SNPSummaryComputationManager::processed_snp(): ploidy = " << m_ploidy.transpose() << "...\n" ;
	#endif

			boost::function< void ( std::string const& value_name, genfile::VariantEntry const& value ) > callback
				= boost::bind( boost::ref( m_result_signal ), snp, _1, _2 ) ;

			Computations::iterator i = m_computations.begin(), end_i = m_computations.end() ;
			for( ; i != end_i; ++i ) {
				i->second->operator()(
					snp,
					m_genotypes,
					m_ploidy,
					data_reader,
					callback
				) ;
			}

			if( snp.number_of_alleles() != 2 ) {
				m_result_signal( snp, "comment", "non-biallelic" ) ;
			}
		}
		catch( genfile::MalformedInputError const& e ) {
			m_result_signal( snp, "comment", "!! Error reading data for variant " + genfile::string_utils::to_string( snp ) + ": " + e.format_message() ) ;
		}
		++m_snp_index ;
	}

	void SNPSummaryComputationManager::end_processing_snps() {
		Computations::iterator i = m_computations.begin(), end_i = m_computations.end() ;
		for( ; i != end_i; ++i ) {
			i->second->end_processing_snps(
				boost::bind( boost::ref( m_per_sample_result_signal ), _1, _2, _3 )
			) ;
		}
	}

	void SNPSummaryComputationManager::stratify_by( StrataMembers const& strata, std::string const& stratification_name ) {
		Computations::iterator computation_i = m_computations.begin(), end_computations = m_computations.end() ;
		for( ; computation_i != end_computations; ) {
			Computations::iterator this_computation_i = computation_i++ ;
			std::string const this_name = this_computation_i->first ;
			SNPSummaryComputation::UniquePtr computation( m_computations.release( this_computation_i ).release() ) ;
			m_computations.insert(
				this_name,
				SNPSummaryComputation::UniquePtr( new StratifyingSNPSummaryComputation( computation, stratification_name, strata ))
			) ;
		}
	}

	std::string SNPSummaryComputationManager::get_summary( std::string const& prefix, std::size_t column_width ) {
		std::string result ;
		Computations::const_iterator i = m_computations.begin(), end_i = m_computations.end() ;
		for( ; i != end_i; ++i ) {
			result += i->second->get_summary( prefix, column_width ) + "\n";
		}
		return result ;
	}
}
