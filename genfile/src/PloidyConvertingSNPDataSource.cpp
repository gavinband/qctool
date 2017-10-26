
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <vector>
#include <numeric>
#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/PloidyConvertingSNPDataSource.hpp"
#include "genfile/Error.hpp"

//#define DEBUG 1
#if DEBUG
#include <iostream>
#include <iomanip>
#endif

namespace genfile {
	PloidyConvertingSNPDataSource::UniquePtr PloidyConvertingSNPDataSource::create(
		SNPDataSource::UniquePtr source,
		GetPloidy get_ploidy
	) {
		return UniquePtr(
			new PloidyConvertingSNPDataSource( source, get_ploidy )
		) ;
	}

	PloidyConvertingSNPDataSource::PloidyConvertingSNPDataSource(
		SNPDataSource::UniquePtr source,
		GetPloidy get_ploidy
	): 
		m_source( source ),
		m_get_ploidy( get_ploidy )
	{}

	namespace {
		struct PloidyConvertingSetter: public VariantDataReader::PerSampleSetter
		{
			PloidyConvertingSetter(
				VariantDataReader::PerSampleSetter& setter,
				std::vector< int > const& ploidy
			):
				m_setter( setter ),
				m_ploidy( ploidy ),
				m_number_of_alleles( 0 ),
				m_convert( false ),
				m_original_ploidy(0),
				m_target_ploidy(0),
				m_original_order_type(eUnknownOrderType)
			{
				using genfile::string_utils::to_string ;
				for( std::size_t i = 0; i < m_ploidy.size(); ++i ) {
					if( m_ploidy[i] > 2 ) {
						throw genfile::BadArgumentError(
							"genfile::PloidyConvertingSetter::PloidyConvertingSetter()",
							"m_ploidy[" + to_string(i) + "]=" + to_string( m_ploidy[i] ),
							"I only handle conversion to zero-, haploid, or diploid state"
						) ;
					}
				}
			}
			
			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				m_setter.initialise( nSamples, nAlleles ) ;
				m_number_of_alleles = nAlleles ;
				m_convert = false ;
			}

			bool set_sample( std::size_t n ) {
				m_target_ploidy = m_ploidy[n] ;
				return m_setter.set_sample( n ) ;
			}
			
			void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType order_type, ValueType value_type ) {
				m_original_ploidy = ploidy ;
				m_original_order_type = order_type ;
				m_original_number_of_entries = n ;
				if( m_original_ploidy != 2 ) {
					throw genfile::BadArgumentError(
						"genfile::PloidyConvertingSetter::set_number_of_entries()",
						"ploidy=" + genfile::string_utils::to_string( ploidy ),
						"Expected diploidy"
					) ;
				}
				
				// Only convert diploid calls to haploid or 0-ploid calls
				// And only genotype calls or probabilities
				m_convert = (m_target_ploidy < 2) && (
					// unphased genotype probs (GP)
					(order_type == ePerUnorderedGenotype && value_type == eProbability) 
					// phased haplotype probs
					|| (order_type == ePerPhasedHaplotypePerAllele && value_type == eProbability)
					// phased genotype calls
					|| (order_type == ePerOrderedHaplotype && value_type == eAlleleIndex) 
					// unphased genotype calls
					|| (order_type == ePerUnorderedHaplotype && value_type == eAlleleIndex)
				) ;
				if( m_convert ) {
					// Store values locally, then deal with them in finalise()
					m_values.resize( n ) ;
				} else {
					// No conversion
					m_setter.set_number_of_entries( ploidy, n, order_type, value_type ) ;
				}
			}

			void set_value( std::size_t entry_i, MissingValue const value ) {
				if( m_convert ) {
					m_values[ entry_i ] = -1 ;
					if( (entry_i+1) == m_original_number_of_entries) {
						finalise_entries() ;
					}
				} else {
					m_setter.set_value( entry_i, value ) ;
				}
			}
			void set_value( std::size_t entry_i, std::string& value ) {
				if( m_convert ) {
					assert(0) ;
				} else {
					m_setter.set_value( entry_i, value ) ;
				}
			}
			void set_value( std::size_t entry_i, Integer const value ) {
				if( m_convert ) {
					m_values[ entry_i ] = value ;
					if( (entry_i+1) == m_original_number_of_entries) {
						finalise_entries() ;
					}
				} else {
					m_setter.set_value( entry_i, value ) ;
				}
			}
			void set_value( std::size_t entry_i, double const value ) {
				if( m_convert ) {
					m_values[ entry_i ] = value ;
					if( (entry_i+1) == m_original_number_of_entries) {
						finalise_entries() ;
					}
				} else {
					m_setter.set_value( entry_i, value ) ;
				}
			}

			void finalise_entries() {
#if DEBUG
				std::cerr << "Original order type: " << m_original_order_type << "\n"
					<< "Original ploidy: " << m_original_ploidy << "\n"
					<< "target ploidy: " << m_target_ploidy << "\n"
					<< (m_convert ? "...converting\n" : "...not converting.\n") ;
#endif
				if( m_convert ) {
					assert( m_original_ploidy == 2 && m_target_ploidy < 2 ) ;
					if( m_target_ploidy == -1 ) {
						// Just set values to missing, whatever they were.
						m_setter.set_number_of_entries( 2, m_values.size(), ePerUnorderedGenotype, eProbability ) ;
						for( std::size_t i = 0; i < m_values.size(); ++i ) {
							m_setter.set_value( i, genfile::MissingValue() ) ;
						}
					} else if( m_original_order_type == ePerUnorderedGenotype ) {
						if( m_target_ploidy == 0 ) {
							finalise_zeroploid_genotype_probs() ;
						} else {
							finalise_haploid_genotype_probs() ;
						}
					} else if( m_original_order_type == ePerPhasedHaplotypePerAllele ) {
						if( m_target_ploidy == 0 ) {
							finalise_zeroploid_genotype_probs() ;
						} else {
							finalise_haploid_haplotype_probs() ;
						}
					} else if(
						m_original_order_type == ePerOrderedHaplotype
						|| m_original_order_type == ePerUnorderedHaplotype
					) {
						if( m_target_ploidy == 0 ) {
							finalise_zeroploid_genotype_calls() ;
						} else {
							finalise_haploid_genotype_calls() ;
						}
					} else {
						// not supported.
						assert(0) ;
					}
				}
			}
			
			void finalise() {
				m_setter.finalise() ;
			}
			
		private:
			VariantDataReader::PerSampleSetter& m_setter ;
			std::vector< int > const& m_ploidy ;
			std::size_t m_number_of_alleles ;
			bool m_convert ;
			uint32_t m_original_ploidy ;
			int m_target_ploidy ;
			OrderType m_original_order_type ;
			std::size_t m_original_number_of_entries ;
			std::vector< double > m_values ;

		private:
			void finalise_zeroploid_genotype_probs() {
				m_setter.set_number_of_entries( 1, 1, ePerUnorderedGenotype, eProbability ) ;
				m_setter.set_value( 0, 1.0 ) ;
			}

			void finalise_zeroploid_genotype_calls() {
				m_setter.set_number_of_entries( 0, 0, ePerOrderedHaplotype, eAlleleIndex ) ;
			}

			void finalise_haploid_genotype_probs() {
				assert( m_original_ploidy == 2 ) ;
				
#if DEBUG
				std::cerr << "m_values[1] = " << std::setprecision( 100 ) << m_values[1] << "\n" ;
				std::cerr << "0.0 + m_values[1] = " << std::setprecision( 100 ) << (0.0 + m_values[1]) << "\n" ;
#endif
				double totalHet = 0.0 ;
				std::size_t indexOfHom = 0 ;
				for( std::size_t i = 0; i < m_number_of_alleles; ++i ) {
					// For ploidy = 2,  hom calls come at indices (triangular numbers-1)
					// e.g. (1, 3, 6, ...) - 1.
					// e.g.
					// 0   1   2   3   4   5   6   7   8   9
					// AA, AB, BB, AC, BC, CC, AD, BD, CD, DD
					//
					// This loop moves the hom call probabilities
					// to the start of the array, and at the same time
					// accumulates het calls at invervening indices to test
					// if any are nonzero.
					if( (i+1) < m_number_of_alleles ) {
						totalHet += std::accumulate(
							m_values.begin() + indexOfHom + 1,
							m_values.begin() + indexOfHom + i+2,
							0.0
						) ;
#if DEBUG
						std::cerr << "Added values [" << indexOfHom + 1 << ", " << indexOfHom + i+2 << ")\n" ;
#endif
					}
					
					// move hom call probability to front of values
					m_values[i] = m_values[indexOfHom] ;
					indexOfHom += (i+2) ;
				}
#if DEBUG
				std::cerr << "m_values[1] = " << std::setprecision( 100 ) << m_values[1] << "\n" ;
				std::cerr << "totalHet = " << std::setprecision( 100 ) << totalHet << "\n" ;
				std::cerr << "....it is " << (( totalHet > 0.0 ) ? "greater than" : "not greater than") << " 0.\n" ;
#endif
				m_setter.set_number_of_entries( 1, m_number_of_alleles, ePerUnorderedGenotype, eProbability ) ;
				if( totalHet > 0.0 ) {
					for( std::size_t i = 0; i < m_number_of_alleles; ++i ) {
						m_setter.set_value( i, genfile::MissingValue() ) ;
					}
				} else {
					for( std::size_t i = 0; i < m_number_of_alleles; ++i ) {
						if( m_values[i] == -1 ) {
							m_setter.set_value( i, genfile::MissingValue() ) ;
						} else {
							m_setter.set_value( i, m_values[i] ) ;
						}
					}
				}
			}

			void finalise_haploid_haplotype_probs() {
				assert( m_original_ploidy == 2 ) ;
				// We can't handle this currently so set all values to missing
				// (Note that because this reports per-haplotype probabilities independently,
				// the only case we could handle would be where there is no uncertainty i.e.
				// all haps certainly equal to the same value.
				m_setter.set_number_of_entries( 1, m_number_of_alleles, ePerUnorderedGenotype, eProbability ) ;
				for( std::size_t i = 0; i < m_number_of_alleles; ++i ) {
					m_setter.set_value( i, genfile::MissingValue() ) ;
				}
			}
			
			void finalise_haploid_genotype_calls() {
				m_setter.set_number_of_entries( 1, 1, ePerOrderedHaplotype, eAlleleIndex ) ;
				// Check all values are equal, i.e. homozygous call
				double const value = m_values.front() ;
				bool missing = (value == -1) ;
				for( std::size_t i = 1; i < m_values.size() && !missing; ++i ) {
					if( m_values[i] != value ) {
						missing = true ;
						break ;
					}
				}
				if( missing || value == -1 ) {
					m_setter.set_value( 0, genfile::MissingValue() ) ;
				} else {
					m_setter.set_value( 0, Integer( value ) ) ;
				}
			}
		} ;

		struct PloidyConvertingVariantDataReader: public VariantDataReader
		{
			static UniquePtr create(
				VariantDataReader::UniquePtr reader,
				std::vector< int > const& ploidy
			) {
				return UniquePtr(
					new PloidyConvertingVariantDataReader( reader, ploidy )
				) ;
			}

			PloidyConvertingVariantDataReader(
				VariantDataReader::UniquePtr reader,
				std::vector< int > const& ploidy
			):
				m_reader( reader ),
				m_ploidy( ploidy )
			{}

			PloidyConvertingVariantDataReader& get( std::string const& spec, VariantDataReader::PerSampleSetter& setter ) {
				PloidyConvertingSetter convertingSetter( setter, m_ploidy ) ;
				m_reader->get( spec, convertingSetter ) ;
				return *this ;
			}
			
			std::size_t get_number_of_samples() const { return m_reader->get_number_of_samples() ; }
			
			bool supports( std::string const& spec ) const {
				return m_reader->supports( spec ) ;
			}

			void get_supported_specs( SpecSetter setter ) const {
				m_reader->get_supported_specs( setter ) ;
			}
			
		private:
			VariantDataReader::UniquePtr m_reader ;
			std::vector< int > m_ploidy ;
		} ;
	}

	VariantDataReader::UniquePtr PloidyConvertingSNPDataSource::read_variant_data_impl() {
		if( m_variant.get_position().chromosome().is_sex_determining() ) {
			// compute haploid samples
			m_ploidy.resize( m_source->number_of_samples() ) ;
			for( std::size_t i = 0; i < m_ploidy.size(); ++i ) {
				m_ploidy[i] = m_get_ploidy( m_variant, i ) ;
			}
			return PloidyConvertingVariantDataReader::create(
				m_source->read_variant_data(),
				m_ploidy
			) ;
		} else {
			return m_source->read_variant_data() ;
		}
	}
}
