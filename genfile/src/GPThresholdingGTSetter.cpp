
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <cassert>
#include "genfile/VariantDataReader.hpp"
#include "genfile/Error.hpp"
#include "genfile/GPThresholdingGTSetter.hpp"


namespace genfile {
	GPThresholdingGTSetter::GPThresholdingGTSetter( VariantDataReader::PerSampleSetter& target, double const threshhold ):
		m_target( target ),
		m_threshhold( threshhold ),
		m_number_of_alleles(0),
		m_order_type( eUnknownOrderType ),
		m_ploidy(0)
	{
		if( m_threshhold <= 0.5 ) {
			throw BadArgumentError(
				"genfile::ThreshholdingDataReader::GPThresholdingGTSetter::GPThresholdingGTSetter()",
				"threshhold",
				"Only threshholds greater than 0.5 are supported."
			) ;
		}
		assert( m_threshhold > 0.5 ) ;
	}
	
	GPThresholdingGTSetter::~GPThresholdingGTSetter() throw() {}
		
	void GPThresholdingGTSetter::initialise( std::size_t nSamples, std::size_t nAlleles ) {
		m_number_of_alleles = nAlleles ;
		m_target.initialise( nSamples, nAlleles ) ;
	} ;

	bool GPThresholdingGTSetter::set_sample( std::size_t i ) {
		m_target.set_sample( i ) ;
		m_missing = false ;
		m_order_type = eUnknownOrderType ;
		m_ploidy = 0 ;
		m_entry_i = 0 ;
		return true ;
	}
	
	void GPThresholdingGTSetter::set_number_of_entries( uint32_t ploidy, std::size_t number_of_entries, OrderType const order_type, ValueType const value_type ) {
		assert( order_type == ePerUnorderedGenotype || ePerPhasedHaplotypePerAllele ) ;
		assert( value_type == eProbability ) ;
		m_order_type = order_type ;
		
		m_ploidy = ploidy ;
		m_calls.setZero( number_of_entries ) ;
		m_target.set_number_of_entries(
			m_ploidy, m_ploidy,
			( m_order_type == ePerUnorderedGenotype ) ? ePerUnorderedHaplotype : ePerOrderedHaplotype,
			eAlleleIndex
		) ;
	}
	
	void GPThresholdingGTSetter::set_value( std::size_t, MissingValue const value ) {
		m_missing = true ;
		m_entry_i++ ;
		if( m_entry_i == m_calls.size() ) { 
			send_results() ;
		}
	}

	void GPThresholdingGTSetter::set_value( std::size_t, std::string& value ) {
		throw BadArgumentError(
			"genfile::ThreshholdingDataReader::GPThresholdingGTSetter::set_value",
			"value",
			"Expected a floating-point value (got a string)."
		) ;
	}

	void GPThresholdingGTSetter::set_value( std::size_t, Integer const value ) {
		throw BadArgumentError(
			"genfile::ThreshholdingDataReader::GPThresholdingGTSetter::set_value",
			"value",
			"Expected a floating-point value (got an integer)."
		) ;
	}

	void GPThresholdingGTSetter::set_value( std::size_t, double const value ) {
		m_calls( m_entry_i++ ) = value ;
		if( m_entry_i == m_calls.size() ) { 
			send_results() ;
		}
	}
	
	void GPThresholdingGTSetter::finalise() {
		m_target.finalise() ;
	}

	void GPThresholdingGTSetter::send_results() {
		if( m_missing || m_calls.array().maxCoeff() < m_threshhold ) {
			for( uint32_t i = 0 ; i < m_ploidy; ++i ) {
				m_target.set_value( i, genfile::MissingValue() ) ;
			}
		} else {
			// How to figure out the order here?
			// Well..it goes like this.
			if( m_order_type == ePerPhasedHaplotypePerAllele ) {
				// for phased data it is simple:
				for( uint32_t i = 0; i < m_ploidy; ++i ) {
					for( uint32_t j = 0; j < m_number_of_alleles; ++j ) {
						if( m_calls[i*m_number_of_alleles+j] >= m_threshhold ) {
							m_target.set_value( i, Integer( j ) ) ;
							break ;
						}
					}
				}
			} else if( m_order_type == ePerUnorderedGenotype ) {
				// For unphased data it is more complex.
				// We have m_ploidy = n chromosomes in total.
				// Genotypes are all ways to put n_alleles = k alleles into
				// those chromosomes.
				// Represent these as k-vectors that sum to n (i.e. v_i is the count of allele i)
				// Or (k-1)-vectors that sum to at most n.
				// The order is chosen so that lower indices are always used first.
				// 3,0,0 = AAA
				// 2,1,0 = AAB
				// 1,2,0 = ABB
				// 0,3,0 = BBB
				// 2,0,1 = BBC
				// 1,1,1 = ABC
				// 0,2,1 = BBC
				// 0,1,2 = BCC
				// 0,0,3 = CCC
			

				// We fast-path the common (diploid, biallelic) case.
				if( m_number_of_alleles == 2 && m_ploidy == 2 ) {
					if( m_calls[0] >= m_threshhold ) {
						m_target.set_value( 0, Integer(0) ) ;
						m_target.set_value( 1, Integer(0) ) ;
					} else if( m_calls[1] >= m_threshhold ) {
						m_target.set_value( 0, Integer(0) ) ;
						m_target.set_value( 1, Integer(1) ) ;
					} else if( m_calls[2] >= m_threshhold ) {
						m_target.set_value( 0, Integer(1) ) ;
						m_target.set_value( 1, Integer(1) ) ;
					} else {
						m_target.set_value( 0, genfile::MissingValue() ) ;
						m_target.set_value( 1, genfile::MissingValue() ) ;
					}
				} else {
					// We enumerate vectors in this order and stop when we meet a probabiltiy
					// meeting the desired threshhold.
					std::vector< uint16_t > limits( (m_number_of_alleles-1), m_ploidy ) ;
					std::vector< uint16_t > allele_counts( m_number_of_alleles, 0 ) ;
					allele_counts[0] = m_ploidy ;
					bool finished = false ;
					for( std::size_t index = 0; !finished; ++index ) {
						if( m_calls[index] >= m_threshhold ) {
							std::size_t this_index = 0 ;
							for( std::size_t allele = 0; allele < m_number_of_alleles; ++allele ) {
								for( uint16_t count = 0; count < allele_counts[allele]; ++count ) {
									m_target.set_value( ++this_index, Integer( allele )) ;
								}
							}
							break ;
						}
						
						// Move to next possible genotype
						std::size_t j = 0 ;
						for( ; j < (m_number_of_alleles-1); ++j ) {
							uint16_t value = allele_counts[j+1] ;
							if( value < limits[ j ] ) {
								++allele_counts[j+1] ;
								--allele_counts[0] ;
								for( std::size_t k = 0; k < j; ++k ) {
									--limits[k] ;
								}
								break ;
							} else {
								// allele count has reached its limit.
								// Reset it to zero.
								// Note that to get here all lower-order counts must be zero.
								allele_counts[j+1] = 0 ;
								allele_counts[0] += value ;
								for( std::size_t k = 0; k < j; ++k ) {
									limits[k] += value ;
								}
							}
						}
						if( j == (m_number_of_alleles-1) ) {
							finished = true ;
						}
					}
				}
			}
		}
	}
} ;

