
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <queue>
#include <memory>
#include <boost/format.hpp>
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/VariantDataReader.hpp"
#include "genfile/ThreshholdingSNPDataSource.hpp"

namespace genfile {
// This SNPDataSource applies a threshhold to genotype calls.
	ThreshholdingSNPDataSource::ThreshholdingSNPDataSource( SNPDataSource::UniquePtr source, double const threshhold ):
		m_source( source ),
		m_threshhold( threshhold )
	{}
	
	ThreshholdingSNPDataSource::operator bool() const {
		return *m_source ;
	}

	SNPDataSource::Metadata ThreshholdingSNPDataSource::get_metadata() const {
		// Remove GT from metadata
		Metadata metadata = m_source->get_metadata() ;
		for( Metadata::iterator i = metadata.begin(); i != metadata.end(); ) {
			Metadata::iterator erase_i = i++ ;
			if( erase_i->first == "FORMAT" ) {
				if( erase_i->second.find( "ID" ) != erase_i->second.end() && erase_i->second[ "ID" ] == "GT" ) {
					metadata.erase( erase_i ) ;
				}
			}
		}
		std::map< std::string, std::string > format ;
		format[ "ID" ] = "GT" ;
		format[ "Number" ] = "1" ;
		format[ "Type" ] = "String" ;
		format[ "Description" ] = ( boost::format( "Genotype call probabilities, threshholded at %.2f" ) % m_threshhold ).str() ;
		metadata.insert( std::make_pair( "FORMAT", format )) ;
		return metadata ;
	}

	unsigned int ThreshholdingSNPDataSource::number_of_samples() const {
		return m_source->number_of_samples() ;
	}
	void ThreshholdingSNPDataSource::get_sample_ids( GetSampleIds getter ) const {
		return m_source->get_sample_ids( getter ) ;
	}
	SNPDataSource::OptionalSnpCount ThreshholdingSNPDataSource::total_number_of_snps() const {
		return m_source->total_number_of_snps() ;
	}
	std::string ThreshholdingSNPDataSource::get_source_spec() const {
		return "ThreshholdingSNPDataSource(" + m_source->get_source_spec() + ")" ;
	}
	SNPDataSource const& ThreshholdingSNPDataSource::get_parent_source() const {
		return *m_source ;
	}
	SNPDataSource const& ThreshholdingSNPDataSource::get_base_source() const {
		return m_source->get_base_source() ;
	}
	std::string ThreshholdingSNPDataSource::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return "ThreshholdingSNPDataSource(" + m_source->get_source_spec() + ")" ;
	}
	void ThreshholdingSNPDataSource::get_snp_identifying_data_impl( VariantIdentifyingData* variant	) {
		m_source->get_snp_identifying_data( variant ) ;
	}

	namespace {
		struct ThreshholdingDataReader: public VariantDataReader {
			ThreshholdingDataReader( VariantDataReader::UniquePtr source, double const threshhold ):
				m_source( source ),
				m_threshhold( threshhold )
			{}

			~ThreshholdingDataReader() {}

			struct Threshholder: public VariantDataReader::PerSampleSetter {
				Threshholder( PerSampleSetter& target, double const threshhold ):
					m_target( target ),
					m_threshhold( threshhold ),
					m_number_of_alleles(0),
					m_order_type( eUnknownOrderType ),
					m_ploidy(0)
				{
					if( m_threshhold <= 0.5 ) {
						throw BadArgumentError(
							"genfile::ThreshholdingDataReader::Threshholder::Threshholder()",
							"threshhold",
							"Only threshholds greater than 0.5 are supported."
						) ;
					}
					assert( m_threshhold > 0.5 ) ;
				}
				
				~Threshholder() throw() {}
					
				void initialise( std::size_t nSamples, std::size_t nAlleles ) {
					m_number_of_alleles = nAlleles ;
					m_target.initialise( nSamples, nAlleles ) ;
				} ;
				bool set_sample( std::size_t i ) {
					m_target.set_sample( i ) ;
					m_missing = false ;
					m_order_type = eUnknownOrderType ;
					m_ploidy = 0 ;
					m_entry_i = 0 ;
					return true ;
				}
				void set_number_of_entries( uint32_t ploidy, std::size_t number_of_entries, OrderType const order_type, ValueType const value_type ) {
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
				void set_value( std::size_t, MissingValue const value ) {
					m_missing = true ;
					m_entry_i++ ;
					if( m_entry_i == m_calls.size() ) { 
						send_results() ;
					}
				}
				void set_value( std::size_t, std::string& value ) {
					throw BadArgumentError(
						"genfile::ThreshholdingDataReader::Threshholder::set_value",
						"value",
						"Expected a floating-point value (got a string)."
					) ;
				}
				void set_value( std::size_t, Integer const value ) {
					throw BadArgumentError(
						"genfile::ThreshholdingDataReader::Threshholder::set_value",
						"value",
						"Expected a floating-point value (got an integer)."
					) ;
				}
				void set_value( std::size_t, double const value ) {
					m_calls( m_entry_i++ ) = value ;
					if( m_entry_i == m_calls.size() ) { 
						send_results() ;
					}
				}
				
				void finalise() {}
			private:
				PerSampleSetter& m_target ;
				double const m_threshhold ;
				std::size_t m_number_of_alleles ;
				OrderType m_order_type ;
				uint32_t m_ploidy ;
				int m_entry_i ;
				Eigen::VectorXd m_calls ;
				bool m_missing ;
				
			private:
				void send_results() {
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
									if( m_calls[i*m_number_of_alleles+j] > m_threshhold ) {
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
							
							// We enumerate vectors in this order and stop when we meet a probabiltiy
							// meeting the desired threshhold.
							
							// This code handles up to 
							std::vector< uint16_t > limits( (m_number_of_alleles-1), m_ploidy ) ;
							std::vector< uint16_t > allele_counts( m_number_of_alleles, 0 ) ;
							allele_counts[0] = m_ploidy ;
							bool finished = false ;
							for( std::size_t index = 0; !finished; ++index ) {
								if( m_calls[index] > m_threshhold ) {
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
			} ;

			VariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				if( spec == ":genotypes:" || spec == "GT" ) {
					if( !m_source->supports( "GP" )) {
						throw genfile::BadArgumentError(
							"genfile::ThreshholdingDataReader::Threshholder::get()",
							"setter",
							"Underlying source must support genotype probabilities (GP) field."
						) ;
					}
					Threshholder threshholder( setter, m_threshhold ) ;
					m_source->get( "GP", threshholder ) ;
				} else {
					m_source->get( spec, setter ) ;
				}
				return *this ;
			}

			bool supports( std::string const& spec ) const {
				return spec == ":genotypes:" || spec == "GT" || m_source->supports( spec ) ;
			}
			void get_supported_specs( SpecSetter setter ) const {
				setter( ":genotypes:", "Integer" ) ;
				setter( "GT", "Integer" ) ;
				return m_source->get_supported_specs( setter ) ;
			}
			std::size_t get_number_of_samples() const {
				return m_source->get_number_of_samples() ;
			}
		private:
			VariantDataReader::UniquePtr m_source ;
			std::string const m_field ;
			double const m_threshhold ;
		} ;
	}

	VariantDataReader::UniquePtr ThreshholdingSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr(
			new ThreshholdingDataReader(
				m_source->read_variant_data(),
				m_threshhold
			)
		) ;
	}

	void ThreshholdingSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}

	void ThreshholdingSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}
}

