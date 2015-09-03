
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <queue>
#include <memory>
#include <boost/format.hpp>
#include "genfile/SNPIdentifyingData.hpp"
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
	void ThreshholdingSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		m_source->get_snp_identifying_data(
			set_number_of_samples, set_SNPID, set_RSID, set_chromosome, set_SNP_position, set_allele1, set_allele2
		) ;
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
					m_threshhold( threshhold )
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
					m_target.initialise( nSamples, nAlleles ) ;
					m_calls.setZero( ( nAlleles * (nAlleles+1))/2 ) ;
				} ;
				bool set_sample( std::size_t i ) {
					m_target.set_sample( i ) ;
					m_missing = false ;
					m_calls.setZero( m_calls.size() ) ;
					m_entry_i = 0 ;
					return true ;
				}
				void set_number_of_entries( std::size_t n, OrderType const order_type, ValueType const value_type ) {
					assert( n == m_calls.size() ) ;
					assert( order_type == ePerUnorderedGenotype ) ;
					assert( value_type == eProbability ) ;
					// assume diploid
					m_target.set_number_of_entries( 2, ePerUnorderedHaplotype, eAlleleIndex ) ;
				}
				void operator()( MissingValue const value ) {
					m_missing = true ;
					m_entry_i++ ;
					if( m_entry_i == m_calls.size() ) { 
						send_results() ;
					}
				}
				void operator()( std::string& value ) {
					throw BadArgumentError(
						"genfile::ThreshholdingDataReader::Threshholder::operator()",
						"value",
						"Expected a floating-point value (got a string)."
					) ;
				}
				void operator()( Integer const value ) {
					throw BadArgumentError(
						"genfile::ThreshholdingDataReader::Threshholder::operator()",
						"value",
						"Expected a floating-point value (got an integer)."
					) ;
				}
				void operator()( double const value ) {
					m_calls( m_entry_i++ ) = value ;
					if( m_entry_i == m_calls.size() ) { 
						send_results() ;
					}
				}
			private:
				PerSampleSetter& m_target ;
				double const m_threshhold ;
				int m_entry_i ;
				Eigen::VectorXd m_calls ;
				bool m_missing ;
				
			private:
				void send_results() {
					if( m_missing || m_calls.array().maxCoeff() < m_threshhold ) {
						m_target( genfile::MissingValue() ) ;
						m_target( genfile::MissingValue() ) ;
					} else {
						if( m_calls(0) >= m_threshhold ) {
							m_target( Integer( 0 ) ) ;
							m_target( Integer( 0 ) ) ;
						}
						else if( m_calls(1) >= m_threshhold ) {
							m_target( Integer( 0 ) ) ;
							m_target( Integer( 1 ) ) ;
						}
						else if( m_calls(2) >= m_threshhold ) {
							m_target( Integer( 1 ) ) ;
							m_target( Integer( 1 ) ) ;
						}
						else {
							m_target( genfile::MissingValue() ) ;
							m_target( genfile::MissingValue() ) ;
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

