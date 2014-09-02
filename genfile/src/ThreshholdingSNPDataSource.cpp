
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <queue>
#include <memory>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/ptr_container/ptr_deque.hpp>
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

	unsigned int ThreshholdingSNPDataSource::number_of_samples() const {
		return m_source->number_of_samples() ;
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
				{}
					
				void set_number_of_samples( std::size_t n ) {
					m_target.set_number_of_samples( n ) ;
				} ;
				void set_sample( std::size_t i ) {
					m_target.set_sample( i ) ;
				}
				void set_number_of_entries( std::size_t n ) {
					m_target.set_number_of_entries( n ) ;
				}
				void set_order_type( OrderType const type ) {
					m_target.set_order_type( type ) ;
				}
				void operator()( MissingValue const value ) {
					m_target( value ) ;
				}
				void operator()( std::string& value ) {
					m_target( value ) ;
				}
				void operator()( Integer const value ) {
					m_target( value ) ;
				}
				void operator()( double const value ) {
					m_target( ( value >= m_threshhold) ? 1.0 : 0.0 ) ;
				}
			private:
				PerSampleSetter& m_target ;
				double const m_threshhold ;
			} ;

			VariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				if( spec != ":genotypes:" ) {
					m_source->get( spec, setter ) ;
				} else {
					Threshholder threshholder( setter, m_threshhold ) ;
					m_source->get( spec, threshholder ) ;
				}
				return *this ;
			}

			bool supports( std::string const& spec ) const {
				return m_source->supports( spec ) ;
			}
			void get_supported_specs( SpecSetter setter ) const {
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

