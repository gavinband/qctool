
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SampleMappingSNPDataSource.hpp"

namespace genfile {
	namespace impl {
		struct SampleMapping {
			typedef std::auto_ptr< SampleMapping > UniquePtr ;
			static UniquePtr create(
				CohortIndividualSource const& source_samples,
				std::string const& source_sample_column,
				CohortIndividualSource const& target_samples,
				std::string const& target_sample_column
			) ;
			
			SampleMapping(
				CohortIndividualSource const& reference_samples,
				std::string const& reference_sample_column,
				CohortIndividualSource const& samples,
				std::string const& mapped_sample_column
			) ;
			
			std::size_t number_of_source_samples() const ;
			std::size_t number_of_target_samples() const ;
			boost::optional< std::size_t > find_source_sample_for( std::size_t n ) const ;
			boost::optional< std::size_t > find_target_sample_for( std::size_t n ) const ;

		private:
			
		} ;
		
		struct SampleMappingPerSampleSetter: public VariantDataReader::PerSampleSetter {
			SampleMappingPerSampleSetter(
				VariantDataReader::PerSampleSetter& setter,
				SampleMapping const& mapped_to_reference_sample_mapping
			):
					m_setter( setter ),
					m_sample_mapping( mapped_to_reference_sample_mapping ),
					m_set_this_sample( false )
			{}
			
			~SampleMappingPerSampleSetter() throw() {} ;
			
			void set_number_of_samples( std::size_t n ) {
				m_setter.set_number_of_samples( m_sample_mapping.number_of_source_samples() ) ;
			}

			void set_sample( std::size_t n ) {
				boost::optional< std::size_t > const source_sample = m_sample_mapping.find_source_sample_for( n ) ;
				if( source_sample ) {
					m_setter.set_sample( source_sample.get() ) ;
					m_set_this_sample = true ;
				} else {
					m_set_this_sample = false ;
				}
			}

			void set_number_of_entries( std::size_t n ) {
				if( m_set_this_sample ) {
					m_setter.set_number_of_entries( n ) ;
				}
			}

			void operator()( MissingValue const value ) {
				if( m_set_this_sample ) {
					m_setter( value ) ;
				}
			}
			
			void operator()( std::string& value ) {
				if( m_set_this_sample ) {
					m_setter( value ) ;
				}
			}

			void operator()( VariantEntry::Integer const value ) {
				if( m_set_this_sample ) {
					m_setter( value ) ;
				}
			}
			void operator()( double const value ) {
				if( m_set_this_sample ) {
					m_setter( value ) ;
				}
			}

		private:
			VariantDataReader::PerSampleSetter& m_setter ;
			SampleMapping const& m_sample_mapping ;
			bool m_set_this_sample ;
		} ;
		
		struct SampleMappingVariantDataReader: public VariantDataReader {
		public:
			static UniquePtr create(
				VariantDataReader::UniquePtr data_reader,
				SampleMapping const& sample_mapping
			) {
				return UniquePtr(
					new SampleMappingVariantDataReader( data_reader, sample_mapping )
				) ;
			}

			SampleMappingVariantDataReader(
				VariantDataReader::UniquePtr data_reader,
				SampleMapping const& sample_mapping
			):
				m_data_reader( data_reader ),
				m_sample_mapping( sample_mapping )
			{}

		public:
			VariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				SampleMappingPerSampleSetter mapping_setter( setter, m_sample_mapping ) ;
				m_data_reader->get(
					spec,
					mapping_setter
				) ;
				return *this ;
			}
			
			std::size_t get_number_of_samples() const {
				return m_sample_mapping.number_of_source_samples() ;
			}

			bool supports( std::string const& spec ) const {
				return m_data_reader->supports( spec ) ;
			}

			void get_supported_specs( SpecSetter setter ) const {
				return m_data_reader->get_supported_specs( setter ) ;
			}

		private:
			VariantDataReader::UniquePtr m_data_reader ;
			SampleMapping const& m_sample_mapping ;
		} ;
	}

	SampleMappingSNPDataSource::SampleMappingSNPDataSource(
		CohortIndividualSource const& reference_samples,
		std::string const& reference_sample_column,
		CohortIndividualSource const& samples,
		std::string const& mapped_sample_column,
		SNPDataSource::UniquePtr source
	):
		m_sample_mapping(
			impl::SampleMapping::create(
				reference_samples, reference_sample_column,
				samples, mapped_sample_column
			)
		),
		m_source( source )
	{
	}

	SampleMappingSNPDataSource::operator bool() const {
		return (*m_source) ;
	}

	unsigned int SampleMappingSNPDataSource::number_of_samples() const {
		return m_sample_mapping->number_of_source_samples() ;
	}

	SNPDataSource::OptionalSnpCount SampleMappingSNPDataSource::total_number_of_snps() const {
		return m_source->total_number_of_snps() ;
	}

	std::string SampleMappingSNPDataSource::get_source_spec() const {
		return "SampleMappingSNPDataSource(" + m_source->get_source_spec() + ")" ;
	}

	SNPDataSource const& SampleMappingSNPDataSource::get_parent_source() const {
		return *m_source ;
	}

	void SampleMappingSNPDataSource::get_snp_identifying_data_impl( 
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

	VariantDataReader::UniquePtr SampleMappingSNPDataSource::read_variant_data_impl() {
		return impl::SampleMappingVariantDataReader::create(
			m_source->read_variant_data(),
			*m_sample_mapping
		) ;
	}

	void SampleMappingSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}
	
	void SampleMappingSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}
}
