
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <sstream>
#include <set>
#include <boost/bind.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SampleFilteringSNPDataSource.hpp"
#include "genfile/get_set.hpp"

namespace genfile {
	// Create a chain of SNPDataSources taking data from the specified files.
	std::auto_ptr< SampleFilteringSNPDataSource > SampleFilteringSNPDataSource::create(
		std::auto_ptr< SNPDataSource > source,
		std::set< std::size_t > const& indices_of_samples_to_filter_out
	) {
		return std::auto_ptr< SampleFilteringSNPDataSource >( new SampleFilteringSNPDataSource( source, indices_of_samples_to_filter_out )) ;
	}
	
	SampleFilteringSNPDataSource::SampleFilteringSNPDataSource(
		std::auto_ptr< SNPDataSource > source,
		std::set< std::size_t > const& indices_of_samples_to_filter_out
	)
		: m_source( source ),
		  m_indices_of_samples_to_filter_out( indices_of_samples_to_filter_out.begin(), indices_of_samples_to_filter_out.end() ),
		  m_genotype_data( m_source->number_of_samples() * 3 )
	{
		std::sort( m_indices_of_samples_to_filter_out.begin(), m_indices_of_samples_to_filter_out.end() ) ;
		if( m_indices_of_samples_to_filter_out.size() > 0 && m_indices_of_samples_to_filter_out.back() >= m_source->number_of_samples() ) {
			throw SampleIndexOutOfRangeError( m_indices_of_samples_to_filter_out.back(), m_source->number_of_samples() ) ;
		}
	}
	
	SampleFilteringSNPDataSource::~SampleFilteringSNPDataSource() {}

	SNPDataSource::Metadata SampleFilteringSNPDataSource::get_metadata() const {
		return m_source->get_metadata() ;
	}

	unsigned int SampleFilteringSNPDataSource::number_of_samples() const {
		return m_source->number_of_samples() - m_indices_of_samples_to_filter_out.size() ;
	}

	bool SampleFilteringSNPDataSource::has_sample_ids() const {
		return m_source->has_sample_ids() ;
	}

	namespace {
		void set_vector_entry(
			std::size_t i,
			std::string const& id,
			std::vector< std::string >* result
		) {
			assert( i < result->size() ) ;
			(*result)[i] = id ;
		}
	}

	void SampleFilteringSNPDataSource::get_sample_ids( GetSampleIds callback ) const {
		std::vector< std::string > ids( number_of_samples() ) ;
		m_source->get_sample_ids(
			boost::bind( &set_vector_entry, _1, _2, &ids )
		) ;
		std::vector< std::size_t >::const_iterator
			i( m_indices_of_samples_to_filter_out.begin() ),
			end_i( m_indices_of_samples_to_filter_out.end() ) ;

		std::size_t index_of_last_unfiltered_sample = 0 ;
		std::size_t index_in_output = 0 ;
		for( ; i != end_i; ++i ) {
			for( std::size_t sample_i = index_of_last_unfiltered_sample; sample_i < *i; ++sample_i ) {
				callback( index_in_output++, ids[sample_i] ) ;
			}
			index_of_last_unfiltered_sample = (*i) + 1 ;
		}
		for( std::size_t sample_i = index_of_last_unfiltered_sample; sample_i < m_source->number_of_samples() ; ++sample_i ) {
			callback( index_in_output++, ids[sample_i] ) ;
		}
	}

	SNPDataSource::OptionalSnpCount SampleFilteringSNPDataSource::total_number_of_snps() const {
		return m_source->total_number_of_snps() ;
	}

	std::string SampleFilteringSNPDataSource::get_source_spec() const {
		return "sample-filtered:" + m_source->get_source_spec() ;
	}

	SNPDataSource const& SampleFilteringSNPDataSource::get_parent_source() const {
		return *m_source ;
	}

	SNPDataSource const& SampleFilteringSNPDataSource::get_base_source() const {
		return m_source->get_base_source() ;
	}

	std::string SampleFilteringSNPDataSource::get_summary( std::string const& prefix, std::size_t width ) const {
		std::ostringstream ostr ;
		ostr << m_source->get_summary( prefix, width ) ;
		ostr << prefix << std::setw( width ) << "Number of samples (post-filter):" << " " << number_of_samples() << "\n" ;
		return ostr.str() ;
	}

	SampleFilteringSNPDataSource::operator bool() const {
		return *m_source ;
	}

	void SampleFilteringSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}

	void SampleFilteringSNPDataSource::get_snp_identifying_data_impl( 
		VariantIdentifyingData* result
	) {
		m_source->get_snp_identifying_data( result ) ;
	}

	void SampleFilteringSNPDataSource::read_source_probability_data() {
		m_source->read_snp_probability_data(
			boost::bind( &SampleFilteringSNPDataSource::set_unfiltered_genotype_probabilities, this, _1, _2, _3, _4 )
		) ;
	}

	void SampleFilteringSNPDataSource::return_filtered_genotype_probabilities(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		std::vector< std::size_t >::const_iterator
			i( m_indices_of_samples_to_filter_out.begin() ),
			end_i( m_indices_of_samples_to_filter_out.end() ) ;

		std::size_t index_of_last_unfiltered_sample = 0 ;
		std::size_t index_of_filtered_sample = 0 ;

		for( ; i != end_i; ++i ) {
			for( std::size_t sample_i = index_of_last_unfiltered_sample; sample_i < *i; ++sample_i ) {
				std::size_t genotype_data_index = sample_i * 3 ;
				set_genotype_probabilities(
					index_of_filtered_sample++,
					m_genotype_data[ genotype_data_index ],
					m_genotype_data[ genotype_data_index + 1 ],
					m_genotype_data[ genotype_data_index + 2 ]
				) ;
			}
			index_of_last_unfiltered_sample = (*i) + 1 ;
		}

		for( std::size_t sample_i = index_of_last_unfiltered_sample; sample_i < m_source->number_of_samples() ; ++sample_i ) {
			std::size_t genotype_data_index = sample_i * 3 ;
			set_genotype_probabilities(
				index_of_filtered_sample++,
				m_genotype_data[ genotype_data_index ],
				m_genotype_data[ genotype_data_index + 1 ],
				m_genotype_data[ genotype_data_index + 2 ]
			) ;
		}
	}
	
	void SampleFilteringSNPDataSource::set_unfiltered_genotype_probabilities( std::size_t i, double aa, double ab, double bb ) {
		i *= 3 ;
		assert( i + 2 < m_genotype_data.size() ) ;
		m_genotype_data[i++] = aa ;
		m_genotype_data[i++] = ab ;
		m_genotype_data[i] = bb ;
	}

	void SampleFilteringSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}

	namespace impl {
		
		struct SampleFilteringPerSampleSetter: public VariantDataReader::PerSampleSetter {
			SampleFilteringPerSampleSetter(
				VariantDataReader::PerSampleSetter& setter,
				std::vector< std::size_t > const& indices_of_samples_to_filter_out
			):
					m_setter( setter ),
					m_indices_of_samples_to_filter_out( indices_of_samples_to_filter_out ),
					m_next_filtered_out( m_indices_of_samples_to_filter_out.size() ),
					m_filter_out_this_sample( false )
			{}
			
			~SampleFilteringPerSampleSetter() throw() {} ;
			
			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				m_setter.initialise( nSamples - m_indices_of_samples_to_filter_out.size(), nAlleles ) ;
				m_number_filtered_out = 0 ;
				m_next_filtered_out = ( m_indices_of_samples_to_filter_out.size() > 0 ) ? m_indices_of_samples_to_filter_out.front() : nSamples ;
				m_filter_out_this_sample = false ;
			}

			bool set_sample( std::size_t n ) {
			 	if( n == m_next_filtered_out ) {
					m_filter_out_this_sample = true ;
					++m_number_filtered_out ;
					m_next_filtered_out = ( m_number_filtered_out < m_indices_of_samples_to_filter_out.size() ) ? m_indices_of_samples_to_filter_out[ m_number_filtered_out ] : 0 ;
					return false ;
				} else {
					m_filter_out_this_sample = false ;
					m_setter.set_sample( n - m_number_filtered_out ) ;
					return true ;
				}
			}

			void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
				if( !m_filter_out_this_sample ) {
					m_setter.set_number_of_entries( ploidy,  n, order_type, value_type ) ;
				}
			}

			void set_value( std::size_t i, MissingValue const value ) { if( !m_filter_out_this_sample ) { m_setter.set_value( i, value ) ; } }
			void set_value( std::size_t i, std::string& value ) { if( !m_filter_out_this_sample ) { m_setter.set_value( i, value ) ; } }
			void set_value( std::size_t i, VariantEntry::Integer const value ) { if( !m_filter_out_this_sample ) { m_setter.set_value( i, value ) ; } }
			void set_value( std::size_t i, double const value ) { if( !m_filter_out_this_sample ) { m_setter.set_value( i, value ) ; } }

			void finalise() {
				m_setter.finalise() ;
			}

		private:
			VariantDataReader::PerSampleSetter& m_setter ;
			std::vector< std::size_t > const& m_indices_of_samples_to_filter_out ;
			std::size_t m_number_filtered_out ;
			std::size_t m_next_filtered_out ;
			bool m_filter_out_this_sample ;
		} ;
		
		struct SampleFilteringVariantDataReader: public VariantDataReader {
			SampleFilteringVariantDataReader(
				VariantDataReader::UniquePtr data_reader,
				std::vector< std::size_t > const& indices_of_samples_to_filter_out
			):
				m_data_reader( data_reader ),
				m_indices_of_samples_to_filter_out( indices_of_samples_to_filter_out )
			{}

			VariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				SampleFilteringPerSampleSetter filtering_setter( setter, m_indices_of_samples_to_filter_out ) ;
				m_data_reader->get(
					spec,
					filtering_setter
				) ;
				return *this ;
			}
			
			std::size_t get_number_of_samples() const { return m_data_reader->get_number_of_samples() - m_indices_of_samples_to_filter_out.size() ; }

			bool supports( std::string const& spec ) const {
				return m_data_reader->supports( spec ) ;
			}

			void get_supported_specs( SpecSetter setter ) const {
				return m_data_reader->get_supported_specs( setter ) ;
			}

		private:
			VariantDataReader::UniquePtr m_data_reader ;
			std::vector< std::size_t > const& m_indices_of_samples_to_filter_out ;
		} ;
	}

	VariantDataReader::UniquePtr SampleFilteringSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr(
			new impl::SampleFilteringVariantDataReader(
				m_source->read_variant_data(),
				m_indices_of_samples_to_filter_out
			)
		) ;
	}
}

