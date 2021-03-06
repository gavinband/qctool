
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/vector_of.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap.hpp>
#include "genfile/Error.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/SampleMappingSNPDataSource.hpp"

namespace genfile {
	namespace impl {
		struct SampleMapping {
			typedef std::auto_ptr< SampleMapping > UniquePtr ;
			typedef boost::bimaps::bimap<
				boost::bimaps::unordered_set_of< std::size_t >,
				boost::bimaps::unordered_set_of< std::size_t >
			> Map ;

		public:
			static UniquePtr create(
				CohortIndividualSource const& source_samples,
				std::string const& source_sample_column,
				CohortIndividualSource const& target_samples,
				std::string const& target_sample_column
			) ;
			
			SampleMapping(
				CohortIndividualSource const& source_samples,
				std::string const& source_sample_column,
				CohortIndividualSource const& target_samples,
				std::string const& target_sample_column
			) {
				setup( source_samples, source_sample_column, target_samples, target_sample_column ) ;
			}
			
			std::size_t number_of_source_samples() const {
				return m_source_ids.size() ;
			}

			std::size_t number_of_target_samples() const {
				return m_target_ids.size() ;
			}

			boost::optional< std::size_t > find_source_sample_for( std::size_t n ) const {
				typedef Map::right_map::const_iterator right_const_iterator ;
				boost::optional< std::size_t > result ;
				right_const_iterator where = m_map.right.find( n ) ;
				if( where != m_map.right.end() ) {
					result = where->second ;
				}
				return result ;
			}

			boost::optional< std::size_t > find_target_sample_for( std::size_t n ) const {
				typedef Map::left_map::const_iterator left_const_iterator ;
				boost::optional< std::size_t > result ;
				left_const_iterator where = m_map.left.find( n ) ;
				if( where != m_map.left.end() ) {
					result = where->second ;
				}
				return result ;
			}

		private:
			typedef std::vector< VariantEntry > SampleIdList ;
			SampleIdList m_source_ids ;
			SampleIdList m_target_ids ;
			Map m_map ;
		private:
			void setup(
				CohortIndividualSource const& source_samples,
				std::string const& source_sample_column,
				CohortIndividualSource const& target_samples,
				std::string const& target_sample_column
			) {
				//std::cerr << "setup()\n" ;
				void(SampleIdList::*push_back)(VariantEntry const&) = &SampleIdList::push_back;	
				source_samples.get_column_values( source_sample_column, boost::bind( push_back, &m_source_ids, _2 ) ) ;
				target_samples.get_column_values( target_sample_column, boost::bind( push_back, &m_target_ids, _2 ) ) ;
				// typedef boost::unordered_map< VariantEntry, std::size_t > TargetMap ;
				typedef std::map< VariantEntry, std::size_t > TargetMap ;
				TargetMap target_map ;
				for( std::size_t i = 0; i < m_target_ids.size(); ++i ) {
					if( !target_map.insert( std::make_pair( m_target_ids[i], i )).second ) {
						throw BadArgumentError(
							"genfile::impl::SampleMapping::setup()",
							"target_sample_column=\"" + target_sample_column + "\"",
							"Target sample " + source_sample_column + "=\"" + m_target_ids[i].as< std::string >() + "\" occurs more than once."
						) ;
					} ;
				}
				for( std::size_t i = 0; i < m_source_ids.size(); ++i ) {
				//	std::cerr << i << ".. " ;
					TargetMap::const_iterator where = target_map.find( m_source_ids[i] ) ;
					if( where != target_map.end() ) {
						// for debug purposes
						//std::vector< VariantEntry >::const_iterator whereold = std::find( m_target_ids.begin(), m_target_ids.end(), m_source_ids[i] ) ;
						//assert(( whereold - m_target_ids.begin() ) == where->second ) ;
						if( !m_map.insert( Map::value_type( i, where->second ) ).second ) {
							throw BadArgumentError(
								"genfile::impl::SampleMapping::setup()",
								"target_sample_column=\"" + target_sample_column + "\"",
								"Target sample \"" + where->first.as< std::string >() + "\" matches multiple source samples."
							) ;
						}
					}
				}
				//std::cerr << "end setup()\n" ;
			}
		} ;
		
		SampleMapping::UniquePtr SampleMapping::create(
			CohortIndividualSource const& source_samples,
			std::string const& source_sample_column,
			CohortIndividualSource const& target_samples,
			std::string const& target_sample_column
		) {
			return SampleMapping::UniquePtr(
				new SampleMapping(
					source_samples, source_sample_column,
					target_samples, target_sample_column
				)
			) ;
		}

		struct SampleMappingPerSampleSetter: public VariantDataReader::PerSampleSetter {
			SampleMappingPerSampleSetter(
				VariantDataReader::PerSampleSetter& setter,
				SampleMapping const& mapped_to_reference_sample_mapping
			):
					m_setter( setter ),
					m_sample_mapping( mapped_to_reference_sample_mapping ),
					m_set_this_sample( false ),
					m_last_source_sample( 0 )
			{}
			
			~SampleMappingPerSampleSetter() throw() {} ;
			
			void initialise( std::size_t nSamples, std::size_t nAlleles ) {
				m_setter.initialise( m_sample_mapping.number_of_source_samples(), nAlleles ) ;
				m_last_source_sample = 0 ;
			}

			bool set_sample( std::size_t n ) {
				boost::optional< std::size_t > const source_sample = m_sample_mapping.find_source_sample_for( n ) ;
				if( source_sample ) {
					for( int s = m_last_source_sample; s < source_sample.get(); ++s ) {
						m_setter.set_sample( s ) ;
					}
					m_setter.set_sample( source_sample.get() ) ;
					m_set_this_sample = true ;
					m_last_source_sample = source_sample.get() + 1 ;
				} else {
					m_set_this_sample = false ;
				}
				return m_set_this_sample ;
			}

			void set_number_of_entries( uint32_t ploidy, std::size_t n, OrderType const order_type, ValueType const value_type ) {
				if( m_set_this_sample ) {
					m_setter.set_number_of_entries( ploidy, n, order_type, value_type ) ;
				}
			}

			void set_value( std::size_t i, MissingValue const value ) {
				if( m_set_this_sample ) {
					m_setter.set_value( i, value ) ;
				}
			}
			
			void set_value( std::size_t i, std::string& value ) {
				if( m_set_this_sample ) {
					m_setter.set_value( i, value ) ;
				}
			}

			void set_value( std::size_t i, VariantEntry::Integer const value ) {
				if( m_set_this_sample ) {
					m_setter.set_value( i, value ) ;
				}
			}
			void set_value( std::size_t i, double const value ) {
				if( m_set_this_sample ) {
					m_setter.set_value( i, value ) ;
				}
			}
			
			void finalise() {
				m_setter.finalise() ;
			}

		private:
			VariantDataReader::PerSampleSetter& m_setter ;
			SampleMapping const& m_sample_mapping ;
			bool m_set_this_sample ;
			std::size_t m_last_source_sample ;
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

	SNPDataSource::Metadata SampleMappingSNPDataSource::get_metadata() const {
		return m_source->get_metadata() ;
	}

	unsigned int SampleMappingSNPDataSource::number_of_samples() const {
		return m_sample_mapping->number_of_source_samples() ;
	}

	void SampleMappingSNPDataSource::get_sample_ids( GetSampleIds ) const {
		assert(0) ; // requires implementation
		return ;	
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

	void SampleMappingSNPDataSource::get_snp_identifying_data_impl( VariantIdentifyingData* variant	) {
		m_source->get_snp_identifying_data( variant ) ;
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
