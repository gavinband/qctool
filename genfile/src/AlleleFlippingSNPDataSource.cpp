
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/bind.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/AlleleFlippingSNPDataSource.hpp"
#include "genfile/get_set.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	std::pair< std::vector< SNPIdentifyingData >, AlleleFlippingSNPDataSource::AlleleFlipSpec >
	AlleleFlippingSNPDataSource::get_allele_flip_spec(
		std::vector< SNPIdentifyingData > reference_snps,
		std::vector< SNPIdentifyingData > snps_to_match,
		SNPIdentifyingData::CompareFields const& comparator
	) {
		std::sort( reference_snps.begin(), reference_snps.end(), comparator ) ;
		AlleleFlipSpec allele_flips( snps_to_match.size() ) ;
		
		for( std::size_t i = 0; i < snps_to_match.size(); ++i ) {
			SNPIdentifyingData swapped_snp = snps_to_match[i] ;
			std::swap( swapped_snp.first_allele(), swapped_snp.second_allele() ) ;
			if( snps_to_match[i].first_allele() == "?" || snps_to_match[i].second_allele() == "?" ) {
				allele_flips[i] = eUnknownFlip ;
			}
			else if(
				std::binary_search(
					reference_snps.begin(),
					reference_snps.end(),
					snps_to_match[i],
					comparator
				)
			) {
				allele_flips[i] = eNoFlip ;
			}
			else if(
				std::binary_search(
					reference_snps.begin(),
					reference_snps.end(),
					swapped_snp,
					comparator
				)
			) {
				snps_to_match[i] = swapped_snp ;
				allele_flips[i] = eFlip ;
			}
			else {
				snps_to_match[i].first_allele() = '?' ;
				snps_to_match[i].second_allele() = '?' ;
				allele_flips[i] = eUnknownFlip ;
			}
		}
		return std::make_pair( snps_to_match, allele_flips ) ;
	}
	
	AlleleFlippingSNPDataSource::UniquePtr AlleleFlippingSNPDataSource::create(
		SNPDataSource::UniquePtr source,
		AlleleFlipSpec const& allele_flips
	) {
		return UniquePtr( new AlleleFlippingSNPDataSource( source, allele_flips )) ;
	}
	
	AlleleFlippingSNPDataSource::AlleleFlippingSNPDataSource(
		SNPDataSource::UniquePtr source,
		AlleleFlipSpec const& allele_flips
	):
		m_source( source ),
		m_allele_flips( allele_flips )
	{
		assert( !m_source->total_number_of_snps() || m_allele_flips.size() == *m_source->total_number_of_snps() ) ;
	}
	
	std::string AlleleFlippingSNPDataSource::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return m_source->get_summary( prefix + "  " ) ;
	}

	void AlleleFlippingSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		char allele_flip = eUnknownFlip ;
		if( number_of_snps_read() < m_allele_flips.size() ) {
			allele_flip = m_allele_flips[ number_of_snps_read() ] ;
		}
		switch( allele_flip ) {
			case eNoFlip:
				m_source->get_snp_identifying_data(
					set_number_of_samples,
					set_SNPID,
					set_RSID,
					set_chromosome,
					set_SNP_position,
					set_allele1,
					set_allele2
				) ;
				break ;
			case eFlip:
				m_source->get_snp_identifying_data(
					set_number_of_samples,
					set_SNPID,
					set_RSID,
					set_chromosome,
					set_SNP_position,
					set_allele2,
					set_allele1
				) ;
				break ;
			case eUnknownFlip:
				m_source->get_snp_identifying_data(
					set_number_of_samples,
					set_SNPID,
					set_RSID,
					set_chromosome,
					set_SNP_position,
					ignore(),
					ignore()
				) ;
				set_allele1( "?" ) ;
				set_allele2( "?" ) ;
				break ;
			default:
				assert(0) ;
		}
	}

	namespace impl {
		struct FlippedAlleleSetter: public VariantDataReader::PerSampleSetter {
			~FlippedAlleleSetter() throw() {}
			
			FlippedAlleleSetter( VariantDataReader::PerSampleSetter& setter ):
				m_setter( setter ),
				m_values( 3 ),
				m_number_of_entries( 0 ),
				m_entry_i( 0 )
			{}
			void set_number_of_samples( std::size_t n ) { m_setter.set_number_of_samples( n ) ; }
			void set_sample( std::size_t n ) { m_setter.set_sample( n ) ; }
			void set_number_of_entries( std::size_t n ) {
				m_number_of_entries = n ;
				m_entry_i = 0 ;
				m_setter.set_number_of_entries( n ) ;
			}

			void operator()( MissingValue const value ) { store( value ) ; }
			void operator()( std::string& value ) { store( value ) ; }
			void operator()( Integer const value ) { store( value ) ; }
			void operator()( double const value ) { store( value ) ; }

		private:
			VariantDataReader::PerSampleSetter& m_setter ;
			std::vector< VariantEntry > m_values ;
			std::size_t m_number_of_entries ;
			std::size_t m_entry_i ;

			template< typename T >
			void store( T value ) {
				if( m_values.size() < ( m_entry_i + 1 ) ) {
					m_values.resize( m_entry_i + 1 ) ;
				}
				m_values[ m_entry_i++ ] = value ;
				if( m_entry_i == m_number_of_entries ) {
					set_values() ;
				}
			}

			void set_values() {
				for( std::size_t i = 0; i < m_number_of_entries; ++i ) {
					VariantEntry const& entry = m_values[ m_number_of_entries - 1 - i ] ;
					if( entry.is_missing() ) {
						m_setter( MissingValue() ) ;
					} else if( entry.is_string() ) {
						std::string value = entry.as< std::string >() ;
						m_setter( value ) ;
					} else if( entry.is_int() ) {
						m_setter( entry.as< VariantEntry::Integer >() ) ;
					} else if( entry.is_double() ) {
						m_setter( entry.as< double >() ) ;
					} else {
						assert(0) ;
					}
				}
			}
		} ;
		
		struct UnknownAlleleSetter: public VariantDataReader::PerSampleSetter {
			UnknownAlleleSetter( VariantDataReader::PerSampleSetter& setter ):
				m_setter( setter )
			{}
			
			void set_number_of_samples( std::size_t n ) { m_setter.set_number_of_samples( n ) ; }
			void set_sample( std::size_t n ) { m_setter.set_sample( n ) ; }
			void set_number_of_entries( std::size_t n ) { m_setter.set_number_of_entries( n ) ; }
			void operator()( MissingValue const value ) { m_setter( value ) ; }
			void operator()( Integer const value ) { m_setter( MissingValue() ) ; }
			void operator()( double const value ) { m_setter( MissingValue() ) ; }
		private:
			VariantDataReader::PerSampleSetter& m_setter ;
		} ;

		class AlleleFlippingSNPDataReader: public VariantDataReader {
		public:
			AlleleFlippingSNPDataReader(
				AlleleFlippingSNPDataSource& source,
				VariantDataReader::UniquePtr base_reader,
				AlleleFlippingSNPDataSource::AlleleFlipSpec const& allele_flips
			):
				m_source( source ),
				m_base_reader( base_reader ),
				m_allele_flips( allele_flips )
			{}
			
			std::size_t get_number_of_samples() const { return m_source.number_of_samples() ; }
			
			AlleleFlippingSNPDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				assert( m_source.number_of_snps_read() > 0 ) ;
				switch( m_allele_flips[ m_source.number_of_snps_read() - 1 ] ) {
					case AlleleFlippingSNPDataSource::eNoFlip:
						m_base_reader->get( spec, setter ) ;
						break ;
					case AlleleFlippingSNPDataSource::eFlip:
						if( spec == "genotypes" || spec == "intensities" ) {
							FlippedAlleleSetter flipped_setter( setter ) ;
							m_base_reader->get( spec, flipped_setter ) ;
						}
						else {
							// Pass through to base reader.
							m_base_reader->get( spec, setter ) ;
						}
						break ;
					case AlleleFlippingSNPDataSource::eUnknownFlip:
						if( spec == "genotypes" || spec == "intensities" ) {
							UnknownAlleleSetter unknown_allele_setter( setter ) ;
							m_base_reader->get( spec, unknown_allele_setter ) ;
						}
						else {
							m_base_reader->get( spec, setter ) ;
						}
						break ; 
					default:
						assert(0) ;
				}
				return *this ;
			}
			
			bool supports( std::string const& spec ) const {
				return m_base_reader->supports( spec ) ;
			}
			
			void get_supported_specs( SpecSetter setter ) const {
				return m_base_reader->get_supported_specs( setter ) ;
			}

		private:
			AlleleFlippingSNPDataSource& m_source ;
			VariantDataReader::UniquePtr m_base_reader ;
			AlleleFlippingSNPDataSource::AlleleFlipSpec const& m_allele_flips ;
		} ;
	}

	VariantDataReader::UniquePtr AlleleFlippingSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr(
			new impl::AlleleFlippingSNPDataReader( *this, m_source->read_variant_data(), m_allele_flips )
		) ;
	}
	
	AlleleFlippingSNPDataSource::AlleleFlippingGenotypeProbabilitySetter::AlleleFlippingGenotypeProbabilitySetter( GenotypeProbabilitySetter setter ):
		m_setter( setter )
	{}

	void AlleleFlippingSNPDataSource::AlleleFlippingGenotypeProbabilitySetter::operator()( std::size_t i, double AA, double AB, double BB ) const {
		m_setter( i, BB, AB, AA ) ;
	}

	void AlleleFlippingSNPDataSource::ignore_snp_probability_data_impl() {
		m_source->ignore_snp_probability_data() ;
	}

	void AlleleFlippingSNPDataSource::reset_to_start_impl() {
		m_source->reset_to_start() ;
	}
}
