
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <set>
#include <boost/bind.hpp>
#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceRack.hpp"
#include "genfile/VariantIdentifyingData.hpp"
#include "genfile/OffsetFlippedAlleleSetter.hpp"
#include "genfile/get_set.hpp"
#include "genfile/get_list_of_snps_in_source.hpp"
#include "genfile/string_utils/slice.hpp"

// #define DEBUG_SNP_DATA_SOURCE_RACK 1

namespace genfile {
	namespace {
		static char const eUnknownFlip = OffsetFlippedAlleleSetter::eUnknownFlip ;
		static char const eNoFlip = OffsetFlippedAlleleSetter::eNoFlip ;
		static char const eFlip  = OffsetFlippedAlleleSetter::eFlip ;
	}

	std::auto_ptr< SNPDataSourceRack > SNPDataSourceRack::create( std::vector< wildcard::FilenameMatch > const& filenames ) {
		std::auto_ptr< SNPDataSourceRack > rack( new SNPDataSourceRack() ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			rack->add_source( SNPDataSource::create( filenames[i].filename(), filenames[i].match() )) ;
		}
		return rack ;
	}

	std::auto_ptr< SNPDataSourceRack > SNPDataSourceRack::create(
		std::vector< wildcard::FilenameMatch > const& filenames,
		std::string const& snp_match_fields
	) {
		std::auto_ptr< SNPDataSourceRack > rack( new SNPDataSourceRack( snp_match_fields ) ) ;
		for( std::size_t i = 0; i < filenames.size(); ++i ) {
			rack->add_source( SNPDataSource::create( filenames[i].filename(), filenames[i].match() )) ;
		}
		return rack ;
	}
	
	SNPDataSourceRack::SNPDataSourceRack()
		: m_number_of_samples(0),
		  m_read_past_end( false ),
		  m_comparator( "position,rsid,SNPID,alleles" )
	{
	}

	SNPDataSourceRack::SNPDataSourceRack( std::string const& snp_match_fields )
		: m_number_of_samples(0),
		  m_read_past_end( false ),
		  m_comparator( snp_match_fields )
	{
		// First match must be on position, since this is implicitly assumed below.
		assert( snp_match_fields.find( "position" ) == 0 ) ;
	}

	SNPDataSourceRack::SNPDataSourceRack( genfile::VariantIdentifyingData::CompareFields const& comparator )
		: m_number_of_samples(0),
		  m_read_past_end( false ),
		  m_comparator( comparator )
	{
		assert(
			(m_comparator.get_compared_fields().size() > 0)
			&&
			(m_comparator.get_compared_fields()[0] == VariantIdentifyingData::CompareFields::ePosition)
		) ;
	}

	SNPDataSourceRack::~SNPDataSourceRack() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			delete m_sources[i] ;
		}
	}

	void SNPDataSourceRack::add_source(
		std::auto_ptr< SNPDataSource > source
	) {
		m_sources.push_back( source.release() ) ;
		m_flips.push_back( eNoFlip ) ;
		m_number_of_samples += m_sources.back()->number_of_samples() ;
		m_sources.back()->reset_to_start() ;
		// update metadata
		if( m_sources.size() == 1 ) {
			m_metadata = m_sources.back()->get_metadata() ;
		} else {
			// We just take the first source's metadata, and assume the user knows what they're doing.
		}
	}
	
	void SNPDataSourceRack::check_snps_are_sorted_by_position(
		std::vector< VariantIdentifyingData > const& snps,
		std::size_t cohort_index
	) {
		for( std::size_t i = 1; i < snps.size(); ++i ) {
			if( snps[i].get_position() < snps[i-1].get_position() ) {
				throw BadArgumentError(
					"genfile::SNPDataSourceRack::add_source()",
					"snps in cohort "
						+ string_utils::to_string( cohort_index + 1 )
						+ " must be in nondecreasing order of position, at SNP "
						+ string_utils::to_string( i )
						+ " with position "
						+ string_utils::to_string( snps[i].get_position() )
				) ;
			}
		}
	}

	// Implicit conversion to bool.  Return true if there have been no errors so far.
	SNPDataSourceRack::operator bool() const {
		if( m_sources.empty() ) {
			return 0 ;
		}
		for( std::size_t i = 0; i < m_sources.size(); ++ i ) {
			if( !*m_sources[i] ) {
				return false ;
			}
		}
		return true ;
	}

	SNPDataSource::Metadata SNPDataSourceRack::get_metadata() const {
		return m_metadata ;
	}

	namespace impl {
		struct OffsetSetter {
		public:
			OffsetSetter( SNPDataSource::GetSampleIds setter, std::size_t offset ):
				m_setter( setter ),
				m_offset( offset )
			{}

			void operator()( std::size_t i, std::string const& value ) const {
				m_setter( i + m_offset, value ) ;
			}
		private:
			SNPDataSource::GetSampleIds m_setter ;
			std::size_t m_offset ;
		} ;
	}
	void SNPDataSourceRack::get_sample_ids( GetSampleIds getter ) const {
		std::size_t offset = 0 ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			impl::OffsetSetter setter( getter, offset ) ;
			m_sources[i]->get_sample_ids( setter ) ;
			offset += m_sources[i]->number_of_samples() ;
		}
	}
	
	std::string SNPDataSourceRack::get_source_spec() const {
		std::string result = "rack:" ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			result += "," + m_sources[i]->get_source_spec() ;
		}
		return result ;
	}
	
	std::string SNPDataSourceRack::get_summary( std::string const& prefix, std::size_t width ) const {
		std::ostringstream ostr ;
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			ostr << prefix << std::setw( width ) << "cohort " << (i+1) << ":\n" ;
			ostr << m_sources[i]->get_summary( prefix, width ) ;
			ostr << "\n" ;
		}
		ostr << prefix << std::setw( width ) << "Total all cohorts:" << " " << number_of_samples() << " samples.\n" ;
		return ostr.str() ;
	}

	std::size_t SNPDataSourceRack::number_of_sources() const {
		return m_sources.size() ;
	}

	SNPDataSource& SNPDataSourceRack::get_source( std::size_t index ) const {
		assert( index < m_sources.size() ) ;
		return *m_sources[ index ] ;
	}

	// Return the number of samples represented in the snps in this source.
	unsigned int SNPDataSourceRack::number_of_samples() const {
		return m_number_of_samples ;
	}

	// Return the total number of snps the source contains.
	SNPDataSource::OptionalSnpCount SNPDataSourceRack::total_number_of_snps() const {
		if( m_sources.size() == 0 ) {
			return std::size_t( 0 ) ;
		}
		else if( m_sources.size() == 1 ) {
			return m_sources[0]->total_number_of_snps() ;
		}
		else {
			return OptionalSnpCount() ;
		}
	}

	void SNPDataSourceRack::reset_to_start_impl() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			m_sources[i]->reset_to_start() ;
		}
		m_read_past_end = false ;
	}

	namespace {
		void combine_ids( std::set< std::string >& target, std::string& combined, std::string const& value ) {
			if( target.insert( value ).second ) {
				combined += ( value.size() > 0 ? "," : "" ) + value ;
			}
		}
	}
	
	void SNPDataSourceRack::get_snp_identifying_data_impl( 
		VariantIdentifyingData* result
	) {
		if( m_sources.size() == 0 ) {
			return ;
		}
		
		VariantIdentifyingData this_snp ;
		std::string merged_SNPID, merged_rsid, merged_allele1, merged_allele2 ;
		std::set< std::string > rsids ;
		std::set< std::string > SNPIDs ;
		std::vector< char > flips( m_sources.size(), eNoFlip ) ;
		std::size_t source_i = 0 ;
		if( m_sources[0]->get_snp_identifying_data( &this_snp ) ) {
			// we use source_i == m_source.size() to indicate a successful match across cohorts.
			while( (*this) && source_i < m_sources.size() ) {
				SNPIDs.clear() ;
				rsids.clear() ;

				this_snp.get_alternative_identifiers( boost::bind( &combine_ids, &SNPIDs, &merged_SNPID, _1 ) ) ;
				rsids.insert( this_snp.get_rsid() ) ;
				merged_rsid = this_snp.get_rsid() ;
				merged_allele1 = this_snp.get_first_allele() ;
				merged_allele2 = this_snp.get_second_allele() ;

				VariantIdentifyingData this_source_snp ;
				source_i = 1 ;
				for( ;
					source_i < m_sources.size() && move_source_to_snp_matching( source_i, this_snp, &this_source_snp );
					++source_i
				) {
#if DEBUG_SNP_DATA_SOURCE_RACK
					std::cerr << "SNP " << this_snp << " is in source " << source_i << ".\n" ;
#endif
					this_source_snp.get_alternative_identifiers( boost::bind( &combine_ids, &SNPIDs, &merged_SNPID, _1 ) ) ;
					if( rsids.insert( this_source_snp.get_rsid() ).second ) {
						merged_rsid += "," + std::string( this_source_snp.get_rsid() ) ;
					}
					if( m_comparator.get_flip_alleles_if_necessary() ) {
						if( this_source_snp.get_first_allele() == this_snp.get_second_allele() && this_source_snp.get_second_allele() == this_snp.get_first_allele() ) {
							flips[ source_i ] = eFlip ;
						} else if( this_source_snp.get_first_allele() == this_snp.get_first_allele() && this_source_snp.get_second_allele() == this_snp.get_second_allele() ) {
							// do nothing
						} else {
							merged_allele1 += "/" + std::string( this_source_snp.get_first_allele() ) ;
							merged_allele2 += "/" + std::string( this_source_snp.get_second_allele() ) ;
							// Don't know how to flip.
							flips[ source_i ] = eUnknownFlip ;
						}
					} else {
						if( this_source_snp.get_first_allele() != this_snp.get_first_allele() || this_source_snp.get_second_allele() != this_snp.get_second_allele() ) {
							merged_allele1 += "/" + std::string( this_source_snp.get_first_allele() ) ;
							merged_allele2 += "/" + std::string( this_source_snp.get_second_allele() ) ;
							// Set genotypes to missing for this cohort.
							flips[ source_i ] = eUnknownFlip ;
						}
					}
				}
				if( (*this) && source_i < m_sources.size() ) {
					m_sources[0]->ignore_snp_probability_data() ;
					m_sources[0]->get_snp_identifying_data( &this_snp ) ;
				}
			}
		}
		if( *this ) {
			*result = VariantIdentifyingData( merged_SNPID, merged_rsid, this_snp.get_position(), merged_allele1, merged_allele2 ) ;
			// Remember the flips.
			m_flips = flips ;
		}
	}

	// Read the identifying info of the next SNP in the given source
	// which has the given chromosome and SNP position.
	// If no such SNP is found, throw MissingSNPError.
	// If such a SNP is found but it has different SNPID, RSID, or alleles, throw SNPMismatchError.
	bool SNPDataSourceRack::move_source_to_snp_matching(
		std::size_t source_i,
		VariantIdentifyingData const& reference_snp,
		VariantIdentifyingData* result
	) {
#if DEBUG_SNP_DATA_SOURCE_RACK
		std::cerr << "Moving source " << source_i << " to SNP " << reference_snp << ".\n" ;
#endif
		VariantIdentifyingData this_snp ;
		
		while(
			m_sources[source_i]->get_next_snp_with_specified_position(
				&this_snp,
				reference_snp.get_position().chromosome(),
				reference_snp.get_position().position()
			)
		) {
#if DEBUG_SNP_DATA_SOURCE_RACK
			std::cerr << "Found candidate: " << this_snp << ".\n" ;
#endif
			if( m_comparator.are_equal( this_snp, reference_snp )) {
				*result = this_snp ;
#if DEBUG_SNP_DATA_SOURCE_RACK
				std::cerr << "Match!\n" ;
#endif
				return true ;
			} else {
#if DEBUG_SNP_DATA_SOURCE_RACK
				std::cerr << "No match!\n" ;
#endif
				m_sources[source_i]->ignore_snp_probability_data() ;
			}
		}
		return false ;
	}
		
	namespace impl {
		void add_spec_to_map( std::map< std::string, std::string >* map, std::string const& name, std::string const& type ) {
			(*map)[ name ] = type ;
		}

		struct RackVariantDataReader: public VariantDataReader
		{
			RackVariantDataReader( SNPDataSourceRack& rack ):
				m_rack( rack )
			{
				for( std::size_t i = 0; i < m_rack.m_sources.size(); ++i ) {
					m_data_readers.push_back( m_rack.m_sources[i]->read_variant_data().release() ) ;
				}
			}

			~RackVariantDataReader() {
				for( std::size_t i = 0; i < m_data_readers.size(); ++i ) {
					delete m_data_readers[i] ;
				}
			}

			RackVariantDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				OffsetFlippedAlleleSetter offset_flip_sample_setter( setter, m_rack.number_of_samples(), eNoFlip, 0 ) ;
				std::size_t sample_offset = 0 ;
				for( std::size_t i = 0; i < m_rack.m_sources.size(); ++i ) {
					offset_flip_sample_setter.set_offset( sample_offset ) ;
					offset_flip_sample_setter.set_flip( m_rack.get_flip( i ) ) ;
					m_data_readers[i]->get( spec, offset_flip_sample_setter ) ;
					sample_offset += m_rack.m_sources[i]->number_of_samples() ;
				}
				return *this ;
			}

			std::size_t get_number_of_samples() const {
				return m_rack.number_of_samples() ;
			}

			// True if all the constituent sources support the spec.
			bool supports( std::string const& spec ) const {
				for( std::size_t i = 0; i < m_data_readers.size(); ++i ) {
					if( !m_data_readers[i]->supports( spec )) {
						return false ;
					}
				}
				return true ;
			}
			
			void get_supported_specs( SpecSetter setter ) const {
				if( m_data_readers.empty() ) {
					return ;
				}

				typedef std::map< std::string, std::string > SetOfThings ;
				SetOfThings things ;
				for( std::size_t i = 0; i < m_data_readers.size(); ++i ) {
					SetOfThings these_things ;
					m_data_readers[0]->get_supported_specs(
						boost::bind(
							&add_spec_to_map,
							&these_things,
							_1,
							_2
						)
					) ;
					if( i == 0 ) {
						things = these_things ;
					}
					else {
						for( SetOfThings::iterator thing_i = things.begin(); thing_i != things.end(); ) {
							SetOfThings::iterator where = these_things.find( thing_i->first ) ;

							if( where == these_things.end() ) {
								SetOfThings::iterator this_thing = thing_i++ ;
								things.erase( this_thing ) ;
							}
							else{
								++thing_i ;
							}
						}
					}
				}
				for( SetOfThings::const_iterator thing_i = things.begin(); thing_i != things.end(); ++thing_i ) {
					setter( thing_i->first, thing_i->second ) ;
				}
			}

			private:
				std::vector< VariantDataReader* > m_data_readers ;
				SNPDataSourceRack& m_rack ;
		} ;
	}

	VariantDataReader::UniquePtr SNPDataSourceRack::read_variant_data_impl() {
		VariantDataReader::UniquePtr result( new impl::RackVariantDataReader( *this ) ) ;
		return result ;
	}

	void SNPDataSourceRack::ignore_snp_probability_data_impl() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			m_sources[i]->ignore_snp_probability_data() ;
		}
	}
	
	SNPDataSourceRack::RackGenotypeProbabilitySetter::RackGenotypeProbabilitySetter(
		GenotypeProbabilitySetter const& base_setter,
		uint32_t index_of_first_sample
	) 
		: m_base_setter( base_setter ),
	  	  m_index_of_first_sample( index_of_first_sample )
	{}
	
	void SNPDataSourceRack::RackGenotypeProbabilitySetter::operator()(
		std::size_t i, double AA, double AB, double BB
	) const {
		m_base_setter( i + m_index_of_first_sample, AA, AB, BB ) ;
	}
}
