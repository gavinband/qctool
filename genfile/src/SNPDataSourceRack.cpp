#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <set>
#include <boost/bind.hpp>
#include "genfile/Error.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/SNPDataSourceRack.hpp"
#include "genfile/SNPIdentifyingData.hpp"
#include "genfile/get_set.hpp"
#include "genfile/get_list_of_snps_in_source.hpp"

namespace genfile {
	
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

	SNPDataSourceRack::~SNPDataSourceRack() {
		for( std::size_t i = 0; i < m_sources.size(); ++i ) {
			delete m_sources[i] ;
		}
	}

	void SNPDataSourceRack::add_source(
		std::auto_ptr< SNPDataSource > source
	) {
		m_sources.push_back( source.release() ) ;
		m_number_of_samples += m_sources.back()->number_of_samples() ;
		m_sources.back()->reset_to_start() ;
	}

	void SNPDataSourceRack::check_snps_are_sorted_by_position(
		std::vector< SNPIdentifyingData > const& snps,
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

	std::vector< SNPIdentifyingData > SNPDataSourceRack::get_intersected_snps(
		std::vector< SNPIdentifyingData > const& snps1,
		std::vector< SNPIdentifyingData > const& snps2
	) const {
		//
		// The algorithm here must match that for get_snp_identifying_data_impl.
		// This has the following feature: each list can only be traversed once,
		// forwards.
		// For each snp in the first list we find the first SNP not already considered in
		// the second list which (has the same position and) matches according to our comparator.
		// 
		std::vector< SNPIdentifyingData >::const_iterator
			i1 = snps1.begin(),
			i1_end = snps1.end(),
			i2 = snps2.begin(),
			i2_end = snps2.end() ;

		SNPIdentifyingData::CompareFields position_comparator( "position" ) ;

		std::vector< SNPIdentifyingData > intersected_snps ;

		for( ; i1 != i1_end && i2 != i2_end ; ++i1 ) {
			// Find next SNP with matching position.
			std::vector< SNPIdentifyingData >::const_iterator
				i2_pos = std::lower_bound( i2, i2_end, *i1, position_comparator ) ;
			// Find next SNP with all fields matching, including position.
			for( ; i2_pos != i2_end && !m_comparator.are_equal( *i2_pos, *i1 ) && i2_pos->get_position() == i1->get_position(); ++i2_pos ) ;

			if( i2_pos != i2_end && i2_pos->get_position() == i1->get_position() ) {
				intersected_snps.push_back( *i1 ) ;
				++i2_pos ;
			}
			i2 = i2_pos ;
		}

		return intersected_snps ;
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

	void SNPDataSourceRack::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		if( m_sources.size() == 0 ) {
			return ;
		}
		
		SNPIdentifyingData this_snp ;
		std::size_t last_matching_source = 0 ;
		while( m_sources[0]->get_snp_identifying_data( this_snp ) && last_matching_source < m_sources.size() ) {
			// We will report ?'s in any field that differs between cohorts.
			SNPIdentifyingData this_source_snp ;
			last_matching_source = 1 ;
			for( ;
				last_matching_source < m_sources.size() && move_source_to_snp_matching( last_matching_source, this_snp, &this_source_snp );
				++last_matching_source
			) {
				if( this_source_snp.get_SNPID() != this_snp.get_SNPID() ) {
					this_snp.SNPID() += "/" + this_source_snp.get_SNPID() ;
				}
				if( this_source_snp.get_rsid() != this_snp.get_rsid() ) {
					this_snp.rsid() += "/" + this_source_snp.get_SNPID() ;
				}
				if( this_source_snp.get_first_allele() != this_snp.get_first_allele() ) {
					this_snp.first_allele() = '?' ;
				}
				if( this_source_snp.get_second_allele() != this_snp.get_second_allele() ) {
					this_snp.second_allele() = '?' ;
				}
			}
		}
		if( *this ) {
			set_number_of_samples( m_number_of_samples ) ;
			set_SNPID( this_snp.get_SNPID() ) ;
			set_RSID( this_snp.get_rsid() ) ;
			set_chromosome( this_snp.get_position().chromosome() ) ;
			set_SNP_position( this_snp.get_position().position() ) ;
			set_allele1( this_snp.get_first_allele() ) ;
			set_allele2( this_snp.get_second_allele() ) ;
		}
	}

	// Read the identifying info of the next SNP in the given source
	// which has the given chromosome and SNP position.
	// If no such SNP is found, throw MissingSNPError.
	// If such a SNP is found but it has different SNPID, RSID, or alleles, throw SNPMismatchError.
	bool SNPDataSourceRack::move_source_to_snp_matching(
		std::size_t source_i,
		SNPIdentifyingData const& reference_snp,
		SNPIdentifyingData* result
	) {
		SNPIdentifyingData this_snp ;
		
		while( m_sources[source_i]->get_next_snp_with_specified_position(
			ignore(),
			set_value( this_snp.SNPID() ),
			set_value( this_snp.rsid() ),
			set_value( this_snp.position().chromosome() ),
			set_value( this_snp.position().position() ),
			set_value( this_snp.first_allele() ),
			set_value( this_snp.second_allele() ),
			reference_snp.get_position().chromosome(),
			reference_snp.get_position().position()
		) ) {
			if( m_comparator.are_equal( this_snp, reference_snp )) {
				*result = this_snp ;
				return true ;
			}

			m_sources[source_i]->ignore_snp_probability_data() ;
		}
		return false ;
	}
		

	namespace impl {
		struct OffsetSampleSetter: public VariantDataReader::PerSampleSetter {
			OffsetSampleSetter(
				VariantDataReader::PerSampleSetter& setter,
				std::size_t offset,
				std::size_t number_of_samples
			):
				m_setter( setter ),
				m_offset( offset ),
				m_number_of_samples( number_of_samples )
			{
				m_setter.set_number_of_samples( m_number_of_samples ) ;
			}

			void set_number_of_samples( std::size_t n ) { /* do nothing. */ }
			void set_sample( std::size_t n ) {
				assert( ( n + m_offset ) < m_number_of_samples ) ;
				m_setter.set_sample( n + m_offset ) ;
			}

			void set_number_of_entries( std::size_t n ) {
				m_setter.set_number_of_entries( n ) ;
			}

			void operator()( MissingValue const value ) { m_setter( value ) ; }
			void operator()( std::string& value ) { m_setter( value ) ; }
			void operator()( Integer const value ) { m_setter( value ) ; }
			void operator()( double const value ) { m_setter( value ) ; }

			void set_offset( std::size_t offset ) { m_offset = offset ; }
			std::size_t get_offset() const { return m_offset ; }

		private:
			VariantDataReader::PerSampleSetter& m_setter ;
			std::size_t m_offset ;
			std::size_t m_number_of_samples ;
		} ;
		
		void add_spec_to_map( std::map< std::string, std::string >* map, std::string const& name, std::string const& type ) {
			(*map)[ name ] = type ;
		}

		struct RackVariantDataReader: public VariantDataReader
		{
			RackVariantDataReader( SNPDataSourceRack& rack )
				: m_rack( rack )
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
				OffsetSampleSetter offset_sample_setter( setter, 0 , m_rack.number_of_samples() ) ;
				for( std::size_t i = 0; i < m_rack.m_sources.size(); ++i ) {
					m_data_readers[i]->get( spec, offset_sample_setter ) ;
					offset_sample_setter.set_offset( offset_sample_setter.get_offset() + m_rack.m_sources[i]->number_of_samples() ) ;
				}
				return *this ;
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
