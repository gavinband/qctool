#ifndef GENFILE_VCF_FORMAT_METADATA_ACCESSOR
#define GENFILE_VCF_FORMAT_METADATA_ACCESSOR

#include <map>
#include <string>
#include <boost/variant.hpp>
#include "genfile/MissingValue.hpp"
#include "genfile/VariantEntry.hpp"

namespace genfile {
	namespace vcf {
		struct FormatSpec {
			FormatSpec( std::map< std::string, std::string > const& values ) ;
			FormatSpec( FormatSpec const& other ) ;
			Entry parse_entry( std::string const& value ) const ;
		private:
			std::string const m_ID ;
			std::string const m_description ;
			Entry const m_number ;
			Type const m_type ;

			InfoSpec& operator=( FormatSpec const& other ) ;
		} ;
		
		typedef std::pair< std::size_t, FormatSpec > LocatedFormatSpec ;

		typedef genfile::VariantEntry Entry ;

		typedef boost::function( std::size_t i, Entry const& entry ) Setter ;
		
		void read_individual_data(
			std::istream& stream,
			std::vector< LocatedFormatSpec > const& format_specs,
			std::vector< Setter > const& setters
		) ;
	}
	
	struct IDDataReader
	{
		typedef boost::function< void ( GenomePosition const& ) > PositionSetter ;
		typedef boost::function< void ( std::vector< std::string > const& ) > IDSetter ;
		typedef boost::function< void ( std::vector< std::string > const& ) > AlleleSetter ;
		typedef boost::function< void ( double ) > QualSetter ;
		void operator()( std::istream&, PositionSetter, IDSetter, AlleleSetter, QualSetter ) const ;
	} ;
	
	struct InfoDataAccessor
	{
		typedef std::multimap< std::string, std::map< std::string, std::string > > Metadata ;
		typedef boost::function< void ( std::size_t i, Entry const& ) > Setter ;
		InfoDataAccessor( Metadata const& metadata ) ;
		
		void read_data(
			std::istream& stream,
			std::string const& field_selection,
			Setter setter
		) ;
		
	private:
		Metadata m_metadata ;
	} ;
	
	struct IndividualDataAccessor
	{
		typedef std::multimap< std::string, std::map< std::string, std::string > > Metadata ;
		typedef boost::function< void ( std::size_t i, Entry const& ) > Setter ;
		IndividualDataAccessor( Metadata const& metadata ) ;
		
		template< typename Setters >
		void read_data(
			std::istream& stream,
			std::string const& field_selection,
			Setters setter
		) const {
			boost::io::ios_flags_saver state_saver( istream ) ;
			stream.exceptions( std::ios::eofbit | std::ios::failbit | std::ios::badbit ) ;
			
			std::size_t const N = boost::tuples::length< ValueTuple >::value ;
			std::string format ;
			stream >> format ;
			if( !stream ) {
				return ;
			}

			std::vector< vcf::LocatedFormatSpec > specs = get_format_specs( std::required_elts, format ) ;
		}
		
		template< typename Setters >
		void read_data(
			std::istream& stream,
			std::vector< vcf::LocatedFormatSpec > const& elts,
			Setters setter
		) const {
			
		}
		
	private:
		std::vector< LocatedFormatSpec > get_format_specs(
			std::string const& format,
			std::vector< std::string > const& required_elts
		) {
			return get_format_specs(
				string_utils::split( format, ":" ),
				string_utils::split( required_elts, ":" )
			) ;
		}
		
		std::vector< LocatedFormatSpec > get_format_specs(
			std::vector< std::string > const& format,
			std::vector< std::string > const& required_elts
		) {
			std::vector< LocatedFormatSpec > result ;
			for( std::size_t i = 0; i < required_elts.size(); ++i ) {
				result.push_back(
					LocatedFormatSpec(
						std::find( format.begin(), format.end(), required_elts.begin[i] ) - format.begin(),
						get_format_spec( required_elts[i] )
					)
				) ;
			}
			return result ;
		}
		
		FormatSpec get_format_spec( std::string const& id ) const {
			std::map< std::string, vcf::FormatSpec >::const_iterator where = m_format_specs.find( id ) ;
			if( id == m_format_specs.end() ) {
				throw BadArgumentError( "genfile::vcf::IndividualDataAccessor::get_format_spec()", "id = \"" + id "\"." ) ;
			}
			return where.second ;
		}
		
	private:
		Metadata m_metadata ;
		std::map< std::string, vcf::FormatSpec > m_format_specs ;
	} ;
}

#endif
