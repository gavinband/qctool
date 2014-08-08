
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <string>
#include <map>
#include <set>
#include <memory>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <boost/bimap.hpp>
#include <boost/filesystem/operations.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/snp_data_utils.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/VCFFormatSNPDataSource.hpp"
#include "genfile/vcf/StrictMetadataParser.hpp"
#include "genfile/vcf/TrivialMetadataParser.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/vcf/InfoReader.hpp"
#include "genfile/vcf/CallReader.hpp"
#include "genfile/vcf/get_set.hpp"

using boost::tuples::tie ;

namespace genfile {
	namespace {
		typedef boost::bimap< std::string, std::string > FieldMapping ;
		FieldMapping get_field_mapping(
			boost::ptr_map< std::string, vcf::VCFEntryType > const& format_types
		) {
			boost::ptr_map< std::string, vcf::VCFEntryType >::const_iterator
				begin = format_types.begin(),
				end = format_types.end() ;
			FieldMapping result ;
			for( ; begin != end; ++begin ) {
				result.insert( FieldMapping::value_type( begin->first, begin->first )) ;
			}
			return result ;
		}
	}

	VCFFormatSNPDataSource::VCFFormatSNPDataSource(
		std::auto_ptr< std::istream > stream_ptr,
		boost::optional< Metadata > metadata
	):
		m_spec( "(unnamed stream)" ),
		m_compression_type( "no_compression" ),
		m_stream_ptr( stream_ptr ),
		m_metadata_parser( new vcf::StrictMetadataParser( m_spec, *m_stream_ptr ) ),
		m_metadata( ( metadata ? *metadata : m_metadata_parser->get_metadata() ) ),
		m_info_types( vcf::get_entry_types( m_metadata, "INFO" )),
		m_format_types( vcf::get_entry_types( m_metadata, "FORMAT" )),
		m_field_mapping( get_field_mapping( m_format_types )),
		m_column_names( read_column_names( *m_stream_ptr )),
		m_number_of_samples( m_column_names.size() - 9 ),
		m_have_id_data( false )
	{
		setup() ;
		reset_stream() ;
	}

	VCFFormatSNPDataSource::VCFFormatSNPDataSource(
		std::string const& filename,
		boost::optional< Metadata > metadata
	):
		m_spec( filename ),
		m_compression_type( get_compression_type_indicated_by_filename( filename )),
		m_stream_ptr( open_text_file_for_input( filename, m_compression_type ) ),
		m_metadata_parser( new vcf::StrictMetadataParser( m_spec, *m_stream_ptr ) ),
		m_metadata( ( metadata ? *metadata : m_metadata_parser->get_metadata() ) ),
		m_info_types( vcf::get_entry_types( m_metadata, "INFO" )),
		m_format_types( vcf::get_entry_types( m_metadata, "FORMAT" )),
		m_field_mapping( get_field_mapping( m_format_types )),
		m_column_names( read_column_names( *m_stream_ptr )),
		m_number_of_samples( m_column_names.size() - 9 ),
		m_have_id_data( false )
	{
		setup() ;
		reset_stream() ;
	}

	void VCFFormatSNPDataSource::setup() {}

	void VCFFormatSNPDataSource::check_genotype_probability_field( std::string const& field ) const {
		EntryTypeMap::const_iterator where = m_format_types.find( field ) ;
		if(
			where == m_format_types.end()
		) {
			throw KeyNotFoundError( "FORMAT definitions of " + get_source_spec(), field ) ;
		}
		else if( dynamic_cast< vcf::GenotypeCallVCFEntryType const* >( &(*where->second) ) != 0  ) {
			// fine, do nothing.
		}
		else {
			vcf::ListVCFEntryType const* type = dynamic_cast< vcf::ListVCFEntryType const* >( &(*where->second) ) ;
			if(
				type == 0
				|| 
				// We allow either fields with 3 or 4 entries.
				// (In the latter case the 4th entry is treated as a NULL genotype call).
				// Here we just check that the entry type's range does not exclude this number of entries.
				type->get_value_count_range( 2, 2 ).first > 4
				||
				type->get_value_count_range( 2, 2 ).second < 3
			) {
				throw BadArgumentError(
					"genfile::VCFFormatSNPDataSource::check_genotype_probability_field()",
					"field = \"" + field + "\""
				) ;
			}
		}
	}

	void VCFFormatSNPDataSource::reset_stream() {
		m_stream_ptr->clear() ;
		if( m_compression_type == "no_compression" ) {
			m_stream_ptr->seekg( 0 ) ;
		}
		else if( m_spec != "(unnamed stream)" ) {
			m_stream_ptr = open_text_file_for_input( m_spec, m_compression_type ) ;
		}
		else {
			throw OperationUnsupportedError( "void VCFFormatSNPDataSource::reset_stream()", "open stream", m_spec ) ;
		}

		m_stream_ptr->exceptions( std::ios::eofbit | std::ios::failbit | std::ios::badbit ) ;

		// Find our way back to the start of data.
		for( std::size_t i = 0; i < ( get_index_of_first_data_line() ); ++i ) {
			std::string line ;
			std::getline( *m_stream_ptr, line ) ;
		}

		assert( *m_stream_ptr ) ;
		m_have_id_data = false ;
	}

	std::vector< std::string > VCFFormatSNPDataSource::read_column_names( std::istream& stream ) const {
		std::string line ;
		if( !std::getline( stream, line ) ) {
			throw MalformedInputError( m_spec, m_metadata_parser->get_number_of_lines() ) ;
		}
		std::vector< std::string > elts = string_utils::split( line, "\t" ) ;
		if( elts.size() < 8 ) {
			throw MalformedInputError( m_spec, m_metadata_parser->get_number_of_lines() ) ;
		}
		else if(
			elts[0] != "#CHROM"
			|| elts[1] != "POS"
			|| elts[2] != "ID"
			|| elts[3] != "REF"
			|| elts[4] != "ALT"
			|| elts[5] != "QUAL"
			|| elts[6] != "FILTER"
			|| elts[7] != "INFO"
			|| elts[8] != "FORMAT"
		) {
			throw MalformedInputError( m_spec, m_metadata_parser->get_number_of_lines() ) ;
		}
		return elts ;
	}

	VCFFormatSNPDataSource::operator bool() const {
		return *m_stream_ptr ;
	}

	unsigned int VCFFormatSNPDataSource::number_of_samples() const {
		return m_column_names.size() - 9 ;
	}

	SNPDataSource::OptionalSnpCount VCFFormatSNPDataSource::total_number_of_snps() const {
		return m_number_of_lines ;
	}

	std::string VCFFormatSNPDataSource::get_source_spec() const {
		return m_spec ;
	}

	std::string VCFFormatSNPDataSource::get_summary( std::string const& prefix, std::size_t column_width ) const {
		return prefix + m_spec ;
	}
	
	void VCFFormatSNPDataSource::set_field_mapping( std::string const& key, std::string const& value ) {
		FieldMapping::right_iterator where = m_field_mapping.right.find( value ) ;
		if( where == m_field_mapping.right.end() ) {
			throw BadArgumentError( "genfile::VCFFormatSNPDataSource::set_field_mapping()", "value = \"" + value + "\"" ) ;
		}
		m_field_mapping.right.replace_data( where, key ) ;
	}
	
	
	namespace impl {
		void read_element( std::istream& stream, std::string& elt, char delim, std::size_t column ) {
			std::getline( stream, elt, delim ) ;
			// ensure element contains no whitespace.
			if( elt.find_first_of( " \t\n\r" ) != std::string::npos ) {
				throw MalformedInputError( "(unnamed stream)", 0, column ) ;
			}
		}

		char read_format_and_get_trailing_char( std::istream& stream, std::string& format, std::size_t column ) {
			char after_format ;
			for( after_format = stream.get(); after_format != '\n' && after_format != '\t'; after_format = stream.get() ) {
				format.append( 1, after_format ) ;
			}
			if( format.find_first_of( " \t\n\r" ) != std::string::npos ) {
				throw MalformedInputError( "(unnamed stream)", 0, column ) ;
			}
			return after_format ;
		}
	}
	
	void VCFFormatSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		return get_snp_identifying_data_impl(
			m_stream_ptr.get(),
			set_number_of_samples,
			set_SNPID,
			set_RSID,
			set_chromosome,
			set_SNP_position,
			set_allele1,
			set_allele2
		) ;
	}

	void VCFFormatSNPDataSource::get_snp_identifying_data_impl( 
		std::istream* stream_ptr,
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		if( !m_have_id_data ) {
			std::size_t entry_count = 0 ;
			try {
				impl::read_element( *stream_ptr, m_CHROM, '\t', entry_count++ ) ;
				if( m_number_of_lines && number_of_snps_read() == *m_number_of_lines ) {
					throw MalformedInputError( m_spec, number_of_snps_read() + get_index_of_first_data_line() ) ;
				}
			}
			catch( std::ios_base::failure const& ) {
				// end of data, this is only an error if number of lines did not match.
				if( m_number_of_lines && number_of_snps_read() != *m_number_of_lines ) {
					throw MalformedInputError( m_spec, number_of_snps_read() + get_index_of_first_data_line() ) ;
				}
				return ;
			}
			catch( MalformedInputError const& e ) {
				throw MalformedInputError( m_spec, number_of_snps_read() + get_index_of_first_data_line(), e.column() ) ;
			}

			try {
				impl::read_element( *stream_ptr, m_POS, '\t', entry_count++ ) ;
				impl::read_element( *stream_ptr, m_ID, '\t', entry_count++ ) ;
				impl::read_element( *stream_ptr, m_REF, '\t', entry_count++ ) ;
				impl::read_element( *stream_ptr, m_ALT, '\t', entry_count++ ) ;
				impl::read_element( *stream_ptr, m_QUAL, '\t', entry_count++ ) ;
				impl::read_element( *stream_ptr, m_FILTER, '\t', entry_count++ ) ;
				impl::read_element( *stream_ptr, m_INFO, '\t', entry_count++ ) ;
			}
			catch( std::ios_base::failure const& ) {
				throw MalformedInputError( get_source_spec(), number_of_snps_read() + get_index_of_first_data_line(), entry_count ) ;
			}
			catch( MalformedInputError const& e ) {
				throw MalformedInputError( m_spec, number_of_snps_read() + get_index_of_first_data_line(), e.column() ) ;
			}
			
			m_have_id_data = true ;
		}
		// If we get here reading was successful.
		set_number_of_samples( number_of_samples() ) ;

		std::vector< std::string > ids = string_utils::split( m_ID, "," ) ;
		if( ids.size() == 0 ) {
			set_RSID( "?" ) ;
			set_SNPID( "?" ) ;
		}
		else if( ids.size() == 1 ) {
			set_RSID( ids[0] ) ;
			set_SNPID( ids[0] ) ;
		}
		else {
			set_RSID( ids[0] ) ;
			set_SNPID( ids[1] ) ;
		}
		
		set_chromosome( Chromosome( m_CHROM ) ) ;
		set_SNP_position( string_utils::to_repr< Position >( m_POS ) ) ;

		m_variant_alleles = string_utils::split( m_ALT, "," ) ;
		m_variant_alleles.insert( m_variant_alleles.begin(), m_REF ) ;

	 	// We ignore the INFO, but may as well verify that it's sane
		// insofar as we can easily do so.
		vcf::InfoReader( m_variant_alleles.size(), m_INFO, m_info_types ) ;

		if( m_variant_alleles.size() >= 1 && m_variant_alleles[0] != "." ) {
			set_allele1( m_variant_alleles[0] ) ;
		}
		else {
			set_allele1( "?" ) ;
		}

		if( m_variant_alleles.size() == 2 && m_variant_alleles[1] != "." ) {
			set_allele2( m_variant_alleles[1] ) ;
		}
		else {
			set_allele2( "?" ) ;
		}
	}

	namespace impl {
		// Wrap a vcf::CallReader with its VCFFormatSNPDataSource in such a way that sensible error messages are returned.
		struct VCFFormatDataReader: public VariantDataReader {
			VCFFormatDataReader(
				VCFFormatSNPDataSource const& source,
				std::vector< std::string > const& variant_alleles,
				std::string const& FORMAT,
				boost::ptr_map< std::string, vcf::VCFEntryType > const& format_types,
				VCFFormatSNPDataSource::FieldMapping field_mapping
			):
			 	m_source( source ),
				m_format_elts( genfile::string_utils::split( FORMAT, ":" )),
				m_format_types( format_types ),
				m_field_mapping( field_mapping )
			{
				if( m_source.number_of_samples() > 0 ) {
					std::string data ;
					std::getline( *(m_source.m_stream_ptr), data ) ;
					m_data_reader.reset( new vcf::CallReader( m_source.number_of_samples(), variant_alleles.size(), FORMAT, data, format_types ) ) ;
				}
			}

		public:
			VCFFormatDataReader& get( std::string const& spec, VariantDataReader::PerSampleSetter& setter ) {
				if( m_source.number_of_samples() > 0 ) {
					try {
						FieldMapping::left_const_iterator where = m_field_mapping.left.find( spec ) ;
						if( where != m_field_mapping.left.end() ) {
							m_data_reader->get( where->second, setter ) ;
						}
						else {
							throw OperationUnsupportedError( "genfile::impl::VCFFormatDataReader::get()", "get \"" + spec + "\"", m_source.get_source_spec() ) ;
						}
					}
					catch( MalformedInputError const& e ) {
						// problem with entry.
						if( e.has_column() ) {
							// error column is the individual index (starting from 0), we add 9 to get the column number.
							throw MalformedInputError( m_source.get_source_spec(), m_source.number_of_snps_read() - 1 + m_source.get_index_of_first_data_line(), e.column() + m_source.get_index_of_first_data_column() ) ;
						}
						else {
							throw MalformedInputError( m_source.get_source_spec(), m_source.number_of_snps_read() -1 + m_source.get_index_of_first_data_line() ) ;
						}
					}
				}
				return *this ;
			}
			
			std::size_t get_number_of_samples() const { return m_source.number_of_samples() ; }
			
			bool supports( std::string const& spec ) const {
				return m_field_mapping.left.find( spec ) != m_field_mapping.left.end() ;
			}

			void get_supported_specs( SpecSetter setter ) const {
				for( std::size_t i = 0; i < m_format_elts.size(); ++i ) {
					setter( m_format_elts[i], m_format_types.find( m_format_elts[i] )->second->get_type().to_string() ) ;
				}
			}

		private:
			VCFFormatSNPDataSource const& m_source ;
			boost::ptr_map< std::string, vcf::VCFEntryType > const& m_format_types ;
			std::vector< std::string > const m_format_elts ;
			typedef VCFFormatSNPDataSource::FieldMapping FieldMapping ;
			FieldMapping m_field_mapping ;
			vcf::CallReader::UniquePtr m_data_reader ;
		} ;
	}
	
	std::string VCFFormatSNPDataSource::read_format() {
		std::string FORMAT ;
		std::size_t count = 8 ;
		char after_format ;
		try {
			after_format = impl::read_format_and_get_trailing_char( *m_stream_ptr,FORMAT, count++ ) ;
		}
		catch( std::ios_base::failure const& ) {
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + get_index_of_first_data_line(), count ) ;
		}

		if(( m_number_of_samples == 0 && after_format != '\n' ) || ( m_number_of_samples > 0 && after_format != '\t' )) {
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + get_index_of_first_data_line(), count ) ;
		}

		return FORMAT ;
	}

	VariantDataReader::UniquePtr VCFFormatSNPDataSource::read_variant_data_impl() {
		std::string FORMAT = read_format() ;
		VariantDataReader::UniquePtr result ;
		if( m_number_of_samples > 0 ) {
			try {
				result.reset( new impl::VCFFormatDataReader( *this, m_variant_alleles, FORMAT, m_format_types, m_field_mapping )) ;
			}
			catch( BadArgumentError const& ) {
				// problem with FORMAT
				throw MalformedInputError( get_source_spec(), number_of_snps_read() + get_index_of_first_data_line(), 8 ) ;
			}
		}
		m_have_id_data = false ;
		return result ;
	}
	
	void VCFFormatSNPDataSource::ignore_snp_probability_data_impl() {
		std::string FORMAT ;
		std::string data ;
		std::size_t count = 7 ;
		try {
			impl::read_element( *m_stream_ptr, FORMAT, '\t', count++ ) ;
			std::getline( *m_stream_ptr, data ) ;
			++count ;
		}
		catch( std::ios_base::failure const& ) {
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + get_index_of_first_data_line(), count ) ;
		}
		// We ignore the data, but parse the FORMAT spec anyway.
		try {
			vcf::CallReader( m_number_of_samples, m_variant_alleles.size(), FORMAT, data, m_format_types ) ;
		}
		catch( BadArgumentError const& ) {
			// problem with FORMAT
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + get_index_of_first_data_line(), 8 ) ;
		}
		m_have_id_data = false ;
	}

	void VCFFormatSNPDataSource::reset_to_start_impl() {
		// seek back to the start and skip over metadata and column names.
		reset_stream() ;
	}	
}
