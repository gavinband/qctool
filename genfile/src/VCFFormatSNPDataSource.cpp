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
		std::auto_ptr< std::istream > stream_ptr
	):
		m_spec( "(unnamed stream)" ),
		m_compression_type( "no_compression" ),
		m_stream_ptr( stream_ptr ),
		m_metadata_parser( new vcf::StrictMetadataParser( m_spec, *m_stream_ptr ) ),
		m_metadata( m_metadata_parser->get_metadata() ),
		m_info_types( vcf::get_entry_types( m_metadata, "INFO" )),
		m_format_types( vcf::get_entry_types( m_metadata, "FORMAT" )),
		m_field_mapping( get_field_mapping( m_format_types )),
		m_column_names( read_column_names( *m_stream_ptr )),
		m_number_of_samples( m_column_names.size() - 9 )
	{
		setup() ;
		m_number_of_lines = determine_number_of_lines( *m_stream_ptr, m_metadata ) ;
	}

	VCFFormatSNPDataSource::VCFFormatSNPDataSource(
		std::string const& filename
	):
		m_spec( filename ),
		m_compression_type( get_compression_type_indicated_by_filename( filename )),
		m_stream_ptr( open_text_file_for_input( filename, m_compression_type ) ),
		m_metadata_parser( new vcf::StrictMetadataParser( m_spec, *m_stream_ptr ) ),
		m_metadata( m_metadata_parser->get_metadata() ),
		m_info_types( vcf::get_entry_types( m_metadata, "INFO" )),
		m_format_types( vcf::get_entry_types( m_metadata, "FORMAT" )),
		m_field_mapping( get_field_mapping( m_format_types )),
		m_column_names( read_column_names( *m_stream_ptr )),
		m_number_of_samples( m_column_names.size() - 9 )
	{
		setup() ;
		m_number_of_lines = determine_number_of_lines( *m_stream_ptr, m_metadata ) ;
	}

	VCFFormatSNPDataSource::VCFFormatSNPDataSource(
		std::string const& filename,
		Metadata const& metadata
	):
		m_spec( filename ),
		m_compression_type( get_compression_type_indicated_by_filename( filename )),
		m_stream_ptr( open_text_file_for_input( filename, m_compression_type ) ),
		m_metadata_parser( new vcf::TrivialMetadataParser( m_spec, *m_stream_ptr ) ),
		m_metadata( metadata ),
		m_info_types( vcf::get_entry_types( m_metadata, "INFO" )),
		m_format_types( vcf::get_entry_types( m_metadata, "FORMAT" )),
		m_field_mapping( get_field_mapping( m_format_types )),
		m_column_names( read_column_names( *m_stream_ptr )),
		m_number_of_samples( m_column_names.size() - 9 )
	{
		setup() ;
		m_number_of_lines = determine_number_of_lines( *m_stream_ptr, m_metadata ) ;
	}
	
	void VCFFormatSNPDataSource::setup() {
		std::string index_filename = m_spec ;
		if( index_filename.size() > 3 && index_filename.substr( index_filename.size() - 3, 3 ) == ".gz" ) {
			index_filename = index_filename.substr( 0, index_filename.size() - 3 ) ;
		}
		index_filename += ".index" ;
		if( boost::filesystem::exists( index_filename ) ) {
			m_index_stream_ptr = open_text_file_for_input( index_filename ) ;
		} 
		// check_genotype_probability_field( m_field_mapping ) ;
		reset_stream() ;
	}

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
		for( std::size_t i = 0; i < ( m_metadata_parser->get_number_of_lines() + 1 ); ++i ) {
			std::string line ;
			std::getline( *m_stream_ptr, line ) ;
		}

		assert( *m_stream_ptr ) ;
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

	std::size_t VCFFormatSNPDataSource::determine_number_of_lines(
		std::istream& vcf_file_stream,
		vcf::MetadataParser::Metadata const& metadata
	) const {
		typedef vcf::MetadataParser::Metadata::const_iterator MetadataIterator ;
		std::pair< MetadataIterator, MetadataIterator > range = metadata.equal_range( "number-of-variants" ) ;
		std::size_t result ;
		if( range.first != range.second ) {
			std::map< std::string, std::string >::const_iterator where = range.first->second.find( "" ) ;
			if( where == range.first->second.end() ) {
				throw MalformedInputError( m_spec, std::distance( metadata.begin(), range.first )) ;
			}
			else {
				result = string_utils::to_repr< std::size_t >( where->second ) ;
			}
			if( (++range.first) != range.second ) {
				throw MalformedInputError( m_spec, std::distance( metadata.begin(), range.first )) ;
			}
		}
		else {
			if( m_index_stream_ptr.get() ) {
				result = count_lines( *m_index_stream_ptr ) - 1 ; // ignore header line.
			}
			else {
				result = count_lines( vcf_file_stream ) ;
			}
		}
		return result ;
	}

	std::size_t VCFFormatSNPDataSource::count_lines( std::istream& str ) const {
		std::size_t count = 0 ;
		std::vector< char > buffer( 10000000 ) ;
		std::size_t last_read_size = 0 ;
		do {
			str.read( &(buffer[0]), 10000000 ) ;
			if( str.gcount() > 0 ) {
				last_read_size = str.gcount() ;
			}
			count += std::count( buffer.begin(), buffer.begin() + str.gcount(), '\n' ) ;
			// A vcf file can't contain a blank line.
			// Because popular editors (vim, nano, ..., but not emacs) typically add a trailing newline,
			// we might get in the situation where the file has two trailing newlines thus messing
			// with our count.
			// Therefore we check here for the special case where what we've read ends in two newlines.
			if( (str.gcount() > 1) && (buffer[ str.gcount() - 1] == '\n') && (buffer[ str.gcount() - 2] == '\n') ) {
				throw FileHasTwoTrailingNewlinesError( get_source_spec(), count ) ;
			}
		}
		while( str ) ;
		assert( last_read_size < buffer.size() ) ;
		if( last_read_size == 0 ) {
			throw MalformedInputError( get_source_spec(), 0 ) ;
		}
		else if( buffer[ last_read_size - 1 ] != '\n' ) {
			throw MissingTrailingNewlineError( get_source_spec(), count ) ;
		}
		return count ;
	}
	
	VCFFormatSNPDataSource::operator bool() const {
		return *m_stream_ptr ;
	}

	unsigned int VCFFormatSNPDataSource::number_of_samples() const {
		return m_column_names.size() - 9 ;
	}

	unsigned int VCFFormatSNPDataSource::total_number_of_snps() const {
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
		std::cerr << "Afterwards:\n" ;
		for( FieldMapping::left_iterator where = m_field_mapping.left.begin(); where != m_field_mapping.left.end(); ++where ) {
			std::cerr << "  " << where->first << ":" << where->second << "\n" ;
		}
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
	
	void VCFFormatSNPDataSource::list_snps( SNPSetter setter, ProgressCallback progress_callback ) {
		std::istream* source = m_index_stream_ptr.get() ;
		if( !source ) {
			SNPDataSource::list_snps( setter, progress_callback ) ;
			return ;
		}
		else {
			using string_utils::slice ;

			std::cerr << "Reading from index file.\n" ;

			source->clear() ;
			source->seekg(0) ;
			std::string line ;

			std::getline( *source, line ) ;
			if( line != "SNPID RSID chromosome position allele1 allele2" ) {
				throw MalformedInputError( m_spec + ".index", 0 ) ;
			}

			for( std::size_t count = 1; std::getline( *source, line ); ++count ) {
				std::vector< slice > elts = slice( line ).split( " " ) ;
				if( elts.size() != 6 ) {
					throw MalformedInputError( m_spec + ".index", count, 6 ) ;
				}
				if( elts[4].size() != 1 ) {
					throw MalformedInputError( m_spec + ".index", count, 4 ) ;
				}
				if( elts[5].size() != 1 ) {
					throw MalformedInputError( m_spec + ".index", count, 5 ) ;
				}
				setter(
					SNPIdentifyingData(
						elts[0],
						elts[1],
						GenomePosition(
							Chromosome( elts[2] ),
							string_utils::to_repr< Position >( elts[3] )
						),
						elts[4][0],
						elts[5][0]
					)
			 	) ;
				if( progress_callback ) {
					progress_callback( number_of_snps_read() + 1, total_number_of_snps() ) ;
				}
			}
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
		std::string CHROM ;
		std::string POS ;
		std::string ID ;
		std::string REF ;
		std::string ALT ;
		std::string QUAL ;
		std::string FILTER ;
		std::string INFO ;

		std::size_t entry_count = 0 ;
		try {
			impl::read_element( *stream_ptr, CHROM, '\t', entry_count++ ) ;
			if( number_of_snps_read() == m_number_of_lines ) {
				throw MalformedInputError( m_spec, number_of_snps_read() ) ;
			}
		}
		catch( std::ios_base::failure const& ) {
			// end of data, this is only an error if number of lines did not match.
			if( number_of_snps_read() != m_number_of_lines ) {
				throw MalformedInputError( m_spec, number_of_snps_read() ) ;
			}
			return ;
		}
		catch( MalformedInputError const& e ) {
			throw MalformedInputError( m_spec, number_of_snps_read(), e.column() ) ;
		}

		try {
			impl::read_element( *stream_ptr,POS, '\t', entry_count++ ) ;
			impl::read_element( *stream_ptr,ID, '\t', entry_count++ ) ;
			impl::read_element( *stream_ptr,REF, '\t', entry_count++ ) ;
			impl::read_element( *stream_ptr,ALT, '\t', entry_count++ ) ;
			impl::read_element( *stream_ptr,QUAL, '\t', entry_count++ ) ;
			impl::read_element( *stream_ptr,FILTER, '\t', entry_count++ ) ;
			impl::read_element( *stream_ptr,INFO, '\t', entry_count++ ) ;
		}
		catch( std::ios_base::failure const& ) {
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata_parser->get_number_of_lines() + 1, entry_count ) ;
		}
		catch( MalformedInputError const& e ) {
			throw MalformedInputError( m_spec, number_of_snps_read(), e.column() ) ;
		}
		
		// If we get here reading was successful.
		set_number_of_samples( number_of_samples() ) ;

		std::vector< std::string > ids = string_utils::split( ID, "," ) ;
		if( ids.size() == 0 ) {
			set_RSID( "?" ) ;
			set_SNPID( "?" ) ;
		}
		else if( ids.size() == 1 ) {
			set_RSID( ids[0] ) ;
			set_SNPID( "?" ) ;
		}
		else {
			set_RSID( ids[0] ) ;
			set_SNPID( ids[1] ) ;
		}
		
		set_chromosome( Chromosome( CHROM ) ) ;
		set_SNP_position( string_utils::to_repr< Position >( POS ) ) ;

		m_variant_alleles = string_utils::split( ALT, "," ) ;
		m_variant_alleles.insert( m_variant_alleles.begin(), REF ) ;

	 	// We ignore the INFO, but may as well verify that it's sane
		// insofar as we can easily do so.
		vcf::InfoReader( m_variant_alleles.size(), INFO, m_info_types ) ;

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
							throw MalformedInputError( m_source.get_source_spec(), m_source.number_of_snps_read() -1 + m_source.get_index_of_first_data_line() + m_source.get_index_of_first_data_column() ) ;
						}
					}
				}
				return *this ;
			}
			
			bool supports( std::string const& spec ) const {
				return m_field_mapping.left.find( spec ) != m_field_mapping.left.end() ;
			}

			void get_supported_specs( SpecSetter setter ) const {
				// We support everything on the left of the field mapping.
				FieldMapping::left_map::const_iterator 
					begin_i = m_field_mapping.left.begin(),
					end_i = m_field_mapping.left.end() ;
				for( ; begin_i != end_i; ++begin_i ) {
					setter( begin_i->first, m_format_types.find( begin_i->second )->second->get_type().to_string() ) ;
				}
			}

		private:
			VCFFormatSNPDataSource const& m_source ;
			boost::ptr_map< std::string, vcf::VCFEntryType > const& m_format_types ;
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
				throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata_parser->get_number_of_lines() + 1, 8 ) ;
			}
		}
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
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata_parser->get_number_of_lines() + 1, count ) ;
		}
		// We ignore the data, but parse the FORMAT spec anyway.
		try {
			vcf::CallReader( m_number_of_samples, m_variant_alleles.size(), FORMAT, data, m_format_types ) ;
		}
		catch( BadArgumentError const& ) {
			// problem with FORMAT
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata_parser->get_number_of_lines() + 1, 8 ) ;
		}
	}

	void VCFFormatSNPDataSource::reset_to_start_impl() {
		// seek back to the start and skip over metadata and column names.
		reset_stream() ;
	}	
}
