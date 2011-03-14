#include <string>
#include <map>
#include <memory>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include "genfile/SNPDataSource.hpp"
#include "genfile/snp_data_utils.hpp"
#include "genfile/string_utils.hpp"
#include "genfile/Error.hpp"
#include "genfile/VCFFormatSNPDataSource.hpp"
#include "genfile/vcf/Types.hpp"
#include "genfile/vcf/InfoReader.hpp"
#include "genfile/vcf/CallReader.hpp"
#include "genfile/vcf/get_set.hpp"

using boost::tuples::tie ;

namespace genfile {
	VCFFormatSNPDataSource::VCFFormatSNPDataSource(
		std::auto_ptr< std::istream > stream_ptr,
		std::string const& genotype_probability_field
	):
		m_spec( "(unnamed stream)" ),
		m_stream_ptr( stream_ptr ),
		m_metadata( vcf::MetadataParser( m_spec, *m_stream_ptr ).get_metadata() ),
		m_info_types( vcf::get_entry_types( m_metadata, "INFO" )),
		m_format_types( vcf::get_entry_types( m_metadata, "FORMAT" )),
		m_genotype_probability_field( genotype_probability_field ),
		m_column_names( read_column_names( *m_stream_ptr )),
		m_start_of_data( m_stream_ptr->tellg() ),
		m_number_of_lines( count_lines( *m_stream_ptr ))
	{
		setup() ;
	}

	namespace impl {
		std::string get_file_part_of_spec( std::string const& spec ) {
			std::vector< std::string > elts = string_utils::split( spec, ":" ) ;
			if( !elts.size() == 2 ) {
				throw BadArgumentError( "genfile::impl::get_file_part_of_spec()", "spec = \"" + spec +"\"" ) ;
			}
			return elts[0] ;
		}

		std::string get_field_part_of_spec( std::string const& spec ) {
			std::vector< std::string > elts = string_utils::split( spec, ":" ) ;
			if( !elts.size() == 2 ) {
				throw BadArgumentError( "genfile::impl::get_file_part_of_spec()", "spec = \"" + spec +"\"" ) ;
			}
			return elts[1] ;
		}
	}
	
	VCFFormatSNPDataSource::VCFFormatSNPDataSource(
		std::string const& filename,
		std::string const& genotype_probability_field
	):
		m_spec( "file://" + filename ),
		m_stream_ptr( open_text_file_for_input( filename, get_compression_type_indicated_by_filename( filename ))),
		m_metadata( vcf::MetadataParser( m_spec, *m_stream_ptr ).get_metadata() ),
		m_info_types( vcf::get_entry_types( m_metadata, "INFO" )),
		m_format_types( vcf::get_entry_types( m_metadata, "FORMAT" )),
		m_genotype_probability_field( genotype_probability_field ),
		m_column_names( read_column_names( *m_stream_ptr )),
		m_start_of_data( m_stream_ptr->tellg() ),
		m_number_of_lines( count_lines( *m_stream_ptr ))
	{
		setup() ;
	}
	
	void VCFFormatSNPDataSource::setup() {
		EntryTypeMap::const_iterator where = m_format_types.find( m_genotype_probability_field ) ;
		if(
			where == m_format_types.end()
			|| (
				dynamic_cast< vcf::GenotypeCallVCFEntryType const* >( &(*where->second) ) == 0
			&& dynamic_cast< vcf::OnePerGenotypeVCFEntryType const* >( &(*where->second) ) == 0
			)
		) {
			throw BadArgumentError(
				"genfile::VCFFormatSNPDataSource::setup()",
				"m_genotype_probability_field = \"" + m_genotype_probability_field + "\""
			) ;
		}
		reset_stream() ;
		m_stream_ptr->exceptions( std::ios::eofbit | std::ios::failbit | std::ios::badbit ) ;
	}

	void VCFFormatSNPDataSource::reset_stream() const {
		m_stream_ptr->clear() ;
		m_stream_ptr->seekg( m_start_of_data ) ;
		assert( *m_stream_ptr ) ;
	}

	std::vector< std::string > VCFFormatSNPDataSource::read_column_names( std::istream& stream ) const {
		std::string line ;
		if( !std::getline( stream, line ) ) {
			throw MalformedInputError( m_spec, m_metadata.size() ) ;
		}
		std::vector< std::string > elts = string_utils::split( line, "\t" ) ;
		if( elts.size() < 8 ) {
			throw MalformedInputError( m_spec, m_metadata.size() ) ;
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
			throw MalformedInputError( m_spec, m_metadata.size() ) ;
		}
		return elts ;
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
		return m_column_names.size() - 8 ;
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
	
	void VCFFormatSNPDataSource::set_genotype_probability_field( std::string const& value ) {
		m_genotype_probability_field = value ;
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
			read_element( CHROM, '\t', entry_count++ ) ;
		}
		catch( std::ios_base::failure const& ) {
			// end of data, this is not an error.
			return ;
		}

		try {
			read_element( POS, '\t', entry_count++ ) ;
			read_element( ID, '\t', entry_count++ ) ;
			read_element( REF, '\t', entry_count++ ) ;
			read_element( ALT, '\t', entry_count++ ) ;
			read_element( QUAL, '\t', entry_count++ ) ;
			read_element( FILTER, '\t', entry_count++ ) ;
			read_element( INFO, '\t', entry_count++ ) ;
		}
		catch( std::ios_base::failure const& ) {
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata.size() + 1, entry_count ) ;
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

		if( m_variant_alleles.size() >= 1 && m_variant_alleles[0].size() == 1 && m_variant_alleles[0] != "." ) {
			set_allele1( m_variant_alleles[0][0] ) ;
		}
		else {
			set_allele1( '?' ) ;
		}

		if( m_variant_alleles.size() == 2 && m_variant_alleles[1].size() == 1 && m_variant_alleles[1] != "." ) {
			set_allele2( m_variant_alleles[1][0] ) ;
		}
		else {
			set_allele2( '?' ) ;
		}
	}

	void VCFFormatSNPDataSource::read_element( std::string& elt, char delim, std::size_t column ) const {
		std::getline( *m_stream_ptr, elt, delim ) ;
		// ensure element contains no whitespace.
		if( elt.find_first_of( " \t\n\r" ) != std::string::npos ) {
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata.size() + 1, column ) ;
		}
	}
	
	void VCFFormatSNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		std::string FORMAT ;
		std::string data ;
		std::size_t count = 7 ;
		try {
			read_element( FORMAT, '\t', count++ ) ;
			std::getline( *m_stream_ptr, data ) ;
			++count ;
		}
		catch( std::ios_base::failure const& ) {
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata.size() + 1, count ) ;
		}

		try {
			vcf::CallReader( m_variant_alleles.size(), FORMAT, data, m_format_types )
				( m_genotype_probability_field, vcf::make_genotype_probability_setter( set_genotype_probabilities ) ) ;
		}
		catch( BadArgumentError const& ) {
			// problem with FORMAT
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata.size() + 1, 8 ) ;
		}
		catch( MalformedInputError const& e ) {
			// problem with entry
			if( e.has_column() ) {
				throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata.size() + 1, e.column() + 8 ) ;
			}
			else {
				throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata.size() + 1 ) ;
			}
		}
	}

	void VCFFormatSNPDataSource::ignore_snp_probability_data_impl() {
		std::string FORMAT ;
		std::string data ;
		std::size_t count = 7 ;
		try {
			read_element( FORMAT, '\t', count++ ) ;
			std::getline( *m_stream_ptr, data ) ;
			++count ;
		}
		catch( std::ios_base::failure const& ) {
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata.size() + 1, count ) ;
		}
		// We ignore the data, but parse the FORMAT spec anyway.
		try {
			vcf::CallReader( m_variant_alleles.size(), FORMAT, data, m_format_types ) ;
		}
		catch( BadArgumentError const& ) {
			// problem with FORMAT
			throw MalformedInputError( get_source_spec(), number_of_snps_read() + m_metadata.size() + 1, 8 ) ;
		}
	}

	void VCFFormatSNPDataSource::reset_to_start_impl() {
		// seek back to the start and skip over metadata and column names.
		reset_stream() ;
	}	
}
