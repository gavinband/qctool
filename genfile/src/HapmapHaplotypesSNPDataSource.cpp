
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/HapmapHaplotypesSNPDataSource.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/FileUtils.hpp"
#include "genfile/Error.hpp"

namespace genfile {
	HapmapHaplotypesSNPDataSource::HapmapHaplotypesSNPDataSource( std::auto_ptr< std::istream > stream, Chromosome chromosome ):
		m_filename( "(unnamed stream)" ),
		m_compression_type( "no_compression" ),
		m_number_of_samples( 0 ),
		m_total_number_of_snps( 0 ),
		m_chromosome( chromosome )
	{
		setup( stream ) ;
	}
	
	HapmapHaplotypesSNPDataSource::HapmapHaplotypesSNPDataSource( std::string const& filename, Chromosome chromosome )
		: m_filename( filename ),
		  m_compression_type( get_compression_type_indicated_by_filename( filename ) ),
		  m_number_of_samples( 0 ),
		  m_total_number_of_snps( 0 ),
		  m_chromosome( chromosome )
	{
		setup( filename, m_compression_type ) ; 
	}

	HapmapHaplotypesSNPDataSource::HapmapHaplotypesSNPDataSource( std::string const& filename, Chromosome chromosome, CompressionType compression_type )
		: m_filename( filename ),
		  m_compression_type( compression_type ),
		  m_number_of_samples( 0 ),
		  m_total_number_of_snps( 0 ),
		  m_chromosome( chromosome )
	{
		setup( filename, compression_type ) ;
	}

	void HapmapHaplotypesSNPDataSource::setup( std::string const& filename, CompressionType compression_type ) {
		m_stream_ptr = open_text_file_for_input( filename, compression_type ) ;
		read_header_data() ;				
		reset_to_start() ;
	}
	
	void HapmapHaplotypesSNPDataSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		m_stream_ptr = stream_ptr ;
		read_header_data() ;
		reset_to_start() ;
	}

	void HapmapHaplotypesSNPDataSource::read_header_data() {
		std::string line ;
		std::getline( *m_stream_ptr, line ) ;
		if( !*m_stream_ptr ) {
			throw MalformedInputError( m_filename, 0 ) ;
		}
		// deal with trailing space.
		if( line.size() > 0 && line[ line.size() - 1 ] == ' ' ) {
			line.resize( line.size() - 1 ) ;
		}
		using string_utils::slice ;
		std::vector< slice > elts = slice( line ).split( " " ) ;
		if( elts.size() < 2 || elts.size() % 2 != 0 ) {
			throw MalformedInputError( m_filename, 0 ) ;
		}
		
		else if( elts[0] != "rsID" ) {
			throw MalformedInputError( m_filename, 0, 0 ) ;
		}
		else if( elts[1] != "position_b36" ) {
			throw MalformedInputError( m_filename, 0, 1 ) ;
		}

		m_number_of_samples = ( elts.size() - 2 ) / 2 ;
		m_total_number_of_snps = count_lines_left_in_stream( *m_stream_ptr ) ;
	}

	void HapmapHaplotypesSNPDataSource::reset_to_start_impl() {
		stream().clear() ;
		stream().seekg( 0 ) ;
		if( !stream() ) {
			// oh, dear, seeking failed.  Reopen the file instead
			m_stream_ptr = open_text_file_for_input( m_filename, m_compression_type ) ;
		}
		std::string line ;
		std::getline( *m_stream_ptr, line ) ;
		if( !stream() ) {
			throw OperationFailedError( "genfile::HapmapHaplotypesSNPDataSource::reset_to_start_impl()", get_source_spec(), "reset to start" ) ;
		}
	}
	
	SNPDataSource::Metadata HapmapHaplotypesSNPDataSource::get_metadata() const {
		std::map< std::string, std::string > format ;
		format[ "ID" ] = "GT" ;
		format[ "Number" ] = "A" ;
		format[ "Description" ] = "Phased haplotype calls" ;
		SNPDataSource::Metadata result ;
		result.insert( std::make_pair( "FORMAT", format )) ;
		return result ;
	}

	void HapmapHaplotypesSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		using string_utils::slice ;
		std::getline( *m_stream_ptr, m_line ) ;
		// deal with trailing space.
		if( m_line.size() > 0 && m_line[ m_line.size() - 1 ] == ' ' ) {
			m_line.resize( m_line.size() - 1 ) ;
		}
		if( *this ) {
			m_elts = slice( m_line ).split( " " ) ;
			if( m_elts.size() != ( 2 + 2 * m_number_of_samples ) ) {
				throw MalformedInputError( m_filename, number_of_snps_read() ) ;
			}

			set_SNPID( "?" ) ;
			set_RSID( m_elts[0] ) ;
			set_chromosome( m_chromosome ) ;
			set_SNP_position( string_utils::to_repr< Position >( m_elts[1] ) ) ;
			std::string allele1 = "?" ;
			std::string allele2 = "?" ;

			std::size_t first_allele_i = 2 ;
			std::size_t second_allele_i = 2 ;
			// find the first allele not the same as the first_allele.
			for( ; second_allele_i != m_elts.size(); ++second_allele_i ) {
				if( m_elts[second_allele_i] != m_elts[first_allele_i] ) {
					break ;
				}
			}
			
			if( first_allele_i < m_elts.size() ) {
				m_first_allele = m_elts[ first_allele_i ] ;
				allele1 = m_elts[ first_allele_i ] ;
			}

			if( second_allele_i < m_elts.size() ) {
				m_second_allele = m_elts[ second_allele_i ] ;
				allele2 = m_elts[ second_allele_i ] ;
			}
			
			set_allele1( allele1 ) ;
			set_allele2( allele2 ) ;
		}
	}

	VariantDataReader::UniquePtr HapmapHaplotypesSNPDataSource::read_variant_data_impl() {
		throw OperationUnsupportedError( "genfile::HapmapHaplotypesSNPDataSource::read_variant_data()", "call", get_source_spec() ) ;
	}

	void HapmapHaplotypesSNPDataSource::ignore_snp_probability_data_impl() {
		m_line = "" ;
		m_elts.clear() ;
		m_first_allele = "" ;
		m_second_allele = "" ;
	}

}

