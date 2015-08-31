
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "../config.hpp"
#if HAVE_BOOST_FILESYSTEM
	#include <boost/filesystem.hpp>
	namespace BFS = boost::filesystem ;
#endif
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/ImputeHaplotypesSNPDataSource.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/FileUtils.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	ImputeHaplotypesSNPDataSource::ImputeHaplotypesSNPDataSource( std::string const& filename, Chromosome chromosome )
		: m_legend_filename( find_legend_file( filename ) ),
		  m_haplotypes_filename( filename ),
		  m_compression_type( get_compression_type_indicated_by_filename( m_haplotypes_filename ) ),
		  m_number_of_samples( 0 ),
		  m_chromosome( chromosome ),
		  m_good( true )
	{
		setup( filename, m_compression_type ) ; 
	}

	ImputeHaplotypesSNPDataSource::ImputeHaplotypesSNPDataSource( std::string const& filename, Chromosome chromosome, CompressionType compression_type )
		: m_legend_filename( find_legend_file( filename ) ),
		  m_haplotypes_filename( filename ),
		  m_compression_type( compression_type ),
		  m_number_of_samples( 0 ),
		  m_chromosome( chromosome ),
		  m_good( true )
	{
		setup( filename, m_compression_type ) ;
	}

	namespace impl {
		// borrowed this code from http://stackoverflow.com/questions/874134/find-if-string-endswith-another-string-in-c.
		bool has_ending( std::string const& a_string, std::string const& ending )
		{
			bool result = false ;
		    if( a_string.length() >= ending.length() ) {
		        result = ( 0 == a_string.compare( a_string.length() - ending.length(), ending.length(), ending ) ) ;
		    } 
			return result ;
		}
		
		std::vector< SNPIdentifyingData2 > read_legend( std::string const& filename, Chromosome const& chromosome, std::istream& stream ) {
			using string_utils::slice ;
			std::string line ;
			std::getline( stream, line ) ;
			if( !stream ) {
				throw MalformedInputError( filename, 0 ) ;
			}
			std::vector< slice > header_elts = slice( line ).split( " " ) ;
			if( header_elts.size() < 4 ) {
				throw MalformedInputError( filename, 0 ) ;
			}

			std::vector< SNPIdentifyingData2 > result ;
			while( std::getline( stream, line )) {
				std::vector< slice > elts = slice( line ).split( " " ) ;
				if( elts.size() != header_elts.size() ) {
					throw MalformedInputError( filename, result.size() + 1 ) ;
				}
				std::string allele1 = elts[2] ;
				std::string allele2 = elts[3] ;
				
				result.push_back(
					SNPIdentifyingData2(
						elts[0],
						GenomePosition( chromosome, string_utils::to_repr< Position >( elts[1] ) ),
						allele1,
						allele2
					)
				) ;
			}
			return result ;
		}
	}

	std::string ImputeHaplotypesSNPDataSource::find_legend_file( std::string const& haps_filename ) const {
		std::string result ;
#if HAVE_BOOST_FILESYSTEM
		// look for a similar-named file 
		if( impl::has_ending( haps_filename, ".haps" ) ) {
			result = haps_filename.substr( 0, haps_filename.size() - 5 ) ;
		} else if( impl::has_ending( haps_filename, ".haps.gz" ) ) {
			result = haps_filename.substr( 0, haps_filename.size() - 8 ) ;
		} else if( impl::has_ending( haps_filename, ".hap" ) ) {
			result = haps_filename.substr( 0, haps_filename.size() - 4 ) ;
		} else if( impl::has_ending( haps_filename, ".hap.gz" ) ) {
			result = haps_filename.substr( 0, haps_filename.size() - 7 ) ;
		}
		
		if( BFS::exists( result + ".legend" ) ) {
			result += ".legend" ;
		} else if( BFS::exists( result + ".legend.gz" )) {
			result += ".legend.gz" ;
		}
		else {
			throw genfile::InputError( result + ".legend[.gz]", "I could not find a legend file matching \"" + haps_filename + "\"." ) ;
		}
#else
		result = replace_or_add_extension( filename, ".legend" ) ;
#endif
		return result ;
	}


	void ImputeHaplotypesSNPDataSource::setup( std::string const& filename, CompressionType compression_type ) {
		std::auto_ptr< std::istream > legend_file = open_text_file_for_input( m_legend_filename ) ;
		m_snps = impl::read_legend( m_legend_filename, m_chromosome, *legend_file ) ;
		m_stream_ptr = open_text_file_for_input( m_haplotypes_filename, compression_type ) ;
		read_header_data() ;				
		reset_to_start() ;
	}
	
	void ImputeHaplotypesSNPDataSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		m_stream_ptr = stream_ptr ;
		read_header_data() ;
		reset_to_start() ;
	}

	void ImputeHaplotypesSNPDataSource::read_header_data() {
		std::string line ;
		std::getline( *m_stream_ptr, line ) ;
		if( *m_stream_ptr ) {
			// deal with trailing space which is sometimes found in IMPUTE haplotype files
			if( line.size() > 0 && line[ line.size() - 1 ] == ' ' ) {
				line.resize( line.size() - 1 ) ;
			}
			using string_utils::slice ;
			std::vector< slice > elts = slice( line ).split( " " ) ;
			if( elts.size() % 2 != 0 ) {
				throw MalformedInputError( m_haplotypes_filename, 0 ) ;
			}
		
			m_number_of_samples = elts.size() / 2 ;
		}
		else {
			m_number_of_samples = 0 ;
		}
	}

	void ImputeHaplotypesSNPDataSource::reset_to_start_impl() {
		stream().clear() ;
		stream().seekg( 0 ) ;
		if( !stream() ) {
			// oh, dear, seeking failed.  Reopen the file instead
			m_stream_ptr = open_text_file_for_input( m_haplotypes_filename, m_compression_type ) ;
		}
		if( !stream() ) {
			throw OperationFailedError( "genfile::ImputeHaplotypesSNPDataSource::reset_to_start_impl()", get_source_spec(), "reset to start" ) ;
		}
		m_good = true ;
	}
	
	SNPDataSource::Metadata ImputeHaplotypesSNPDataSource::get_metadata() const {
		std::map< std::string, std::string > format ;
		format[ "ID" ] = "GT" ;
		format[ "Type" ] = "String" ;
		format[ "Number" ] = "1" ;
		format[ "Description" ] = "Genotype" ;
		SNPDataSource::Metadata result ;
		result.insert( std::make_pair( "FORMAT", format )) ;
		return result ;
	}
	
	void ImputeHaplotypesSNPDataSource::get_snp_identifying_data_impl( 
		IntegerSetter const& set_number_of_samples,
		StringSetter const& set_SNPID,
		StringSetter const& set_RSID,
		ChromosomeSetter const& set_chromosome,
		SNPPositionSetter const& set_SNP_position,
		AlleleSetter const& set_allele1,
		AlleleSetter const& set_allele2
	) {
		if( number_of_snps_read() < m_snps.size() ) {
			SNPIdentifyingData2 const& snp = m_snps[ number_of_snps_read() ] ;
			set_number_of_samples( m_number_of_samples ) ;
			set_SNPID( snp.get_rsid() ) ;
			set_RSID( snp.get_rsid() ) ;
			set_chromosome( snp.get_position().chromosome() ) ;
			set_SNP_position( snp.get_position().position() ) ;
			set_allele1( snp.get_first_allele() ) ;
			set_allele2( snp.get_second_allele() ) ;
		} else {
			m_good = false ;
		}
	}

	namespace impl {
		struct ImputeHaplotypesSNPDataSourceReader: public VariantDataReader {
			typedef string_utils::slice slice ;

			ImputeHaplotypesSNPDataSourceReader(
				ImputeHaplotypesSNPDataSource const& source,
				std::size_t snp_index,
				std::string& line
			):
				m_source( source ),
				m_snp_index( snp_index )
			{
				line.swap( m_line ) ;
			}
			
			ImputeHaplotypesSNPDataSourceReader& get( std::string const& spec, PerSampleSetter& setter ) {
				if( m_elts.size() != 2 * m_source.number_of_samples() ) {
					m_elts = slice( m_line ).split( " " ) ;
				}
				if( m_elts.size() != 2 * m_source.number_of_samples() ) {
					throw MalformedInputError( m_source.get_source_spec(), m_snp_index ) ;
				}

				setter.set_number_of_samples( m_source.number_of_samples(), 2 ) ;
				for( std::size_t i = 0; i < m_source.number_of_samples(); ++i ) {
					if( m_elts[2*i].size() != 1 ) {
						throw MalformedInputError( m_source.get_source_spec(), m_snp_index, 2*i ) ;
					}
					if( m_elts[(2*i)+1].size() != 1 ) {
						throw MalformedInputError( m_source.get_source_spec(), m_snp_index, 2*i ) ;
					}
					setter.set_sample( i ) ;
					setter.set_number_of_entries( 2, ePerOrderedHaplotype, eAlleleIndex ) ;
					for( std::size_t j = 0; j < 2; ++j ) {
						try {
							if( m_elts[2*i+j][0] == '.' || m_elts[2*i+j] == "NA" ) {
								setter( MissingValue() ) ;
							} else {
								setter( string_utils::to_repr< Integer >( m_elts[ 2*i+j ] )) ;
							}
						}
						catch( string_utils::StringConversionError const& e ) {
							throw MalformedInputError( m_source.get_source_spec(), m_snp_index, 2*i + j ) ;
						}
					}
				}
				return *this ;
			}
			
			std::size_t get_number_of_samples() const { return m_source.number_of_samples() ; }

			bool supports( std::string const& spec ) const {
				return spec == "GT" || spec == ":genotypes:" ;
			}
			
			void get_supported_specs( SpecSetter setter ) const {
				setter( "GT", "Float" ) ;
				setter( ":genotypes:", "Float" ) ;
			}
			
			private:
				ImputeHaplotypesSNPDataSource const& m_source ;
				std::size_t const m_snp_index ;
				std::string m_line ;
				std::vector< slice > m_elts ;
		} ;
	}

	VariantDataReader::UniquePtr ImputeHaplotypesSNPDataSource::read_variant_data_impl() {
		assert( m_good ) ;
		std::size_t snp_index = number_of_snps_read() ;

		using string_utils::slice ;
		std::string line ;
		std::getline( *m_stream_ptr, line ) ;

		if( !*this ) {
			throw MalformedInputError( get_source_spec(), snp_index ) ;
		}

		return VariantDataReader::UniquePtr(
			new impl::ImputeHaplotypesSNPDataSourceReader(
				*this,
				snp_index,
				line
			)
		) ;
	}

	void ImputeHaplotypesSNPDataSource::ignore_snp_probability_data_impl() {
		std::string line ;
		std::getline( *m_stream_ptr, line ) ;
	}

}

