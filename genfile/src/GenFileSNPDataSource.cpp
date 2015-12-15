
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/GenFileSNPDataSource.hpp"

namespace genfile {
	GenFileSNPDataSource::GenFileSNPDataSource( std::auto_ptr< std::istream > stream, Chromosome chromosome ):
		m_filename( "(unnamed stream)" ),
		m_compression_type( "no_compression" ),
		m_number_of_samples( 0 ),
		m_chromosome( chromosome ),
		m_have_chromosome_column( false )
	{
		setup( stream ) ;
	}
	
	GenFileSNPDataSource::GenFileSNPDataSource( std::string const& filename, Chromosome chromosome )
		: m_filename( filename ),
		  m_compression_type( get_compression_type_indicated_by_filename( filename ) ),
		  m_number_of_samples( 0 ),
		  m_chromosome( chromosome ),
		  m_have_chromosome_column( false )
	{
		setup( filename, m_compression_type ) ; 
	}

	void GenFileSNPDataSource::setup( std::string const& filename, CompressionType compression_type ) {
		m_stream_ptr = open_text_file_for_input( filename, compression_type ) ;
		read_header_data() ;
		reset_to_start() ;
	}

	void GenFileSNPDataSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		m_stream_ptr = stream_ptr ;
		read_header_data() ;
		reset_to_start() ;
	}

	void GenFileSNPDataSource::reset_to_start_impl() {
		stream().clear() ;
		stream().seekg( 0 ) ;
		if( !stream() ) {
			// oh, dear, seeking failed.  Reopen the file instead
			m_stream_ptr = open_text_file_for_input( m_filename, m_compression_type ) ;
		}
		if( !stream() ) {
			throw OperationFailedError( "genfile::GenFileSNPDataSource::reset_to_start_impl()", get_source_spec(), "reset to start" ) ;
		}
	}
	
	SNPDataSource::Metadata GenFileSNPDataSource::get_metadata() const {
		std::map< std::string, std::string > format ;
		format[ "ID" ] = "GP" ;
		format[ "Number" ] = "G" ;
		format[ "Type" ] = "Float" ;
		format[ "Description" ] = "Genotype call probabilities" ;
		SNPDataSource::Metadata result ;
		result.insert( std::make_pair( "FORMAT", format )) ;
		return result ;
	}
	
	void GenFileSNPDataSource::read_snp_identifying_data_impl( 
		uint32_t* number_of_samples,
		std::string* SNPID,
		std::string* RSID,
		Chromosome* chromosome,
		uint32_t* SNP_position,
		std::string* allele1,
		std::string* allele2
	) {
		try {
			if( m_have_chromosome_column ) {
				gen::impl::read_snp_identifying_data( stream(), chromosome, SNPID, RSID, SNP_position, allele1, allele2 ) ;
				if( *this ) {
					*number_of_samples = m_number_of_samples ;
				}
			} else {
				gen::impl::read_snp_identifying_data( stream(), SNPID, RSID, SNP_position, allele1, allele2 ) ;
				if( *this ) {
					*number_of_samples = m_number_of_samples ;
					*chromosome = m_chromosome ;
				}
			}
		} catch( InputError const& e ) {
			throw MalformedInputError(
				m_filename,
				e.message(),
				number_of_snps_read()
			) ;
		}
	}

	namespace impl {
		struct GenFileSNPDataReader: public VariantDataReader {
			GenFileSNPDataReader( GenFileSNPDataSource& source ) {
				try {
					uint32_t this_data_number_of_samples ;
					gen::impl::read_snp_probability_data(
						source.stream(),
						set_value( this_data_number_of_samples ),
						set_genotypes( m_genotypes ),
						source.m_line
					) ;
					if( !source ) {
						throw genfile::MalformedInputError( source.m_filename, source.number_of_snps_read() ) ;
					}
					assert(( m_genotypes.size() % 3 ) == 0 ) ;
				} catch( genfile::InputError const& e ) {
					throw genfile::MalformedInputError( source.m_filename, e.message(), source.number_of_snps_read() ) ;
				}
			}
			
			GenFileSNPDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				if( spec != "GP" && spec != ":genotypes:" ) {
					throw BadArgumentError(
						"genfile::GenFileSNPDataReader::get()",
						"spec=\"" + spec + "\"",
						"Only \"GP\" and \":genotypes:\" are supported in a GEN file."
					) ;
				}
				std::size_t const N = m_genotypes.size() / 3 ;
				setter.initialise( N, 2 ) ;
				uint32_t ploidy = 2 ;
				for( std::size_t i = 0; i < N; ++i ) {
					setter.set_sample( i ) ;
					setter.set_number_of_entries( ploidy,  3, ePerUnorderedGenotype, eProbability ) ;
					for( std::size_t g = 0; g < 3; ++g ) {
						setter.set_value( g, m_genotypes[ 3*i + g ] ) ;
					}
				}
				setter.finalise() ;
				return *this ;
			}
			
			std::size_t get_number_of_samples() const {
				return m_genotypes.size() / 3 ; 
			}
			
			bool supports( std::string const& spec ) const {
				return spec == "GP" || spec == ":genotypes:";
			}
			
			void get_supported_specs( SpecSetter setter ) const {
				setter( "GP", "Float" ) ;
				setter( ":genotypes:", "Float" ) ;
			}
			
		private:
			std::vector< double > m_genotypes ;
			std::string m_buffer ;
		} ;
	}

	VariantDataReader::UniquePtr GenFileSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr( new impl::GenFileSNPDataReader( *this ) ) ;
	}

	void GenFileSNPDataSource::ignore_snp_probability_data_impl() {
		std::string line ;
		std::getline( stream(), line ) ;
	}

	void GenFileSNPDataSource::read_header_data() {
		try {
			// First let's have a look at the file.  If it is empty, we report 0 samples.
			m_stream_ptr->peek() ;
			int flags ;
			if( !m_stream_ptr->eof() ) {
				gen::read_header_information(
					*m_stream_ptr,
					set_value( m_number_of_samples ),
					set_value( flags )
				) ;
				if( !(*m_stream_ptr )) {
					throw MalformedInputError( m_filename, 1 ) ;
				}

				if( flags & 0x1 ) {
					m_have_chromosome_column = true ;
				} else {
					m_have_chromosome_column = false ;
				}
			}
		}
		catch( MalformedInputError const& e ) {
			throw MalformedInputError( m_filename, e.line(), e.column() ) ;
		}
	}
}

