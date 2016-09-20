
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/ImputeHapProbsSNPDataSource.hpp"

namespace genfile {
	ImputeHapProbsSNPDataSource::ImputeHapProbsSNPDataSource( std::auto_ptr< std::istream > stream, Chromosome chromosome ):
		m_filename( "(unnamed stream)" ),
		m_compression_type( "no_compression" ),
		m_number_of_samples( 0 ),
		m_chromosome( chromosome ),
		m_have_chromosome_column( false )
	{
		setup( stream ) ;
	}
	
	ImputeHapProbsSNPDataSource::ImputeHapProbsSNPDataSource( std::string const& filename, Chromosome chromosome )
		: m_filename( filename ),
		  m_compression_type( get_compression_type_indicated_by_filename( filename ) ),
		  m_number_of_samples( 0 ),
		  m_chromosome( chromosome ),
		  m_have_chromosome_column( false )
	{
		setup( filename, m_compression_type ) ; 
	}

	void ImputeHapProbsSNPDataSource::setup( std::string const& filename, CompressionType compression_type ) {
		m_stream_ptr = open_text_file_for_input( filename, compression_type ) ;
		read_header_data() ;
		reset_to_start() ;
	}

	void ImputeHapProbsSNPDataSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		m_stream_ptr = stream_ptr ;
		read_header_data() ;
		reset_to_start() ;
	}

	void ImputeHapProbsSNPDataSource::reset_to_start_impl() {
		stream().clear() ;
		stream().seekg( 0 ) ;
		if( !stream() ) {
			// oh, dear, seeking failed.  Reopen the file instead
			m_stream_ptr = open_text_file_for_input( m_filename, m_compression_type ) ;
		}
		if( !stream() ) {
			throw OperationFailedError( "genfile::ImputeHapProbsSNPDataSource::reset_to_start_impl()", get_source_spec(), "reset to start" ) ;
		}
	}
	
	SNPDataSource::Metadata ImputeHapProbsSNPDataSource::get_metadata() const {
		std::map< std::string, std::string > format ;
		format[ "ID" ] = "HP" ;
		format[ "Number" ] = "2" ;
		format[ "Type" ] = "Float" ;
		format[ "Description" ] = "Per-haplotype allele probabilities" ;
		SNPDataSource::Metadata result ;
		result.insert( std::make_pair( "FORMAT", format )) ;
		return result ;
	}
	
	void ImputeHapProbsSNPDataSource::read_snp_identifying_data_impl( VariantIdentifyingData* result ) {
		using string_utils::slice ;
		using string_utils::to_string ;
		using string_utils::to_repr ;
		if( std::getline( *m_stream_ptr, m_line ) ) {
			m_elts = string_utils::slice( m_line ).split( " " ) ;
			std::size_t const expected = (m_number_of_samples*2) + 5 +  (m_have_chromosome_column ? 1 : 0) ;
			if( m_elts.size() != expected ) {
				throw genfile::MalformedInputError(
					m_filename,
					"Wrong number of entries in line (" + to_string( m_elts.size() )  + ", expected " + to_string( expected ) + ")",
					number_of_snps_read()
				) ;
			}

			if( m_have_chromosome_column ) {
				*result = VariantIdentifyingData(
					m_elts[1],
					m_elts[2],
					GenomePosition( Chromosome( m_elts[0] ), to_repr< uint32_t >( m_elts[3]) ),
					m_elts[4],
					m_elts[5]
				) ;
			} else {
				*result = VariantIdentifyingData(
					m_elts[0],
					m_elts[1],
					GenomePosition( m_chromosome, to_repr< uint32_t >( m_elts[2]) ),
					m_elts[3],
					m_elts[4]
				) ;
			}
		}
	}

	namespace impl {
		struct ImputeHapProbsSNPDataReader: public VariantDataReader {
			ImputeHapProbsSNPDataReader( ImputeHapProbsSNPDataSource& source ):
				m_elts( source.m_elts ),
				m_have_chromosome_column( source.m_have_chromosome_column ),
				m_number_of_samples( source.m_number_of_samples )
			{
			}
			
			ImputeHapProbsSNPDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				using string_utils::to_repr ;
				if( spec != "HP" && spec != ":genotypes:" ) {
					throw BadArgumentError(
						"genfile::ImputeHapProbsSNPDataReader::get()",
						"spec=\"" + spec + "\"",
						"Only \"HP\" and \":genotypes:\" are supported in a IMPUTE haplotype probabilities file."
					) ;
				}
				std::size_t const preamble_cols = 5 + (m_have_chromosome_column ? 1 : 0) ;
				std::size_t const N = (m_elts.size() - preamble_cols) / 2 ;
				assert( N == m_number_of_samples ) ;
				setter.initialise( N, 2 ) ;
				uint32_t ploidy = 2 ;
				for( std::size_t i = 0; i < N; ++i ) {
					setter.set_sample( i ) ;
					setter.set_number_of_entries( ploidy, 4, ePerPhasedHaplotypePerAllele, eProbability ) ;
					double a0 = to_repr< double >( m_elts[2*i+preamble_cols] ) ;
					double a1 = to_repr< double >( m_elts[2*i+preamble_cols+1] ) ;
					setter.set_value( 0, 1.0 - a0 ) ;
					setter.set_value( 1, a0 ) ;
					setter.set_value( 2, 1.0 - a1 ) ;
					setter.set_value( 3, a1 ) ;
				}
				setter.finalise() ;
				return *this ;
			}
			
			std::size_t get_number_of_samples() const {
				return m_number_of_samples ;
			}
			
			bool supports( std::string const& spec ) const {
				return spec == "HP" || spec == ":genotypes:";
			}
			
			void get_supported_specs( SpecSetter setter ) const {
				setter( "HP", "Float" ) ;
				setter( ":genotypes:", "Float" ) ;
			}
			
		private:
			std::vector< string_utils::slice > const& m_elts ;
			bool m_have_chromosome_column ;
			std::size_t const m_number_of_samples ;
		} ;
	}

	VariantDataReader::UniquePtr ImputeHapProbsSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr( new impl::ImputeHapProbsSNPDataReader( *this ) ) ;
	}

	void ImputeHapProbsSNPDataSource::ignore_snp_probability_data_impl() {
		std::string line ;
		std::getline( stream(), line ) ;
	}

	void ImputeHapProbsSNPDataSource::read_header_data() {
		using string_utils::slice ;
		using string_utils::to_string ;
		try {
			// First let's have a look at the file.  If it is empty, we report 0 samples.
			m_stream_ptr->peek() ;
			int flags ;
			if( !m_stream_ptr->eof() ) {
				std::getline( *m_stream_ptr, m_line ) ;
				m_elts = slice( m_line ).split( " " ) ;
				if( m_elts.size() < 5 ) {
					throw genfile::MalformedInputError(
						m_filename,
						"Expected at least 5 columns (found " + to_string( m_elts.size() ) + ").",
						0
					) ;
				}
				if( (m_elts.size() - 6) % 2 == 0 ) {
					m_have_chromosome_column = true ;
					m_number_of_samples = (m_elts.size()-6) / 2 ;
				} else if( (m_elts.size() - 5) % 2 == 0 ) {
					m_have_chromosome_column = false ;
					m_number_of_samples = (m_elts.size()-5) / 2 ;
				} else {
					throw genfile::MalformedInputError(
						m_filename,
						"Expected 2n+5 or 2n+6 columns (found " + to_string( m_elts.size() ) + ").",
						0
					) ;
				}
				if( !(*m_stream_ptr )) {
					throw MalformedInputError( m_filename, 1 ) ;
				}
			}
		}
		catch( MalformedInputError const& e ) {
			throw MalformedInputError( m_filename, e.line(), e.column() ) ;
		}
	}
}

