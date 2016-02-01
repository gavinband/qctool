
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <boost/bind.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/BGenFileSNPDataSource.hpp"
#include "genfile/zlib.hpp"

namespace genfile {
	BGenFileSNPDataSource::BGenFileSNPDataSource( std::string const& filename, Chromosome missing_chromosome ):
		m_filename( filename ),
		m_missing_chromosome( missing_chromosome )
	{
		setup( filename, get_compression_type_indicated_by_filename( filename )) ;
	}

	void BGenFileSNPDataSource::reset_to_start_impl() {
		stream().clear() ;
		stream().seekg(0) ;

		// read offset again and skip to first SNP data block
		bgen::uint32_t offset ;
		bgen::read_offset( (*m_stream_ptr), &offset ) ;
		m_stream_ptr->ignore( offset ) ;
	}

	SNPDataSource::Metadata BGenFileSNPDataSource::get_metadata() const {
		std::map< std::string, std::string > format ;
		format[ "ID" ] = "GP" ;
		format[ "Number" ] = "G" ;
		format[ "Type" ] = "Float" ;
		format[ "Description" ] = "Genotype call probabilities" ;
		SNPDataSource::Metadata result ;
		result.insert( std::make_pair( "FORMAT", format )) ;
		return result ;
	}

	void BGenFileSNPDataSource::get_sample_ids( GetSampleIds getter ) const {
		if( m_sample_ids ) {
			for( std::size_t i = 0; i < m_sample_ids->size(); ++i ) {
				getter( i, m_sample_ids->at(i) ) ;
			}
		}
	}

	void BGenFileSNPDataSource::read_snp_identifying_data_impl( VariantIdentifyingData* result ) {
		std::string SNPID, rsid, allele1, allele2 ;
		uint32_t position ;
		std::string chromosome_string ;
		if( bgen::read_snp_identifying_data( stream(), m_bgen_context, &SNPID, &rsid, &chromosome_string, &position, &allele1, &allele2 ) ) {
			Chromosome chr( chromosome_string ) ;
			if( chr.is_missing() ) {
				chr = m_missing_chromosome ;
			}
			*result = VariantIdentifyingData(
				SNPID, rsid, GenomePosition( chr, position ), allele1, allele2
			) ;
		}
	}

	namespace impl {
		struct BGenFileSNPDataReader: public VariantDataReader {
			BGenFileSNPDataReader( BGenFileSNPDataSource& source ):
				m_source( source )
			{
				assert( source ) ;
				bgen::read_genotype_data_block(
					m_source.stream(),
					m_source.bgen_context(),
					&m_source.m_compressed_data_buffer
				) ;
			}
			
			BGenFileSNPDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				assert( spec == "GP" || spec == ":genotypes:" ) ;
				bgen::uncompress_probability_data(
					m_source.bgen_context(),
					m_source.m_compressed_data_buffer,
					&(m_source.m_uncompressed_data_buffer)
				) ;
				bgen::parse_probability_data(
					&(m_source.m_uncompressed_data_buffer)[0],
					&(m_source.m_uncompressed_data_buffer)[0] + m_source.m_uncompressed_data_buffer.size(),
					m_source.bgen_context(),
					setter
				) ;
				return *this ;
			}
			
			bool supports( std::string const& spec ) const {
				return spec == "GP" || spec == ":genotypes:";
			}
			
			std::size_t get_number_of_samples() const {
				return m_source.number_of_samples() ;
			}

			void get_supported_specs( SpecSetter setter ) const {
				setter( "GP", "Float" ) ;
				//setter( ":genotypes:", "Float" ) ;
			}

		private:
			BGenFileSNPDataSource& m_source ;
			std::vector< double > m_genotypes ;
		} ;
	}

	VariantDataReader::UniquePtr BGenFileSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr( new impl::BGenFileSNPDataReader( *this )) ;
	}

	void BGenFileSNPDataSource::ignore_snp_probability_data_impl() {
		bgen::ignore_genotype_data_block(
			stream(),
			m_bgen_context
		) ;
	}

	void BGenFileSNPDataSource::setup( std::string const& filename, CompressionType compression_type ) {
		m_stream_ptr = open_binary_file_for_input( filename, compression_type ) ;
		bgen::uint32_t offset = 0 ;
		bgen::read_offset( *m_stream_ptr, &offset ) ;
		std::size_t bytes_read = bgen::read_header_block( *m_stream_ptr, &m_bgen_context ) ;

		if( offset < m_bgen_context.header_size() ) {
			throw FileStructureInvalidError() ;
		}
		
		if( m_bgen_context.flags & bgen::e_SampleIdentifiers ) {
			m_sample_ids = std::vector< std::string >() ;
			m_sample_ids->reserve( m_bgen_context.number_of_samples ) ;
			bytes_read += bgen::read_sample_identifier_block(
				*m_stream_ptr,
				m_bgen_context,
				boost::bind(
					&std::vector< std::string >::push_back,
					&(*m_sample_ids),
					_1
				)
			) ;
			assert( m_sample_ids->size() == m_bgen_context.number_of_samples ) ;
		}
		if( offset < bytes_read ) {
			throw FileStructureInvalidError() ;
		}
		// skip any remaining bytes before the first snp block
		m_stream_ptr->ignore( offset - bytes_read ) ;
		if( !*m_stream_ptr ) {
			throw FileStructureInvalidError() ;
		}
	}
}	

