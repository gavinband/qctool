#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/bgen.hpp"
#include "genfile/BGenFileSNPDataSource.hpp"

namespace genfile {
	BGenFileSNPDataSource::BGenFileSNPDataSource( std::string const& filename )
		: m_filename( filename )
	{
		setup( filename, get_compression_type_indicated_by_filename( filename )) ;
	}

	BGenFileSNPDataSource::BGenFileSNPDataSource( std::string const& filename, CompressionType compression_type )
		: m_filename( filename )
	{
		setup( filename, compression_type ) ;
	}

	void BGenFileSNPDataSource::reset_to_start_impl() {
		stream().clear() ;
		stream().seekg(0) ;

		// read offset again and skip to first SNP data block
		bgen::uint32_t offset ;
		bgen::read_offset( (*m_stream_ptr), &offset ) ;
		m_stream_ptr->ignore( offset ) ;
	}

	void BGenFileSNPDataSource::read_snp_identifying_data_impl( 
		uint32_t* number_of_samples,
		std::string* SNPID,
		std::string* RSID,
		Chromosome* chromosome,
		uint32_t* SNP_position,
		std::string* allele1,
		std::string* allele2
	) {
		unsigned char chr ;
		bgen::impl::read_snp_identifying_data( stream(), m_flags, number_of_samples, SNPID, RSID, &chr, SNP_position, allele1, allele2 ) ;
		*chromosome = Chromosome( ChromosomeEnum( chr ) ) ;
	}

	namespace impl {
		struct BGenFileSNPDataReader: public VariantDataReader {
			BGenFileSNPDataReader( BGenFileSNPDataSource& source )
			{
				if( source.m_flags & bgen::e_CompressedSNPBlocks ) {
					bgen::impl::read_compressed_snp_probability_data(
						source.stream(),
						source.m_flags,
						source.number_of_samples(),
						set_genotypes( m_genotypes )
					) ;
				}
				else {
					bgen::impl::read_snp_probability_data(
						source.stream(),
						source.m_flags,
						source.number_of_samples(),
						set_genotypes( m_genotypes )
					) ;
				}
				assert( source ) ;
				assert(( m_genotypes.size() % 3 ) == 0 ) ;
			}
			
			BGenFileSNPDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				std::size_t const N = m_genotypes.size() / 3 ;
				setter.set_number_of_samples( N ) ;
				for( std::size_t i = 0; i < N; ++i ) {
					setter.set_sample( i ) ;
					setter.set_number_of_entries( 3 ) ;
					for( std::size_t g = 0; g < 3; ++g ) {
						setter( m_genotypes[ 3*i + g ] ) ;
					}
				}
				return *this ;
			}
			
			bool supports( std::string const& spec ) const {
				return spec == "genotypes" ;
			}

			void get_supported_specs( SpecSetter setter ) const {
				setter( "genotypes", "Float" ) ;
			}

		private:
			std::vector< double > m_genotypes ;
		} ;
	}

	VariantDataReader::UniquePtr BGenFileSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr( new impl::BGenFileSNPDataReader( *this )) ;
	}

	void BGenFileSNPDataSource::ignore_snp_probability_data_impl() {
		if( m_flags & bgen::e_CompressedSNPBlocks ) {
			bgen::impl::read_compressed_snp_probability_data( stream(), m_flags, number_of_samples(), ignore() ) ;
		}
		else {
			bgen::impl::read_snp_probability_data( stream(), m_flags, number_of_samples(), ignore() ) ;
		}
	}

	void BGenFileSNPDataSource::setup( std::string const& filename, CompressionType compression_type ) {
		m_stream_ptr = open_binary_file_for_input( filename, compression_type ) ;
		bgen::uint32_t offset ;
		bgen::read_offset( (*m_stream_ptr), &offset ) ;
		bgen::uint32_t header_size = read_header_data() ;

		if( offset < header_size ) {
			throw FileStructureInvalidError() ;
		}
		// skip any remaining bytes before the first snp block
		m_stream_ptr->ignore( offset - header_size ) ;
		if( !*m_stream_ptr ) {
			throw FileStructureInvalidError() ;
		}
	}

	bgen::uint32_t BGenFileSNPDataSource::read_header_data() {
		bgen::uint32_t header_size ;

		bgen::read_header_block(
			(*m_stream_ptr),
			set_value( header_size ),
			set_value( m_total_number_of_snps ),
			set_value( m_number_of_samples ),
			ignore(),
			set_value( m_flags )
		) ;

		return header_size ;
	}
}	

