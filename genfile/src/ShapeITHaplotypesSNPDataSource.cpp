#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/ShapeITHaplotypesSNPDataSource.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/FileUtils.hpp"
#include "genfile/VariantDataReader.hpp"

namespace genfile {
	ShapeITHaplotypesSNPDataSource::ShapeITHaplotypesSNPDataSource( std::string const& filename, Chromosome chromosome ):
		m_haplotypes_filename( filename ),
		m_compression_type( get_compression_type_indicated_by_filename( m_haplotypes_filename ) ),
		m_number_of_samples( 0 ),
		m_have_chromosome_column( false ),
		m_chromosome( chromosome ),
		m_good( true )
	{
		setup( filename, m_compression_type ) ; 
	}

	ShapeITHaplotypesSNPDataSource::ShapeITHaplotypesSNPDataSource( std::string const& filename, Chromosome chromosome, CompressionType compression_type ):
		m_haplotypes_filename( filename ),
		m_compression_type( compression_type ),
		m_number_of_samples( 0 ),
		m_have_chromosome_column( false ),
		m_chromosome( chromosome ),
		m_good( true )
	{
		setup( filename, m_compression_type ) ;
	}

	void ShapeITHaplotypesSNPDataSource::setup( std::string const& filename, CompressionType compression_type ) {
		m_stream_ptr = open_text_file_for_input( m_haplotypes_filename, compression_type ) ;
		count_samples() ;				
		reset_to_start() ;
	}
	
	void ShapeITHaplotypesSNPDataSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		m_stream_ptr = stream_ptr ;
		count_samples() ;
		reset_to_start() ;
	}

	void ShapeITHaplotypesSNPDataSource::count_samples() {
		std::string line ;
		std::getline( *m_stream_ptr, line ) ;
		std::vector< string_utils::slice > elts = string_utils::slice( line ).split( " \t" ) ;

		if( ( elts.size() - 5 ) % 2 == 0 ) {
			m_have_chromosome_column = false ;
			m_number_of_samples = ( elts.size() - 5 ) / 2 ;
		} else if( ( elts.size() - 6 ) % 2 == 0 ) {
			m_have_chromosome_column = true ;
			m_number_of_samples = ( elts.size() - 6 ) / 2 ;
		} else {
			throw MalformedInputError( get_source_spec(), 0 ) ;
		}
	}

	void ShapeITHaplotypesSNPDataSource::reset_to_start_impl() {
		stream().clear() ;
		stream().seekg( 0 ) ;
		if( !stream() ) {
			// oh, dear, seeking failed.  Reopen the file instead
			m_stream_ptr = open_text_file_for_input( m_haplotypes_filename, m_compression_type ) ;
		}
		if( !stream() ) {
			throw OperationFailedError( "genfile::ShapeITHaplotypesSNPDataSource::reset_to_start_impl()", get_source_spec(), "reset to start" ) ;
		}
		m_good = true ;
	}
	
	void ShapeITHaplotypesSNPDataSource::read_snp_identifying_data_impl( 
		uint32_t* number_of_samples,
		std::string* SNPID,
		std::string* RSID,
		Chromosome* chromosome,
		uint32_t* SNP_position,
		std::string* allele1,
		std::string* allele2
	) {
		SNPIdentifyingData snp ;
		if( m_have_chromosome_column ) {
			std::string chromosome_string ;
			(*m_stream_ptr) >> chromosome_string ;
			snp.position().chromosome() = Chromosome( chromosome_string ) ;
		}
		(*m_stream_ptr) >> snp.SNPID() >> snp.rsid() >> snp.position().position() >> snp.first_allele() >> snp.second_allele() ;

		if( *m_stream_ptr ) {
			*number_of_samples = m_number_of_samples ;
			*SNPID = snp.get_SNPID() ;
			*RSID = snp.get_rsid() ;
			*chromosome = snp.get_position().chromosome() ;
			*SNP_position = snp.get_position().position() ;
			*allele1 = snp.get_first_allele() ;
			*allele2 = snp.get_second_allele() ;
		}
		else {
			m_good = false ;
		}
	}

	namespace impl {
		struct ShapeITHaplotypesSNPDataSourceReader: public VariantDataReader {
			typedef string_utils::slice slice ;

			ShapeITHaplotypesSNPDataSourceReader(
				ShapeITHaplotypesSNPDataSource const& source,
				std::size_t snp_index,
				std::string& line
			):
				m_source( source ),
				m_snp_index( snp_index )
			{
				line.swap( m_line ) ;
			}
			
			ShapeITHaplotypesSNPDataSourceReader& get( std::string const& spec, PerSampleSetter& setter ) {
				std::size_t const N = m_source.number_of_samples() ;
				if( m_elts.size() != 2*N ) {
					m_elts = slice( m_line ).split( " " ) ;
				}
				assert( m_elts[0].size() > 0 ) ;
				if( m_elts.size() != 2*N ) {
					throw MalformedInputError( m_source.get_source_spec(), m_snp_index ) ;
				}

				setter.set_number_of_samples( N ) ;
				setter.set_order_type( vcf::EntriesSetter::eOrderedList ) ;

				for( std::size_t i = 0; i < N; ++i ) {
						setter.set_sample( i ) ;
						setter.set_number_of_entries( 2 ) ;
						for( std::size_t j = 0; j < 2; ++j ) {
							try {
								if( m_elts[ 2*i+j ] == "." || m_elts[ 2*i+j ] == "NA" ) {
									setter( MissingValue() ) ;
								} else {
									setter( string_utils::to_repr< Integer >( m_elts[ 2*i+j ] )) ;
								}
							}
							catch( string_utils::StringConversionError const& e ) {
								throw MalformedInputError( m_source.get_source_spec(), m_snp_index, m_source.get_number_of_id_columns() + 2*i + j ) ;
							}
						}
				}
				return *this ;
			}
			
			bool supports( std::string const& spec ) const {
				return spec == "genotypes" ;
			}
			
			void get_supported_specs( SpecSetter setter ) const {
				setter( "genotypes", "Integer" ) ;
			}
			
			private:
				ShapeITHaplotypesSNPDataSource const& m_source ;
				std::size_t const m_snp_index ;
				std::string m_line ;
				std::vector< slice > m_elts ;
		} ;
	}
	
	std::size_t ShapeITHaplotypesSNPDataSource::get_number_of_id_columns() const {
		return m_have_chromosome_column ? 6 : 5 ;
	}

	VariantDataReader::UniquePtr ShapeITHaplotypesSNPDataSource::read_variant_data_impl() {
		assert( m_good ) ;
		std::size_t snp_index = number_of_snps_read() ;

		std::string line ;
		if( m_stream_ptr->peek() == ' ' ) {
			m_stream_ptr->get() ;
		}
		
		std::getline( *m_stream_ptr, line ) ;

		if( !*this ) {
			throw MalformedInputError( get_source_spec(), snp_index ) ;
		}

		return VariantDataReader::UniquePtr(
			new impl::ShapeITHaplotypesSNPDataSourceReader(
				*this,
				snp_index,
				line
			)
		) ;
	}

	void ShapeITHaplotypesSNPDataSource::ignore_snp_probability_data_impl() {
		std::string line ;
		std::getline( *m_stream_ptr, line ) ;
	}

}

