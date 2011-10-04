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
		m_total_number_of_snps( 0 ),
		m_chromosome( chromosome ),
		m_have_chromosome_column( false )
	{
		setup( stream ) ;
	}
	
	GenFileSNPDataSource::GenFileSNPDataSource( std::string const& filename, Chromosome chromosome )
		: m_filename( filename ),
		  m_compression_type( get_compression_type_indicated_by_filename( filename ) ),
		  m_number_of_samples( 0 ),
		  m_total_number_of_snps( 0 ),
		  m_chromosome( chromosome ),
		  m_have_chromosome_column( false )
	{
		setup( filename, m_compression_type ) ; 
	}

	GenFileSNPDataSource::GenFileSNPDataSource(
		std::string const& filename,
		Chromosome chromosome,
		CompressionType compression_type,
		vcf::MetadataParser::Metadata const& metadata
	)
		: m_filename( filename ),
		  m_compression_type( compression_type ),
		  m_number_of_samples( 0 ),
		  m_total_number_of_snps( 0 ),
		  m_chromosome( chromosome ),
		  m_have_chromosome_column( false )
	{
		setup( filename, compression_type, metadata ) ;
	}

	void GenFileSNPDataSource::setup( std::string const& filename, CompressionType compression_type, vcf::MetadataParser::Metadata const& metadata ) {
		m_stream_ptr = open_text_file_for_input( filename, compression_type ) ;

		typedef vcf::MetadataParser::Metadata::const_iterator MetadataIterator ;
		std::pair< MetadataIterator, MetadataIterator > range = metadata.equal_range( "number-of-variants" ) ;
		if( range.first != range.second ) {
			std::size_t total_number_of_snps = 0 ;
			std::map< std::string, std::string >::const_iterator where = range.first->second.find( "" ) ;
			if( where == range.first->second.end() ) {
				throw MalformedInputError( "metadata", std::distance( metadata.begin(), range.first )) ;
			}
			else {
				total_number_of_snps = string_utils::to_repr< std::size_t >( where->second ) ;
			}
			if( (++range.first) != range.second ) {
				throw MalformedInputError( "metadata", std::distance( metadata.begin(), range.first )) ;
			}
			m_total_number_of_snps = total_number_of_snps ;
			read_header_data( false ) ;
		} else {
			read_header_data( true ) ;
		}
		reset_to_start() ;
	}

	void GenFileSNPDataSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		m_stream_ptr = stream_ptr ;
		read_header_data( true ) ;
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
	
	void GenFileSNPDataSource::read_snp_identifying_data_impl( 
		uint32_t* number_of_samples,
		std::string* SNPID,
		std::string* RSID,
		Chromosome* chromosome,
		uint32_t* SNP_position,
		char* allele1,
		char* allele2
	) {
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
	}

	namespace impl {
		struct GenFileSNPDataReader: public VariantDataReader {
			GenFileSNPDataReader( GenFileSNPDataSource& source )
			{
				uint32_t this_data_number_of_samples ;
				gen::impl::read_snp_probability_data(
					source.stream(),
					set_value( this_data_number_of_samples ),
					set_genotypes( m_genotypes )
				) ;
				assert( source ) ;
				assert(( m_genotypes.size() % 3 ) == 0 ) ;
			}
			
			GenFileSNPDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
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

	VariantDataReader::UniquePtr GenFileSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr( new impl::GenFileSNPDataReader( *this ) ) ;
	}

	void GenFileSNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		uint32_t this_data_number_of_samples ;
		gen::impl::read_snp_probability_data( stream(), set_value( this_data_number_of_samples ), set_genotype_probabilities ) ;
		if( *this ) {
			assert( this_data_number_of_samples == number_of_samples() ) ;
		}
	}

	void GenFileSNPDataSource::ignore_snp_probability_data_impl() {
		std::string line ;
		std::getline( stream(), line ) ;
	}

	void GenFileSNPDataSource::read_header_data( bool count_snps ) {
		try {
			// First let's have a look at the file.  If it is empty, we report 0 samples.
			m_stream_ptr->peek() ;
			int flags ;
			if( !m_stream_ptr->eof() ) {
				gen::read_header_information(
					*m_stream_ptr,
					set_value( m_total_number_of_snps ),
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

				if( count_snps ) {
					m_stream_ptr->clear() ;
					m_stream_ptr->seekg( 0 ) ;
					m_total_number_of_snps = gen::count_snp_blocks( *m_stream_ptr ) ;
				}
			}
		}
		catch( MalformedInputError const& e ) {
			throw MalformedInputError( m_filename, e.line(), e.column() ) ;
		}
	}
}

