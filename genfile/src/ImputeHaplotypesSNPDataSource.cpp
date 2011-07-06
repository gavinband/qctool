#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/ImputeHaplotypesSNPDataSource.hpp"
#include "genfile/string_utils/slice.hpp"
#include "genfile/FileUtils.hpp"

namespace genfile {
	ImputeHaplotypesSNPDataSource::ImputeHaplotypesSNPDataSource( std::string const& filename, Chromosome chromosome )
		: m_legend_filename( replace_or_add_extension( filename, ".legend" ) ),
		  m_haplotypes_filename( filename ),
		  m_compression_type( get_compression_type_indicated_by_filename( m_haplotypes_filename ) ),
		  m_number_of_samples( 0 ),
		  m_chromosome( chromosome ),
		  m_good( true )
	{
		setup( filename, m_compression_type ) ; 
	}

	ImputeHaplotypesSNPDataSource::ImputeHaplotypesSNPDataSource( std::string const& filename, Chromosome chromosome, CompressionType compression_type )
		: m_legend_filename( replace_or_add_extension( filename, ".legend" ) ),
		  m_haplotypes_filename( filename ),
		  m_compression_type( compression_type ),
		  m_number_of_samples( 0 ),
		  m_chromosome( chromosome ),
		  m_good( true )
	{
		setup( filename, m_compression_type ) ;
	}

	namespace impl {
		std::vector< SNPIdentifyingData > read_legend( std::string const& filename, Chromosome const& chromosome, std::istream& stream ) {
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

			std::vector< SNPIdentifyingData > result ;
			while( std::getline( stream, line )) {
				std::vector< slice > elts = slice( line ).split( " " ) ;
				if( elts.size() != header_elts.size() ) {
					throw MalformedInputError( filename, result.size() + 1 ) ;
				}
				char allele1 = '?' ;
				char allele2 = '?' ;
				if( elts[2].size() == 1 ) {
					allele1 = elts[2][0] ;
				}
				if( elts[3].size() == 1 ) {
					allele2 = elts[3][0] ;
				}
				
				result.push_back(
					SNPIdentifyingData(
						"?",
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
		std::size_t total_number_of_snps = 0 ;
		std::string line ;
		std::getline( *m_stream_ptr, line ) ;
		if( *m_stream_ptr ) {
			++total_number_of_snps ;
			// deal with trailing space.
			if( line.size() > 0 && line[ line.size() - 1 ] == ' ' ) {
				line.resize( line.size() - 1 ) ;
			}
			using string_utils::slice ;
			std::vector< slice > elts = slice( line ).split( " " ) ;
			if( elts.size() % 2 != 0 ) {
				throw MalformedInputError( m_haplotypes_filename, 0 ) ;
			}
		
			m_number_of_samples = elts.size() / 2 ;
			total_number_of_snps += count_lines_left_in_stream( *m_stream_ptr ) + 1 ;
		}
		else {
			m_number_of_samples = 0 ;
		}
		if( !total_number_of_snps == m_snps.size() ) {
			throw MalformedInputError( m_haplotypes_filename, m_snps.size() ) ;
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
			throw OperationFailedError() ;
		}
		m_good = true ;
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
			SNPIdentifyingData const& snp = m_snps[ number_of_snps_read() ] ;
			set_number_of_samples( m_number_of_samples ) ;
			set_SNPID( snp.get_SNPID() ) ;
			set_RSID( snp.get_rsid() ) ;
			set_chromosome( snp.get_position().chromosome() ) ;
			set_SNP_position( snp.get_position().position() ) ;
			set_allele1( snp.get_first_allele() ) ;
			set_allele2( snp.get_second_allele() ) ;
		} else {
			m_good = false ;
		}
	}

	VariantDataReader::UniquePtr ImputeHaplotypesSNPDataSource::read_variant_data_impl() {
		throw OperationUnsupportedError( "genfile::ImputeHaplotypesSNPDataSource::read_variant_data()", get_source_spec() ) ;
	}


	void ImputeHaplotypesSNPDataSource::read_snp_probability_data_impl(
		GenotypeProbabilitySetter const& set_genotype_probabilities
	) {
		if( !m_good ) {
			return ;
		}

		using string_utils::slice ;
		std::string line ;
		std::getline( *m_stream_ptr, line ) ;
		if( !*this ) {
			throw MalformedInputError( m_haplotypes_filename, number_of_snps_read() ) ;
		}
		std::vector< slice > elts = slice( line ).split( " " ) ;
		if( elts.size() != 2 * m_number_of_samples ) {
			throw MalformedInputError( m_haplotypes_filename, number_of_snps_read() ) ;
		}

		for( std::size_t i = 0; i < m_number_of_samples; ++i ) {
			if( elts[2*i].size() != 1 ) {
				throw MalformedInputError( m_haplotypes_filename, number_of_snps_read(), 2*i ) ;
			}
			if( elts[(2*i)+1].size() != 1 ) {
				throw MalformedInputError( m_haplotypes_filename, number_of_snps_read(), 2*i ) ;
			}
			char allele1 = elts[2*i][0] ;
			char allele2 = elts[(2*i)+1][0] ;
			std::size_t A_allele_count = std::size_t( allele1 == '0' ) + std::size_t( allele2 == '0' ) ;
			std::size_t B_allele_count = std::size_t( allele1 == '1' ) + std::size_t( allele2 == '1' ) ;
			if( A_allele_count == 2 ) {
				set_genotype_probabilities( i, 1.0, 0.0, 0.0 ) ;
			}
			else if( A_allele_count == 1 && B_allele_count == 1 ) {
				set_genotype_probabilities( i, 0.0, 1.0, 0.0 ) ;
			}
			else if( B_allele_count == 2 ) {
				set_genotype_probabilities( i, 0.0, 0.0, 1.0 ) ;
			}
			else {
				// In any other case, did not recognise the allele, treat it as missing.
				// Actually
				if( elts[ 2*i ][0] != '0' && elts[ 2*i ][0] != '1' ) {
					throw MalformedInputError( m_haplotypes_filename, number_of_snps_read(), (2*i) ) ;
				}
				else {
					throw MalformedInputError( m_haplotypes_filename, number_of_snps_read(), (2*i)+1 ) ;
				}
				// This is what we would do if we didn't throw, which we do, so we don't do this.
				// set_genotype_probabilities( i, 0.0, 0.0, 0.0 ) ;
			}
		}
	}

	void ImputeHaplotypesSNPDataSource::ignore_snp_probability_data_impl() {
		std::string line ;
		std::getline( *m_stream_ptr, line ) ;
	}

}

