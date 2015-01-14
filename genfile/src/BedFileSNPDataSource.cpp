
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <boost/format.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/BedFileSNPDataSource.hpp"

#define DEBUG_BED_FORMAT 1

namespace genfile {
	BedFileSNPDataSource::BedFileSNPDataSource( std::string const& bedFilename, std::string const& bimFilename, std::string const& famFilename ):
		m_bed_filename( bedFilename ),
		m_bim_filename( bimFilename ),
		m_exhausted( false )
	{
		m_genotype_table.push_back( std::make_pair( 0, 0 )) ;
		m_genotype_table.push_back( std::make_pair( -1, -1 )) ;
		m_genotype_table.push_back( std::make_pair( 0, 1 )) ;
		m_genotype_table.push_back( std::make_pair( 1, 1 )) ;
		setup( bedFilename, bimFilename, famFilename ) ;
	}

	void BedFileSNPDataSource::setup( std::string const& bedFilename, std::string const& bimFilename, std::string const& famFilename ) {
		m_bed_stream_ptr = open_binary_file_for_input( bedFilename ) ;
		m_bim_stream_ptr = open_text_file_for_input( bimFilename ) ;

		{
			std::size_t sampleCount = 0 ;
			std::auto_ptr< std::istream > fam = open_text_file_for_input( famFilename ) ;
			std::string line ;
			while( std::getline( *fam, line )) {
				++sampleCount ;
			}
			m_number_of_samples = sampleCount ;
		}

		char header[3] ;
		m_bed_stream_ptr->read( &header[0], 3 ) ;
		if( header[0] != 108 || header[1] != 27 ) {
#if DEBUG_BED_FORMAT
			std::cerr << "BedFileSNPDataSource::setup(): first few bytes of BED file are: " ;
			boost::format fmt( "%02x" ) ;
			std::cerr << fmt % int( header[0] ) << fmt % int( header[1] ) << fmt % int( header[2] ) << ".\n" ; 
#endif
			throw genfile::MalformedInputError(
				bedFilename,
				"File does not appear to be a binary PED file (according to magic number).  See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml.",
				0
			) ;
		}
		if( header[2] != 1 ) {
			throw genfile::MalformedInputError(
				bedFilename,
				"File has mode " + genfile::string_utils::to_string( int( header[2] ) ) + "; only SNP-major mode (mode = 1) is supported.",
				0
			) ;
		}
	}

	void BedFileSNPDataSource::reset_to_start_impl() {
		m_bim_stream_ptr->clear() ;
		m_bim_stream_ptr->seekg(0) ;
		m_bed_stream_ptr->clear() ;
		m_bed_stream_ptr->seekg(0) ;
		// skip the magic number and mode
		m_bed_stream_ptr->ignore( 3 ) ;
		m_exhausted = false ;
	}

	SNPDataSource::Metadata BedFileSNPDataSource::get_metadata() const {
		std::map< std::string, std::string > format ;
		format[ "ID" ] = "GT" ;
		format[ "Number" ] = "G" ;
		format[ "Description" ] = "Genotype calls" ;
		SNPDataSource::Metadata result ;
		result.insert( std::make_pair( "FORMAT", format )) ;
		return result ;
	}

	void BedFileSNPDataSource::read_snp_identifying_data_impl(
		uint32_t* number_of_samples,
		std::string* SNPID,
		std::string* RSID,
		Chromosome* chromosome,
		uint32_t* SNP_position,
		std::string* allele1,
		std::string* allele2
	) {
		std::string chromosome_string ;
		std::string rsid ;
		uint32_t position ;
		std::string alleleA, alleleB ;
		double cM ;
		(*m_bim_stream_ptr) >> chromosome_string >> rsid >> cM >> position >> alleleA >> alleleB ;
		m_bim_stream_ptr->ignore( std::numeric_limits< std::streamsize >::max(), '\n' ) ;
		
		if( *m_bim_stream_ptr ) {
			*SNPID = rsid ;
			*RSID = rsid ;
			*chromosome = Chromosome( chromosome_string ) ;
			*SNP_position = position ;
			*allele1 = alleleA ;
			*allele2 = alleleB ;
		} else {
			m_exhausted = true ;
		}
	}

	namespace impl {
		struct BedFileSNPDataReader: public VariantDataReader {
			BedFileSNPDataReader( BedFileSNPDataSource& source ):
				m_source( source ),
				m_number_of_samples( source.number_of_samples() ),
				m_buffer( ( m_number_of_samples + 3 ) / 4, 0 )
			{
				source.m_bed_stream_ptr->read( &m_buffer[0], m_buffer.size() ) ;
				if( !(*source.m_bed_stream_ptr )) {
					throw MalformedInputError(
						source.m_bed_filename,
						"Unable to read genotypes from file",
						source.number_of_snps_read()
					) ;
				}
			}
			
			BedFileSNPDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
                assert( spec == "GT" || spec == ":genotypes:" ) ;
				setter.set_number_of_samples( m_number_of_samples ) ;
				for( std::size_t i = 0; i < m_number_of_samples; ++i ) {
					std::size_t index = i/4 ;
					std::size_t data = ( m_buffer[ index ] >> (2*(i%4)) ) & 0x3 ;

					setter.set_sample( i ) ;
					setter.set_number_of_entries( 2 ) ;
					setter.set_order_type(
						VariantDataReader::PerSampleSetter::ePerUnorderedHaplotype,
						VariantDataReader::PerSampleSetter::eAlleleIndex
					) ;
					std::pair< int64_t, int64_t > const& genotype = m_source.m_genotype_table[ data ] ;
					if( genotype.first == -1 ) {
						setter( genfile::MissingValue() ) ;
						setter( genfile::MissingValue() ) ;
					} else {
						setter( genotype.first ) ;
						setter( genotype.second ) ;
					}
				}
				return *this ;
			}
			
			std::size_t get_number_of_samples() const { return m_number_of_samples ; }
			
			bool supports( std::string const& spec ) const {
				return spec == "GT" || spec == ":genotypes:";
			}

			void get_supported_specs( SpecSetter setter ) const {
				setter( "GT", "Integer" ) ;
				setter( ":genotypes:", "Integer" ) ;
			}

		private:
			BedFileSNPDataSource& m_source ;
			std::size_t m_number_of_samples ;
			std::vector< char > m_buffer ;
		} ;
	}

	VariantDataReader::UniquePtr BedFileSNPDataSource::read_variant_data_impl() {
		return VariantDataReader::UniquePtr( new impl::BedFileSNPDataReader( *this )) ;
	}

	void BedFileSNPDataSource::ignore_snp_probability_data_impl() {
		m_bed_stream_ptr->ignore( ( m_number_of_samples + 3 ) / 4 ) ;
	}

}	

