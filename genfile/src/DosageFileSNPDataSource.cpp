
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include "genfile/snp_data_utils.hpp"
#include "genfile/gen.hpp"
#include "genfile/SNPDataSource.hpp"
#include "genfile/DosageFileSNPDataSource.hpp"

namespace genfile {
	DosageFileSNPDataSource::DosageFileSNPDataSource( std::string const& filename, Chromosome chromosome )
		: m_filename( filename ),
		  m_compression_type( get_compression_type_indicated_by_filename( filename ) ),
		  m_number_of_samples( 0 ),
		  m_chromosome( chromosome )
	{
		setup( filename, m_compression_type ) ; 
	}

	DosageFileSNPDataSource::DosageFileSNPDataSource(
		std::string const& filename,
		Chromosome chromosome,
		CompressionType compression_type
	)
		: m_filename( filename ),
		  m_compression_type( compression_type ),
		  m_number_of_samples( 0 ),
		  m_chromosome( chromosome )
	{
		setup( filename, compression_type ) ;
	}

	void DosageFileSNPDataSource::setup( std::string const& filename, CompressionType compression_type ) {
		setup( open_text_file_for_input( filename, compression_type ) ) ;
	}

	void DosageFileSNPDataSource::setup( std::auto_ptr< std::istream > stream_ptr ) {
		m_stream_ptr = stream_ptr ;
        std::getline( stream(), m_line ) ;
        bool bad = false ;
        if( !stream() ) {
            bad = true ;
        } else {
            slice( m_line ).split( " ", &m_elts ) ;
            if( m_elts.size() < 6 ) {
                bad = true ;
            } else if(
                m_elts[0] != "chromosome"
                || m_elts[1] != "SNPID"
                || m_elts[2] != "rsid"
                || m_elts[3] != "position"
                || m_elts[4] != "alleleA"
                || m_elts[5] != "alleleB"
            ) {
                bad = true ;
            }
        }
        
        if( bad ) {
            throw genfile::MalformedInputError(
                m_filename,
                "Malformed header line",
                0
            ) ;
        } else {
            m_number_of_samples = m_elts.size() - 6 ;
        }
	}

	void DosageFileSNPDataSource::reset_to_start_impl() {
		stream().clear() ;
		stream().seekg( 0 ) ;
		if( !stream() ) {
			// oh, dear, seeking failed.  Reopen the file instead
			m_stream_ptr = open_text_file_for_input( m_filename, m_compression_type ) ;
		}
		if( !stream() ) {
			throw OperationFailedError( "genfile::DosageFileSNPDataSource::reset_to_start_impl()", get_source_spec(), "reset to start" ) ;
		}
        std::getline( stream(), m_line ) ;
	}
	
	SNPDataSource::Metadata DosageFileSNPDataSource::get_metadata() const {
		std::map< std::string, std::string > format ;
		format[ "ID" ] = "Dosage" ;
		format[ "Number" ] = "1" ;
		format[ "Description" ] = "Genotype dosage" ;
		SNPDataSource::Metadata result ;
		result.insert( std::make_pair( "FORMAT", format )) ;
		return result ;
	}
	
	void DosageFileSNPDataSource::read_snp_identifying_data_impl( 
		uint32_t* number_of_samples,
		std::string* SNPID,
		std::string* RSID,
		Chromosome* chromosome,
		uint32_t* SNP_position,
		std::string* allele1,
		std::string* allele2
	) {
        std::string SNPID_, RSID_, allele1_, allele2_, chromosome_string ;
        uint32_t position_ ;
        stream() >> chromosome_string ;
        if( !stream() && !m_stream_ptr->eof() ) {
            throw genfile::MalformedInputError(
                m_filename,
                "Unable to read chromosome, but not at EOF",
                number_of_snps_read()
            ) ;
        }
		if( stream() ) {
			stream() >> SNPID_ >> RSID_ >> position_ >> allele1_ >> allele2_ ;
			if( !stream() ) {
				throw genfile::MalformedInputError(
					m_filename,
					"Malformed line",
					number_of_snps_read()
				) ;
			} else {
				*SNPID = SNPID_ ;
				*RSID = RSID_ ;
				*chromosome = Chromosome( chromosome_string ) ;
				*SNP_position = position_ ;
				*allele1 = allele1_ ;
				*allele2 = allele2_ ;
				// Prepare for the data...
				if( stream().peek() == ' ' ) {
					stream().get() ;
				}
			}
		}
	}

	namespace impl {
		struct DosageFileSNPDataReader: public VariantDataReader {
			DosageFileSNPDataReader( DosageFileSNPDataSource& source, std::vector< string_utils::slice > const& elts ):
                m_elts( elts )
            {
				//std::cerr << "m_elts.size() == " << m_elts.size() << ", source.number_of_samples() == " << source.number_of_samples() << ".\n" ;
				//std::cerr << "First few elts are: " << m_elts[0] << ", " << m_elts[1] << ", " << m_elts[2] << ",...\n" ;
                assert( elts.size() == source.number_of_samples() ) ;
            }
			
			DosageFileSNPDataReader& get( std::string const& spec, PerSampleSetter& setter ) {
				std::size_t const N = m_elts.size() ;
				setter.set_number_of_samples( N ) ;
				setter.set_number_of_alleles( 2 ) ;
				for( std::size_t i = 0; i < N; ++i ) {
					setter.set_sample( i ) ;
					setter.set_number_of_entries( 1 ) ;
					setter.set_order_type( PerSampleSetter::eBAlleleDosage, PerSampleSetter::eDosage ) ;
					if( m_elts[i] == "NA" ) {
						setter( genfile::MissingValue() ) ;
					} else {
						setter( string_utils::to_repr< Integer >( m_elts[i] )) ;
					}
				}
				return *this ;
			}
			
			std::size_t get_number_of_samples() const {
				return m_elts.size() ; 
			}
			
			bool supports( std::string const& spec ) const {
				return spec == ":genotypes:" ;
			}
			
			void get_supported_specs( SpecSetter setter ) const {
				setter( ":genotypes:", "Genotype" ) ;
			}
			
		private:
            std::vector< string_utils::slice > m_elts ;
			std::vector< double > m_genotypes ;
		} ;
	}

	VariantDataReader::UniquePtr DosageFileSNPDataSource::read_variant_data_impl() {
        m_line.clear() ;
		m_elts.clear() ;
		std::getline( *m_stream_ptr, m_line ) ;
		slice( m_line ).split( " ", &m_elts ) ;
		return VariantDataReader::UniquePtr( new impl::DosageFileSNPDataReader( *this, m_elts ) ) ;
	}

	void DosageFileSNPDataSource::ignore_snp_probability_data_impl() {
        // nothing to do ;
	}
}

