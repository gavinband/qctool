
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <vector>
#include <utility>
#include <fstream>
#include <boost/format.hpp>
#include "genfile/SNPDataSink.hpp"
#include "genfile/CohortIndividualSource.hpp"
#include "genfile/BedFileSNPDataSink.hpp"
#include "genfile/snp_data_utils.hpp"
#include "genfile/string_utils.hpp"

#define DEBUG_BED_FILE_SNP_DATA_SINK 1

namespace genfile {
	BedFileSNPDataSink::BedFileSNPDataSink(
		std::string const& output_bed_filename,
		double call_threshhold
	):
		m_call_threshhold( call_threshhold )
	{
		assert( call_threshhold > 0.5 ) ;
		if( output_bed_filename.size() < 4 || output_bed_filename.substr( output_bed_filename.size() - 4, 4 ) != ".bed" ) {
			throw BadArgumentError( "BedFileSNPDataSink::BedFileSNPDataSink", "output_bed_filename = \"" + output_bed_filename + "\"" ) ;
		}
		m_output_filename_stub = output_bed_filename.substr( 0, output_bed_filename.size() - 4 ) ;
		setup() ;
	}

	std::string BedFileSNPDataSink::get_spec() const {
		return m_output_filename_stub ;
	}

	
	void BedFileSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter getter ) {
#if DEBUG_BED_FILE_SNP_DATA_SINK
		std::cerr << "BedFileSNPDataSink::set_sample_names().\n" ;
#endif
		m_sample_ids.clear() ;
		m_sample_ids.resize( number_of_samples ) ;
		for( std::size_t i = 0; i < number_of_samples; ++i ) {
			m_sample_ids[i] = getter(i).as< std::string >() ;
		}
	}

	void BedFileSNPDataSink::set_metadata_impl( Metadata const& ) {
		// do nothing.  Metadata not supported.
	}
	
	void BedFileSNPDataSink::setup() {
		m_bed_file = open_binary_file_for_output( m_output_filename_stub + ".bed", "no_compression" ) ;
		m_bim_file = open_text_file_for_output( m_output_filename_stub + ".bim", "no_compression" ) ;
		
		std::vector< char > data( 3, 0 ) ;
		// BED magic number
		data[0] = 0x6C ;
		data[1] = 0x1B ;
		// Mode is SNP-major
		data[2] = 0x1 ;
		m_bed_file->write( &data[0], 3 ) ;
	}
	
	
	BedFileSNPDataSink::~BedFileSNPDataSink() {
	}

	void BedFileSNPDataSink::write_snp_impl(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		Chromosome chromosome,
		uint32_t SNP_position,
		std::string first_allele,
		std::string second_allele,
		GenotypeProbabilityGetter const& get_AA_probability,
		GenotypeProbabilityGetter const& get_AB_probability,
		GenotypeProbabilityGetter const& get_BB_probability,
		Info const& info
	) {
		assert( number_of_samples == m_sample_ids.size() ) ;
		/* Write BIM file line. */
		(*m_bim_file ) << chromosome << " " << RSID << " 0 " << SNP_position << " " << first_allele << " " << second_allele << "\n" ;

		/* Write BED file line */
		m_buffer.resize( (m_sample_ids.size()+3) / 4, 0 ) ;
		for( std::size_t group = 0; group < m_sample_ids.size(); group += 4 ) {
			std::size_t c_i = group / 4 ;
			m_buffer[c_i] = 0 ;
			for( std::size_t i = group; i < std::min( group + 4, m_sample_ids.size() ); ++i ) {
				double
					AA = get_AA_probability( i ),
					AB = get_AB_probability( i ),
					BB = get_BB_probability( i ) ;
				char d = 0 ;
				if( AA > m_call_threshhold ) {
					d = 0x0 ;
				}
				else if( AB > m_call_threshhold ) {
					d = 0x1 ;
				}
				else if( BB > m_call_threshhold ) {
					d = 0x3 ;
				}
				else {
					d = 0x2 ;
				}
				d <<= ( 2 * ( i % 4 )) ;
				//std::cerr << "sample " << i << "; " << ( boost::format( "AA=%f AB = %f BB = %f\n") % AA % AB % BB ) ;
				//std::cerr << "d = " << std::size_t(d) << "\n" ;
				m_buffer[c_i] |= d ;
			}
		}

		m_bed_file->write( &m_buffer[0], m_buffer.size() ) ;
	}
}
