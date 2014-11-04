
//          Copyright Gavin Band 2008 - 2012.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <fstream>
#include <cstdio>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/bgen/bgen.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"
#include "genfile/GenomePosition.hpp"
#include "genfile/SortingBGenFileSNPDataSink.hpp"

namespace genfile {
	SortingBGenFileSNPDataSink::SortingBGenFileSNPDataSink(
		std::string const& filename, 
		SNPDataSink::UniquePtr sink
	):
		m_filename( filename ),
		m_sink( sink ),
		m_offset_of_first_snp( 0 )
	{
	}

	std::string SortingBGenFileSNPDataSink::get_spec() const {
		return m_filename ;
	}


	void SortingBGenFileSNPDataSink::write_snp_impl(
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
		SNPIdentifyingData snp(
			SNPID,
			RSID,
			GenomePosition( chromosome, SNP_position ),
			first_allele,
			second_allele
		) ;
		OffsetMap::iterator offset_i = m_file_offsets.insert(
			std::make_pair(
				snp,
				std::make_pair(
					m_sink->get_stream_pos().second,
					m_sink->get_stream_pos().second
				)
			)
		) ;

		m_sink->write_snp(
			number_of_samples,
			SNPID,
			RSID,
			chromosome,
			SNP_position,
			first_allele,
			second_allele,
			get_AA_probability,
			get_AB_probability,
			get_BB_probability,
			info
		) ;
		
		offset_i->second.second = m_sink->get_stream_pos().second ;
	}

	SortingBGenFileSNPDataSink::~SortingBGenFileSNPDataSink() {
		// Ensure temporary file is flushed.
		m_sink.reset() ;
		
		boost::system::error_code ec ;

		// Move file to a similarly-named temporary.
#if BOOST_FILESYSTEM_VERSION < 3
		boost::filesystem::path temp_filename = boost::filesystem::unique_path( m_filename + ".tmp%%%%-%%%%-%%%%-%%%%", ec ) ;
#else
		boost::filesystem::path temp_filename = boost::filesystem::unique_path( m_filename + ".tmp%%%%-%%%%-%%%%-%%%%", ec ) ;
#endif
		// std::cerr << "Renaming \"" << m_filename << "\" to \"" << temp_filename << "\"...\n" ;
		boost::filesystem::rename(
			boost::filesystem::path( m_filename ),
			boost::filesystem::path( temp_filename )
		) ;
		// std::cerr << "Copying file back...\n" ;
		// Copy temporary back to file in the right order.
		std::ifstream input( temp_filename.c_str(), std::ios::binary ) ;
		std::ofstream output( m_filename.c_str(), std::ios::binary | std::ios::trunc ) ;

		std::vector< char > buffer( 1024*1024 ) ;

		if( m_file_offsets.empty() ) {
			while( input.read( &buffer[0], buffer.size() ) ) {
				output.write( &buffer[0], input.gcount() ) ;
			}
		}
		else {
			OffsetMap::const_iterator
				i = m_file_offsets.begin(),
				end_i = m_file_offsets.end()
			;

			// copy header.
			{
				std::pair< std::ostream::streampos, std::ostream::streampos > chunk ;
				chunk.first = 0 ;
				chunk.second = m_offset_of_first_snp ;
				// std::cerr << "copying 0th chunk " << chunk.first << " - " << chunk.second << " of " << m_file_offsets.size() << "...\n" ;
				for( std::ostream::streampos i = 0; i < chunk.second; i += buffer.size() ) {
					std::size_t n = std::min( std::size_t( chunk.second - i ), buffer.size() ) ;
					input.read( &buffer[0], n ) ;
					output.write( &buffer[0], n ) ;
				}
			}
			for( ; i != end_i; ++i ) {
				std::pair< std::ostream::streampos, std::ostream::streampos > const& chunk = i->second ;
				input.seekg( chunk.first ) ;
				// std::cerr << "copying chunk " << chunk.first << " - " << chunk.second << "...\n" ;
				for( std::ostream::streampos i = chunk.first; i < chunk.second; i += buffer.size() ) {
					std::size_t n = std::min( std::size_t( chunk.second - i ), buffer.size() ) ;
					input.read( &buffer[0], n ) ;
					output.write( &buffer[0], n ) ;
				}
			}
		}
		output.close() ;
		input.close() ;

		// remove the temporary file.
		boost::filesystem::remove( temp_filename ) ;
		// ignore the error code.
	}
	
	void SortingBGenFileSNPDataSink::set_sample_names_impl( std::size_t number_of_samples, SampleNameGetter name_getter ) {
		m_sink->set_sample_names( number_of_samples, name_getter ) ;
		m_offset_of_first_snp = m_sink->get_stream_pos().second ;
	}
}
