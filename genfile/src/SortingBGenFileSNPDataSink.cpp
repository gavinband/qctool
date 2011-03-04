#include <iostream>
#include <string>
#include <utility>
#include "genfile/snp_data_utils.hpp"
#include "genfile/SNPDataSink.hpp"
#include "genfile/bgen.hpp"
#include "genfile/BGenFileSNPDataSink.hpp"
#include "genfile/GenomePosition.hpp"
#include "genfile/SortingBGenFileSNPDataSink.hpp"

namespace genfile {
	SortingBGenFileSNPDataSink::SortingBGenFileSNPDataSink(
		std::string const& filename,
		std::string const& free_data,
		bgen::uint32_t flags
	)
	: 	BasicBGenFileSNPDataSink( filename, free_data, e_NoCompression, flags )
	{
		assert( flags & bgen::e_CompressedSNPBlocks ) ;
		// Reserve enough space for lots of data!
		m_rows.reserve( 10000 ) ;
	}

	void SortingBGenFileSNPDataSink::write_snp_impl(
		uint32_t number_of_samples,
		std::string SNPID,
		std::string RSID,
		Chromosome chromosome,
		uint32_t SNP_position,
		char first_allele,
		char second_allele,
		GenotypeProbabilityGetter const& get_AA_probability,
		GenotypeProbabilityGetter const& get_AB_probability,
		GenotypeProbabilityGetter const& get_BB_probability
	) {
		// Add a row to our data store...
		m_rows.push_back( StoredRow( SNPKey( GenomePosition( chromosome, SNP_position ), RSID ), std::string() ) ) ;

		std::size_t id_field_size = std::min( std::max( SNPID.size(), RSID.size() ), static_cast< std::size_t >( 255 )) ;
		std::ostringstream ostr( std::ios::binary ) ;
		bgen::write_compressed_snp_block( ostr, number_of_samples, id_field_size, SNPID, RSID, chromosome, SNP_position, first_allele, second_allele, get_AA_probability, get_AB_probability, get_BB_probability ) ;				

		m_rows.back().second = ostr.str() ;
	}

	SortingBGenFileSNPDataSink::~SortingBGenFileSNPDataSink() {
		// We are about to close the file.
		// First sort and write all our rows.
		std::sort( m_rows.begin(), m_rows.end() ) ;
		// Now write them.
		for( std::size_t i = 0; i < m_rows.size(); ++i ) {
			stream_ptr()->write( m_rows[i].second.data(), m_rows[i].second.size() ) ;
		}

		// Now we must seek back to the beginning and write the header block.
		stream_ptr()->seekp(4, std::ios_base::beg ) ;
		if( stream_ptr()->bad() ) {
			throw OutputError( filename() ) ;
		}

		write_header_data( *stream_ptr() ) ;
	}
}
