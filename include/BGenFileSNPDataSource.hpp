#ifndef BGENFILESNPDATAPROVIDER_HPP
#define BGENFILESNPDATAPROVIDER_HPP

#include <iostream>
#include <string>
#include "snp_data_utils.hpp"
#include "SNPDataSource.hpp"
#include "bgen.hpp"

// This class represents a SNPDataSource which reads its data
// from a BGEN file.
class BGenFileSNPDataSource: public SNPDataSource
{
public:
	BGenFileSNPDataSource( std::string const& filename )
		: m_filename( filename )
	{
		setup( filename, impl::file_is_gzipped( filename )) ;
	}

	BGenFileSNPDataSource( std::string const& filename, bool file_is_gzipped )
		: m_filename( filename )
	{
		setup( filename, file_is_gzipped ) ;
	}

	unsigned int number_of_samples() const { return m_number_of_samples ; }
	unsigned int total_number_of_snps() const { return m_total_number_of_snps ; }
	FormatType format() const { return e_BGenFormat ; }
	std::istream& stream() { return *m_stream_ptr ; }
	std::istream const& stream() const { return *m_stream_ptr ; }

private:

	std::string m_filename ;
	unsigned int m_number_of_samples, m_total_number_of_snps ;
	std::auto_ptr< std::istream > m_stream_ptr ;

	void setup( std::string const& filename, bool file_is_gzipped ) {
		gen::bgen::uint32_t offset ;
		gen::bgen::read_offset( (*m_stream_ptr), &offset ) ;
		
		m_stream_ptr = impl::open_bgen_file( filename, file_is_gzipped ) ;
		gen::bgen::uint32_t header_size = read_header_data() ;

		if( offset < header_size ) {
			throw FileStructureInvalidError() ;
		}

		// skip any remaining bytes before the first snp block
		m_stream_ptr->ignore( offset - header_size ) ;			
	}
	
	gen::bgen::uint32_t read_header_data() {
		gen::bgen::uint32_t header_size ;
		
		gen::bgen::read_header_block(
			(*m_stream_ptr),
			gen::set_value( header_size ),
			gen::set_value( m_total_number_of_snps ),
			gen::set_value( m_number_of_samples ),
			gen::ignore(),
			gen::ignore(),
			gen::ignore()
		) ;

		return header_size ;
	}
} ;
	


#endif