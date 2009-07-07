#ifndef __GTOOL_GENROWFILESOURCE_HPP
#define __GTOOL_GENROWFILESOURCE_HPP

#include <string>
#include <vector>
#include <stddef.h>
#include <boost/bind.hpp>

#include "GenRow.hpp"
#include "SimpleFileObjectSource.hpp"
#include "ChainingFileObjectSource.hpp"
#include "FileUtil.hpp"
#include "SNPDataProvider.hpp"


// Open a GEN input file, returning the appropriate type of Source object.
std::auto_ptr< ObjectSource< GenRow > > get_genrow_source_from_file( std::string filename ) ;
std::auto_ptr< ObjectSource< GenRow > > get_genrow_source_from_files( std::vector< std::string > filenames ) ;


// This class is an adapter between SNPDataProvider and ObjectSource< GenRow >.
struct SNPDataProviderGenRowSource
{
	SNPDataProviderGenRowSource( std::vector< std::string > filenames )
		: m_snp_data_provider( SNPDataProvider::create( filenames ))
	{}

	SNPDataProviderGenRowSource& read( GenRow & row ) {
		m_snp_data_provider->read_snp(
			boost::bind< void >( &GenRow::set_number_of_samples, row, _1 ),
			boost::bind< void >( &GenRow::set_SNPID, row, _1 ),
			boost::bind< void >( &GenRow::set_RSID, row, _1 ),
			boost::bind< void >( &GenRow::set_SNP_position, row, _1 ),
			boost::bind< void >( &GenRow::set_allele1, row, _1 ),
			boost::bind< void >( &GenRow::set_allele2, row, _1 ),
			boost::bind< void >( &GenRow::set_genotype_probabilities, row, _1, _2, _3, _4 )
		) ;
		
		return *this ;
	}

	bool fail() const { return false ; } // TODO: implement this.

	operator bool() {
		return (*m_snp_data_provider) ;
	}

private:

	std::auto_ptr< SNPDataProvider > m_snp_data_provider ;
} ;


typedef SimpleFileObjectSource< GenRow > SimpleTextFileGenRowSource ;

struct SimpleBinaryFileGenRowSource: public SimpleFileObjectSource< GenRow >
{
	typedef SimpleFileObjectSource< GenRow > base_t ;
	
	SimpleBinaryFileGenRowSource( INPUT_FILE_PTR a_stream_ptr )
		: base_t( a_stream_ptr )
	{
		gen::bgen::read_offset( *stream_ptr(), &m_offset ) ;
		stream_ptr()->ignore( m_offset ) ;
		if( !(*stream_ptr())) {
			throw BadFileFormatException( "SimpleBinaryFileGenRowSource: unable to read (or to skip) offset - this can't be a valid binary gen file." ) ;
		}
	}

	SimpleBinaryFileGenRowSource& read( GenRow & row ) {
		row.read_from_binary_stream( *stream_ptr() ) ;
		return *this ;
	}
	
	bool fail() const { return stream_ptr()->fail() ; }
	
	private:
		uint32_t m_offset ;
} ;

struct ChainingFileGenRowSource: public ChainingFileObjectSource< GenRow >
{
	ChainingFileGenRowSource( std::vector< std::string > filenames ) ;
} ;

#endif
