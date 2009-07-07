#include <iostream>
#include <string>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include "snp_data_utils.hpp"
#include "SNPDataProvider.hpp"
#include "GenFileSNPDataProvider.hpp"
#include "BGenFileSNPDataProvider.hpp"
#include "SNPDataProviderChain.hpp"

namespace impl {
	bool gen_file_is_in_binary_format( std::string filename ) {
		return ( filename.find( ".bgen") != std::string::npos ) ;
	}

	bool file_is_gzipped( std::string filename ) {
		return ( filename.find( ".gz") != std::string::npos ) ;
	}
	
	std::auto_ptr< std::istream > open_gen_file( std::string filename, bool file_is_gzipped ) {
		std::auto_ptr< boost::iostreams::filtering_istream > gen_file_ptr( new boost::iostreams::filtering_istream ) ;
	    if (file_is_gzipped) gen_file_ptr->push(boost::iostreams::gzip_decompressor());
		gen_file_ptr->push(boost::iostreams::file_source(filename.c_str())); 
		return std::auto_ptr< std::istream >( gen_file_ptr ) ;
  	}

	std::auto_ptr< std::istream > open_bgen_file( std::string filename, bool file_is_gzipped ) {
		std::auto_ptr< boost::iostreams::filtering_istream > gen_file_ptr( new boost::iostreams::filtering_istream ) ;
	    if (file_is_gzipped) gen_file_ptr->push(boost::iostreams::gzip_decompressor());
		gen_file_ptr->push(boost::iostreams::file_source(filename.c_str()), std::ios_base::binary ); 
		return std::auto_ptr< std::istream >( gen_file_ptr ) ;
  	}
}

namespace gen {
	Ignorer ignore() { return Ignorer() ; }
}

std::auto_ptr< SNPDataProvider > SNPDataProvider::create( std::string const& filename ) {
	return SNPDataProvider::create( filename, impl::file_is_gzipped( filename )) ;
}

std::auto_ptr< SNPDataProvider > SNPDataProvider::create( std::string const& filename, bool file_is_gzipped ) {
	if( impl::gen_file_is_in_binary_format( filename )) {
		return std::auto_ptr< SNPDataProvider >( new BGenFileSNPDataProvider( filename, file_is_gzipped )) ;
	}
	else {
		return std::auto_ptr< SNPDataProvider >( new GenFileSNPDataProvider( filename, file_is_gzipped )) ;
	}
}

std::auto_ptr< SNPDataProvider > SNPDataProvider::create( std::vector< std::string > const& filenames ) {
	std::auto_ptr< SNPDataProviderChain > chain( new SNPDataProviderChain() ) ;
	for( std::size_t i = 0; i < filenames.size(); ++i ) {
		chain->add_provider( SNPDataProvider::create( filenames[i] )) ;
	}
	return std::auto_ptr< SNPDataProvider >( chain.release() ) ;
}
