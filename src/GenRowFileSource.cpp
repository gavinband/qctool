#include <vector>
#include <string>
#include <utility>
#include "GenRow.hpp"
#include "GenRowFileSource.hpp"
#include "FileUtil.hpp"

// Open a GEN input file, returning the appropriate type of Source object.
std::auto_ptr< ObjectSource< GenRow > > get_genrow_source_from_file( std::string filename ) {
	std::auto_ptr< ObjectSource< GenRow > > genrow_source( new NullObjectSource< GenRow >() ) ;
	if( exists( filename ) && is_regular( filename )) {
		if( determine_file_mode( filename ) & e_BinaryMode ) {
			genrow_source.reset( new SimpleBinaryFileGenRowSource( open_file_for_input( filename ))) ;
		} else {
			genrow_source.reset( new SimpleTextFileGenRowSource( open_file_for_input( filename ))) ;
		}
	}

	return genrow_source ;
}

// Open a GEN input file, returning the appropriate type of Source object.
std::auto_ptr< ObjectSource< GenRow > > get_genrow_source_from_files( std::vector< std::string > filenames ) {
	std::auto_ptr< ObjectSource< GenRow > > genrow_source( new SNPDataSourceGenRowSource( filenames )) ;
	return genrow_source  ;
}

ChainingFileGenRowSource:: ChainingFileGenRowSource( std::vector< std::string > filenames ) {
	for( std::size_t i = 0; i < filenames.size(); ++i ) {
		std::auto_ptr< ObjectSource< GenRow > > source = get_genrow_source_from_file( filenames[i] ) ;
		add_source( source ) ;
	}
}

