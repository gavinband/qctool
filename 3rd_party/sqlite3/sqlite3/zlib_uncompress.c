#include <stdlib.h>
#include <sqlite3/sqlite3ext.h>
#include <zlib.h>
#include <assert.h>
#include <stdint.h>
SQLITE_EXTENSION_INIT1

// decompress data compressed with genfile::zlib_compress()
// This sqlite3 function takes two arguments, the uncompressed size and the compressed data.
static void uncompressFunc( sqlite3_context *context, int argc, sqlite3_value **argv){
	assert( argc==2 );

	sqlite3_int64 uncompressed_size = sqlite3_value_int64( argv[0] ) ;
	sqlite3_int64 nIn = sqlite3_value_bytes( argv[1] ) ;

	if( nIn == 0 ) {
		return;
	}

	unsigned char const* inBuf = sqlite3_value_blob( argv[1] ) ;
	unsigned char* outBuf = malloc( uncompressed_size ) ;
	uLongf out_size = uncompressed_size ;
	int rc = uncompress( outBuf, &out_size, &inBuf[0], nIn ) ;
	if( rc != Z_OK ) {
		free( outBuf ) ;
	} else {
		sqlite3_result_blob( context, outBuf, out_size, free );
	}
}


int sqlite3_extension_init(
	sqlite3 *db,
	char **pzErrMsg,
	const sqlite3_api_routines *pApi
){
  SQLITE_EXTENSION_INIT2( pApi )
  sqlite3_create_function( db, "zlib_uncompress", 2, SQLITE_UTF8, 0, uncompressFunc, 0, 0 );
  return 0;
}
