open.hapdb <-
function( filename, get.variants = FALSE ) {
    db = dbConnect( dbDriver( "SQLite" ), filename ) ;
    cat( "open.hapdb(): loading samples and variants..." ) ;
    samples = dbGetQuery( db, "SELECT * FROM Sample ORDER BY analysis_id, index_in_data" ) ;
	result = list(
		db = db,
		samples = samples
	) ;
	if( get.variants ) {
		result$variants = dbGetQuery( db, "SELECT * FROM Variant" ) ;
	}
    return( result ) ;
}
