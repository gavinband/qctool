open.hapdb <-
function( filename, get.variants = FALSE ) {
    require( RSQLite )
    db = dbConnect( dbDriver( "SQLite" ), filename ) ;
    samples = dbGetQuery(
        db,
        paste(
            'SELECT S.analysis_id AS analysis_id, E.name AS analysis, identifier, index_in_data AS "index"',
            'FROM Sample S',
            'INNER JOIN Entity E ON E.id == S.analysis_id',
            'ORDER BY analysis_id, index_in_data',
            sep = " "
        )
    )
	result = list(
		db = db,
		samples = samples
	) ;
	if( get.variants ) {
		result$variants = get.variants( result ) ;
    	cat( sprintf( "open.hapdb(): opened hapdb with %d samples and %d variants.\n", nrow( result$samples ), nrow( result$variants ) ) )
	} else {
    	cat( sprintf( "open.hapdb(): opened hapdb with %d samples.\n", nrow( result$samples ) ) )
	}
    return( result ) ;
}
