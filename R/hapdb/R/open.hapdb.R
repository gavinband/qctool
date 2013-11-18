open.hapdb <-
function( filename, get.variants = FALSE ) {
    db = dbConnect( dbDriver( "SQLite" ), filename ) ;
    cat( "open.hapdb(): loading samples and variants..." ) ;
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
	}
    return( result ) ;
}
