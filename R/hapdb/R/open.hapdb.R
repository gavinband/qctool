open.hapdb <-
function( filename, get.variants = FALSE ) {
    require( RSQLite )
    db = dbConnect( dbDriver( "SQLite" ), filename ) ;
	dbGetQuery( db, "PRAGMA omit_readlock=on" )
    samples = dbGetQuery(
        db,
        paste(
            'SELECT E.name AS analysis, S.*',
            'FROM Sample S',
            'INNER JOIN Analysis E ON E.id == S.analysis_id',
            'ORDER BY analysis_id, index_in_data',
            sep = " "
        )
    )
    wdup = which( duplicated( colnames( samples ) ) )
    if( length( wdup ) > 0 ) {
        samples = samples[,-wdup]
    }
    result = list(
        db = db,
        samples = samples
    ) ;
    if( get.variants ) {
        result$variants = dbGetQuery( db, "SELECT * FROM Variant" ) ;
        cat( sprintf( "open.hapdb(): opened hapdb with %d samples and %d variants.\n", nrow( result$samples ), nrow( result$variants ) ) )
    } else {
        cat( sprintf( "open.hapdb(): opened hapdb with %d samples.\n", nrow( result$samples ) ) )
    }
    return( result ) ;
}
