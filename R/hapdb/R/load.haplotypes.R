load.haplotypes <-
function( hapdb, chromosome, rsid = NULL, range = NULL, analysis = NULL ) {
    require( RSQLite )
    require( Rcompression )
    sql = paste(
        "SELECT analysis_id, E.name AS analysis, V.id AS variant_id, chromosome, position, rsid, alleleA, alleleB, H.N AS N, H.data",
        "FROM Haplotype H",
        "INNER JOIN Entity E ON E.id == H.analysis_id",
        "INNER JOIN Variant V ON V.id == H.variant_id",
        "WHERE 1 == 1",
        sep = " "
    )
    
    if( length( unique( hapdb$samples$analysis_id ) ) > 0 ) {
        if( is.null( analysis )) {
            stop( "This hapdb file has more than one analysis, please specify one." ) ;
        } else {
            analysis_id = hapdb$samples$analysis_id[ match( analysis, hapdb$samples$analysis ) ]
            if( is.na( analysis_id )) {
                stop( paste( "Analysis \"", analysis, "\" was not found in the hapdb file.", sep = "" ) ) ;
            }
            sql = paste(
                sql,
                sprintf( "AND analysis_id == %d", analysis_id ),
                sep = " "
            ) ;
        }
    }
    
    if( !is.null( chromosome ) ) {
        sql = paste(
            sql,
            sprintf( "AND chromosome == '%s'", as.character( chromosome ) ),
            sep = " "
        )
    }
    if( !is.null( range ) ) {
        sql = paste(
            sql,
            sprintf( "AND position BETWEEN %d AND %d", as.integer( range[1] ), as.integer( range[2] ) ),
            sep = " "
        )
    }
    if( !is.null( rsid ) ) {
        sql = paste(
            sql,
            "AND rsid IN (",
            paste( sprintf( '"%s"', rsid ), collapse = ',' ),
            sep = " "
        )
    }

    D = dbGetQuery( db, sql )

    result = list(
        variant = D[,-which( colnames(D) == "data")],
        data = matrix( NA, nrow = nrow(D), ncol = 2 * D$N[1] )
    )

    for( i in 1:nrow(D) ) {
        compressed_data = unlist( D$data[i] )
        uncompressed_data = uncompress( compressed_data, asText = FALSE )
        result$data[i,] = parse_haplotypes( uncompressed_data, D$N[i] )
    }
    return( result )
}
