load.all.haplotypes <-
function( db, analysis = NULL ) {
    if( is.null( analysis ) ) {
        sql = sprintf( "SELECT analysis_id, E.name AS analysis, V.id AS variant_id, chromosome, position, rsid, alleleA, alleleB, H.N AS N, H.data FROM Haplotype H INNER JOIN Entity E ON E.id == H.analysis_id INNER JOIN Variant V ON V.id == H.variant_id" )
    } else {
        sql = sprintf( "SELECT analysis_id, E.name AS analysis, V.id AS variant_id, chromosome, position, rsid, alleleA, alleleB, H.N AS N, H.data FROM Haplotype H INNER JOIN Entity E ON E.id == H.analysis_id INNER JOIN Variant V ON V.id == H.variant_id AND E.name == '%s'", analysis )
    }
    require( Rcompression )
    D = dbGetQuery( db, sql )

    str(D)
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
