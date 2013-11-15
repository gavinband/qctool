load.haplotypes <-
function( db, rsid, analysis = NULL ) {
    if( is.null( analysis ) ) {
        sql = sprintf( "SELECT analysis_id, E.name AS analysis, V.id AS variant_id, chromosome, position, rsid, alleleA, alleleB, H.N AS N, H.data FROM Haplotype H INNER JOIN Entity E ON E.id == H.analysis_id INNER JOIN Variant V ON V.id == H.variant_id AND V.rsid == '%s'", rsid )
    } else {
        sql = sprintf( "SELECT analysis_id, E.name AS analysis, V.id AS variant_id, chromosome, position, rsid, alleleA, alleleB, H.N AS N, H.data FROM Haplotype H INNER JOIN Entity E ON E.id == H.analysis_id INNER JOIN Variant V ON V.id == H.variant_id AND E.name == '%s' AND V.rsid == '%s'", analysis, rsid )
    }
    require( Rcompression )
    D = dbGetQuery( db, sql )

    result = list(
        variant = D[,-which( colnames(D) == "data")],
        data = matrix( NA, nrow = nrow(D), ncol = 2 * D$N[1] )
    )

    for( i in 1:nrow(D) ) {
        compressed_data = unlist( D$data[i] )
        uncompressed_data = uncompress( compressed_data, asText = FALSE )
        result$data[i,] = parse_haplotypes( uncompressed_data[10:length(uncompressed_data)], D$N[i] ) ;
    }
    return( result )
}
