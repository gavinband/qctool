parse_ubjson_haplotypes <-
function( db, N ) {
    n = 0 ;
    stopifnot( haps[1] == '5b' )
    i = 2 ;
    v = as.integer(
        haps[ sort( union( seq( from = 4, by = 6, length = N ), seq( from = 6, by = 6, length = N ) ) ) ]
    )
    return(v)
}
