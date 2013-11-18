parse_bitpack_haplotypes <-
function( haps, N ) {
    stopifnot( length( haps ) == 2*N ) ;
    return( as.integer( haps ))
}
