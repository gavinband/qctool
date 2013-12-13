parse_3float_genotypes <-
function( data, N ) {
    stopifnot( length( data ) == N*3*4 ) ;
    result = readBin( data, what = "numeric", size = 4 )
    return( result )
}
