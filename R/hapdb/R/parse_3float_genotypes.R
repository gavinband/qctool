parse_3float_genotypes <-
function( data, N ) {
    stopifnot( length( data ) == N*3*4 ) ;
    result = readBin( data, what = "numeric", size = 4, n = N*3 )
    result[ which( result == -1 )] = NA
    return( result )
}
