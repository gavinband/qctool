parse_floats <-
function( data, N, number_of_elements ) {
    stopifnot( length( data ) == N*4*number_of_elements ) ;
    result = readBin( data, what = "numeric", size = 4, n = N*number_of_elements )
    result[ which( result == -1 )] = NA
    return( result )
}
