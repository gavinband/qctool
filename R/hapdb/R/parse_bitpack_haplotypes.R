parse_bitpack_haplotypes <-
function( haps, N ) {
   # stopifnot( length( haps ) == 2*N ) ;
   print( haps )
	result = as.integer( haps )
	result[ which( result == 255 ) ] = NA
    return( result )
}
