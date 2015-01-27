parse_haplotypes <-
function( data, N ) {
    if( data[1] == '73' && data[2] == 0 && data[3] == 7 && rawToChar( data[4:10] ) == 'bitpack' ) {
        result = parse_bitpack_haplotypes( data[11:length(data)], N ) ;
    } else if( data[1] == 's' && data[2] == 0 && data[3] == 6 && rawToChar( data[4:9] ) == 'ubjson' ) {
        result = parse_ubjson_haplotypes( data[10:length(data)], N ) ;
    } else {
        stop() ;
    }
    return(result)
}
