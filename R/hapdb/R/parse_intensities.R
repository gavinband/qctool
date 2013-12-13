parse_intensities <-
function( data, N ) {
    if( data[1] == '73' && data[2] == 0 && data[3] == 6 && rawToChar( data[4:9] ) == '2float' ) {
        number_of_samples = readBin( data[11:14], what = "integer", endian = "big" )
        stopifnot( number_of_samples == N ) ;
        number_of_entries = readBin( data[16:17], what = "integer", endian = "big", size = 2 )
        stopifnot( number_of_entries == 2 ) ;
        result = parse_floats( data[18:length(data)], N, 2 ) ;
    }
    else {
        stop() ;
    }
    return(result)
}
