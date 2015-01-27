bgen.read_header <-
function( con, endian = "little" ) {
if( class( con ) == "character" ) { con = file( con, open = "rb" ) }
    data = as.list( readBin( con, integer(), n = 5, size = 4, endian = endian ) )
    names( data ) = c( "offset", "header_length", "number_of_snps", "number_of_samples", "reserved" )
    if( data$header_length > 20 ) {
        data[[ 'free_data' ]] = readChar( con, nchars = ( data$header_length - 20 ) )
    }
    else {
        data[[ 'free_data' ]] = ""
    }
    data[[ 'flags' ]] = readBin( con, integer(), n = 1, size = 4, endian = endian )
    data[[ 'endian' ]] = endian
    data[[ 'type' ]] = 'bgen'
    data$connection = con
    return( data ) ;
}
