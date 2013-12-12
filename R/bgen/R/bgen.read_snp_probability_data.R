bgen.read_snp_probability_data <-
function( header, data ) {
    con = header$connection
    endian = header$endian
    if( header$flags == 0 ) {
        data$genotypes = readBin( con, integer(), size = 2, n = ( data$number_of_samples * 3 ), endian = endian ) ;
    } else {
        library( Rcompression )
        compressed_size = readBin( con, integer(), size = 4, n = 1, endian = endian )
        compressed_genotypes = readBin( con, 'raw', n = compressed_size, endian = endian )
        uncompressed_genotypes = uncompress( content=compressed_genotypes, size = data$number_of_samples * 6, asText = FALSE )
        data$genotypes = readBin( uncompressed_genotypes, integer(), size = 2, n = data$number_of_samples * 3, signed = TRUE, endian = endian )
        data$genotypes = matrix(
            data = data$genotypes / 10000.0,
            ncol = 3,
            nrow = data$number_of_samples,
            byrow = T
        )
    }
    data$ploidy = 2 ;
    return( data ) ;
}
