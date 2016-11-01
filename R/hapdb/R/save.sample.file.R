save.sample.file <- function( data, filename ) {
	con = file( filename, open = "w" )
	write( colnames( data ), con, ncol = ncol( data ))
	write( attr( data, 'types' ), con, append = T, ncol = ncol( data ) )
	write.table( data, con, append = T, col.names = F, row.names = F, quote = F )
}
