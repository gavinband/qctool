load.sample.file <- function( filename, na.strings = "NA", sep = " " ) {
	# Read header and types
	samples = read.table( filename, as.is = T, hea=T, sep = sep )
	sample.header = colnames( samples )
	types = as.character( samples[1,] )
	# Now read the data rows
	samples = read.delim( filename, as.is = T, hea=F, skip = 2, sep = sep, na.strings = na.strings )
	names( samples ) = sample.header
	
	# In case IDs happened to be numerical, make the first two columns character vectors
	samples[,1] = as.character( samples[,1] )
	samples[,2] = as.character( samples[,2] )

	# Convert discrete levels to factors.
	for( i in 1:ncol( samples )) {
		if( types[i] == 'D' | types[i] == 'B' ) {
			samples[,i] = as.factor( samples[,i] )
		}
	}
	
	return( samples ) ;
}

