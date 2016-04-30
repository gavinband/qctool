library( argparse )
parser <- ArgumentParser( description = 'Generate and test randomly generated GEN files through bgen' )
parser$add_argument(
	"--iterations",
	type = "integer",
	nargs = 1,
	help = "Number of iterations",
	default = 100
)
parser$add_argument(
	"--bits",
	type = "integer",
	nargs = 1,
	help = "Number of bits",
	default = 8
)
parser$add_argument(
	"--variants",
	type = "integer",
	nargs = 1,
	help = "Number of variants",
	default = 100
)
parser$add_argument(
	"--max_samples",
	type = "integer",
	nargs = 1,
	help = "Maximum number of samples to simulate",
	default = 100
)
parser$add_argument(
	"--bgen_version",
	type = "character",
	nargs = 1,
	help = "bgen version to use = v11 or v12",
	default = "v12"
)
parser$add_argument(
	"--qctool",
	type = "character",
	nargs = 1,
	help = "Path to qctool executable",
	default = "qctool_v2.0-dev"
)
parser$add_argument(
	"--verbose",
	help = "Say what we're doing",
	default = FALSE,
	action = "store_true"
)

args = parser$parse_args()

if( args$bgen_version == "v11" ) {
	stopifnot( args$bits == 16 )
} else if( args$bgen_version == "v12" ) {
	stopifnot( args$bits <= 32 )
} else {
	stop( sprintf( "Unrecognised bgen version %s.", args$bgen_version ))
}

chromosomes = c( sprintf( "%02d", 1:22 ), "0X", "0Y", "XY", "MT" )

simulate <- function( numberOfSamples, numberOfVariants, filename ) {
	G = matrix( 0, nrow = numberOfVariants, ncol = 3*numberOfSamples )
	G[, seq( from = 1, length = numberOfSamples, by = 3 )] = runif( numberOfSamples * numberOfVariants )
	G[, seq( from = 2, length = numberOfSamples, by = 3 )] = (
		runif( numberOfSamples * numberOfVariants ) * ( 1 - G[, seq( from = 1, length = numberOfSamples, by = 3 )] )
	)
	G[, seq( from = 3, length = numberOfSamples, by = 3 )] = (
		( 1 - G[, seq( from = 1, length = numberOfSamples, by = 3 )] - G[, seq( from = 2, length = numberOfSamples, by = 3 )] )
	)
	SNPID.length = sample( 1:100, 1 )
	rsid.length = sample( 1:100, 1 )
	V = data.frame(
		SNPID = paste( sample( c( LETTERS, letters ), SNPID.length, replace = T ), collapse = '' ),
		rsid = paste( sample( c( LETTERS, letters ), rsid.length, replace = T ), collapse = '' ),
		chromosome = sample( chromosomes, 1 ),
		position = sample( 1:1000000000, 1 ),
		alleleA = paste( sample( c( LETTERS, letters ), sample( 1:10, 1 ) ), collapse = '' ),
		alleleB = paste( sample( c( LETTERS, letters ), sample( 1:10, 1 ) ), collapse = '' )
	)
	
	write.table(
		cbind( V, G ),
		col.names = F,
		row.names = F,
		quote = F,
		file = filename
	)
}

filename=tempfile( tmpdir = "/tmp/" )
bgen.filename = sprintf( "%s.bgen", filename )
gen.filename = sprintf( "%s.bgen.gen", filename )
for( i in 1:args$iterations ) {
	N = sample( 1:args$max_samples, 1 )
	simulate( N, args$variants, filename )
	original.data = read.table( filename, head = F, as.is=T )
	system(
		sprintf(
			'%s -g %s -filetype gen -og %s -ofiletype %s -bgen-bits %d 2>/dev/null',
			args$qctool,
			filename,
			bgen.filename,
			sprintf( "bgen_%s", args$bgen_version ),
			args$bits
		)
	)

	system(
		sprintf(
			'%s -g %s -og %s -precision 100 2>/dev/null',
			args$qctool,
			bgen.filename,
			gen.filename
		)
	)
	
	reconstructed.data = read.table( gen.filename, head = F, as.is=T )

	stopifnot( nrow( original.data ) == nrow( reconstructed.data ))
	stopifnot( ncol( original.data ) == ncol( reconstructed.data ))
	for( col in 1:6 ) {
		stopifnot( length( which( original.data[,col] != reconstructed.data[,col] )) == 0 )
	}
	max.discrepancy = max( original.data[,7:ncol( original.data)] - reconstructed.data[,7:ncol( reconstructed.data)] )
	bound = 1.0/((2^(args$bits))-1)
	if( args$verbose ) {
		cat( sprintf( "%s\n", gen.filename ))
		cat( sprintf( "simulation %d of %d, %d samples, %d variants, max discrepancy is %g, bound is %g.\n", i, args$iterations, N, args$variants, max.discrepancy, bound ))
	}
	stopifnot( max.discrepancy <= 1.0/((2^(args$bits))-1))
}
cat( "Success." )

