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

opts = parser$parse_args()

chromosomes=c( sprintf( "%02d", 1:22 ), sprintf( "%d", 1:22 ), "other" )
V = data.frame(
	chromosome = sample( chromosomes, opts$variants, replace = T ),
	SNPID = sprintf( "SNP%d", 1:opts$variants ),
	rsid = sprintf( "rs%d", 1:opts$variants ),
	position = as.integer( round( runif( opts$variants, min = 0, max = 1E6 ))),
	alleleA = sample( c( 'A', 'C', 'T', 'G' ), opts$variants, replace = T ),
	alleleB = sample( c( 'A', 'C', 'T', 'G' ), opts$variants, replace = T )
)
for( i in 1:opts$iterations ) {
	N = sample( 1:opts$max_samples, 1 )
	G = matrix( NA, nrow = opts$variants, ncol = N*2 )
	G[,] = runif( nrow(G) * ncol(G) )
	omit.chromosome = ( runif(1) > 0.5 )
	cols = 1:6
	if( omit.chromosome ) {
		cols = 2:6
	}
	filename = tempfile()
	write.table(
		cbind( V[2:6], G ),
		file = filename,
		col.names = F,
		row.names = F,
		quote = F,
		sep = " "
	)
	
	cat( sprintf( "Iteration %d, %dx%d...\n", i, nrow(G), ncol(G) ))
	
	# Convert it to bgen
	cmd1 = sprintf(
		paste(
			sep = ' ',
			opts$qctool,
			'-g %s',
			'-og %s.bgen',
			'-filetype impute_allele_probs'
		),
		filename,
		filename
	)
	if( opts$verbose ) {
		cat( sprintf( "Converting %s to bgen...\n", filename ))
	}
	system( sprintf( "%s 2>/dev/null", cmd1 ) )

	# Convert it back to allele probs
	filename2 = tempfile()
	cmd2 = sprintf(
		paste(
			sep = ' ',
			opts$qctool,
			'-g %s.bgen',
			'-og %s',
			'-ofiletype impute_allele_probs'
		),
		filename,
		filename2
	)
	
	if( omit.chromosome ) {
		cmd2 = sprintf( "%s -omit-chromosome", cmd2 )
	}
	if( opts$verbose ) {
		cat( sprintf( "Converting bgen back to %s...\n", filename2 ))
	}
	system( sprintf( "%s 2>/dev/null", cmd2 ))
	
	if( opts$verbose ) {
		cat( sprintf( "Reading %s...\n", filename2 ))
	}
	result = read.table( filename2, header = F, as.is = T, sep = " " )
	stopifnot( length( which( result[,1:length(cols)] != V[,cols] )) == 0 )
	G2 = as.matrix( result[, (length(cols)+1):ncol(result)])
	mode(G2) = "numeric"
	stopifnot( length( which( abs( G2 - G ) > 1E-6 )) > 0 )
}



