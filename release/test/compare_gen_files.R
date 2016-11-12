library( argparse )
parser <- ArgumentParser( description = 'Compare two gen files to given precision' )
parser$add_argument(
	"--file1",
	type = "character",
	required = TRUE
)
parser$add_argument(
	"--file2",
	type = "character",
	required = TRUE
)
parser$add_argument(
	"--precision",
	type = "integer",
	nargs = 1,
	help = "Decimal places of precision",
	default = 5
)

args = parser$parse_args()

X = read.table( args$file1, hea=F, as.is=T )
Y = read.table( args$file2, hea=F, as.is=T )

XG = as.matrix( X[,7:ncol(X)] ); mode( XG ) = "numeric"
YG = as.matrix( Y[,7:ncol(Y)] ); mode( YG ) = "numeric"

w = which( X[,1:6] != Y[,1:6] )
w2 = which( abs( XG - YG ) > 10^-args$precision )

if( length(w) > 0 ) {
	cat( "File identifying data differs.\n" )
}

if( length(w2 ) > 0 ) {
	cat( "File genotype data differs.\n" )
	w2 = which( abs( XG - YG ) > 10^-args$precision, arr.in = TRUE )
	print(w2)
}

stopifnot( length(w) == 0 && length(w2) == 0 )
cat( "Success.\n" )
quit( "no" )
