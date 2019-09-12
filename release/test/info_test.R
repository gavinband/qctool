
compute.info <- function( genotypes, ploidy ) {
	stopifnot(length( which( !ploidy %in% c( 1, 2 ))) == 0 )

	# replace any missing data with 0 probabilities...
	w = which( is.na( rowSums( genotypes )))
	genotypes[w,] == 0

	# Handle missing genotypes by filling in distribution
	X = genotypes
	rowsums = rowSums( genotypes )
	theta = ( sum(X[,2]) + 2 * sum(X[,3]) ) / ( 2 * sum( X ))
	prior = c( (1-theta)^2, 2 * theta * (1-theta), theta^2 )
	for( i in 1:nrow(X)) {
		X[i,] = X[i,] + (1-rowsums[i]) * prior
	}
	stopifnot( length( which( abs( rowSums(X) - 1 ) > 1E-10 )) == 0 )

	impute.expected.variance = 2 * theta * (1 - theta)
	impute.mean = (genotypes[,2] + 2 * genotypes[,3])
	impute.variance = (genotypes[,2] + 4 * genotypes[,3]) - impute.mean^2
	impute.N = sum( genotypes )

	new.expected.variance = rep( 2 * theta * (1 - theta), nrow( X ))
	new.expected.variance[ ploidy == 1 ] = theta * ( 1 - theta )
	new.mean = (X[,2] + 2 * X[,3])
	new.variance = (X[,2] + 4 * X[,3]) - new.mean^2
	new.N = nrow(X)

	return(
		list(
			frequency = theta,
			missingness = sum(1 - rowSums(genotypes)) / nrow( genotypes ),
			impute.info = 1 - sum( impute.variance / ( impute.N * impute.expected.variance )),
			info = 1 - sum( new.variance / ( new.N * new.expected.variance ))
		)
	)
}

simulate.genotypes <- function( number.of.samples = 10, number.of.snps = 10 ) {
	genotypes = matrix( nrow = number.of.snps, ncol = 3 * number.of.samples )
	for( i in 1:number.of.snps ) {
		p = rbeta( 1, 1.2, 1.2 )
		for( j in 1:number.of.samples ) {
			shape1 = 20 * p
			shape2 = 20 * (1-p)
			q = rbeta( 1, shape1 * 5, shape2 * 5 )
			probs = c( q^2, 2 * q * (1-q), (1-q)^2 )
			genotypes[i, (j*3)+(-2:0)] = probs
		}
	}

	missingness = rbeta( nrow( genotypes ) * number.of.samples, 1, 4 )

	missing.proportion = rbinom( length(missingness), 1, p = missingness )

	w = which( runif( length( missing.proportion )) < 0.5 )

	missing.proportion[ w ] = missing.proportion[w] / 2
	for( i in 1:number.of.snps ) {
		for( j in 1:number.of.samples ) {
			genotypes[i,(j*3)+(-2:0)] = genotypes[i,(j*3)+(-2:0)] * (1-missing.proportion[i*j])
		}
	}
	return(
		list(
			samples = sprintf( "sample_%d", 1:number.of.samples ),
			variants = data.frame(
				chromosome = '01',
				SNPID = sprintf( "SNP_%d", 1:number.of.snps ),
				rsid = sprintf( "RS_%d", 1:number.of.snps ),
				position = 1:number.of.snps,
				alleleA = 'A',
				alleleB = 'G'
			),
			genotypes = genotypes
		)
	)
}

library( argparse )
parser <- ArgumentParser( description = 'Generate and test snp and samples stats for randomly generated genotype files')
parser$add_argument(
	"--samples",
	type = "integer",
	nargs = 1,
	help = "Number of samples to simulate",
	default = 100
)
parser$add_argument(
	"--variants",
	type = "integer",
	nargs = 1,
	help = "Number of samples to simulate",
	default = 1000
)
parser$add_argument(
	"--qctool",
	type = "character",
	nargs = 1,
	help = "Path to qctool executable",
	default = "../../build/release/qctool_v2.0-dev"
)
opts = parser$parse_args()

options(width=200)

data = simulate.genotypes( opts$samples, opts$variants )

genfile = tempfile()
write.table(
	cbind( data$variants, data$genotypes ),
	genfile,
	col.names = F,
	row.names = F,
	sep = ' '
)
statsfile = tempfile()
system( sprintf( "%s -g %s -snp-stats -osnp %s -filetype gen", opts$qctool, genfile, statsfile ))
qctool.stats = read.table( statsfile, head = T, as.is=T )

stats = data.frame()
for( i in 1:nrow( data$variants )) {
	genotypes = matrix(
		data$genotypes[i,],
		ncol = 3,
		byrow = T
	)
	info = compute.info( genotypes, rep( 2, length( data$samples ) ))
	stats = rbind(
		stats,
		cbind(
			data$variants[i,],
			frequency = info$frequency,
			missingness = info$missingness,
			info = info$info,
			impute_info = info$impute.info
		)
	)
}

cat( "Comparing frequency...\n" )
stopifnot( length( which( abs( stats$frequency - qctool.stats$alleleB_frequency ) > 1E-4 )) == 0 )
cat( "Comparing IMPUTE info...\n" )
stopifnot( length( which( abs( stats$impute_info - qctool.stats$impute_info ) > 1E-4 )) == 0 )
cat( "Comparing QCTOOL info...\n" )
stopifnot( length( which( abs( stats$info - qctool.stats$info ) > 1E-4 )) == 0 )

cat( "Success.\n" )
