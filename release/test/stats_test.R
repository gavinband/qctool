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

filename=tempfile( tmpdir = "/tmp/" )

genotypes = list(
    diploid = c( '0/0', '0/1', '1/1', './.' ),
    haploid = c( '0', '1', '.' )
)

base.frequency = rbeta( opts$variants, 1, 10 )
missing.rate = rbeta( opts$variants, 1, 10 )

base.frequencies = list(
    diploid = matrix(
        nrow = opts$variants, ncol = length( genotypes$diploid ),
        dimnames = list( sprintf( "V%d", 1:opts$variants ), genotypes$diploid  )
    ),
    haploid = matrix(
        nrow = opts$variants, ncol = length( genotypes$haploid ),
        dimnames = list( sprintf( "V%d", 1:opts$variants ), genotypes$haploid  )
    )
)
base.frequencies[['haploid']][,'0'] = (1-base.frequency) * (1-missing.rate)
base.frequencies[['haploid']][,'1'] = (1-base.frequency) * (1-missing.rate)
base.frequencies[['haploid']][,'.'] = missing.rate
base.frequencies[['diploid']][,'0/0'] = (1-base.frequency)^2 * (1-missing.rate)
base.frequencies[['diploid']][,'0/1'] = 2*base.frequency*(1-base.frequency) * (1-missing.rate)
base.frequencies[['diploid']][,'1/1'] = base.frequency^2 * (1-missing.rate)
base.frequencies[['diploid']][,'./.'] = missing.rate

#
gender = c( 'M', 'F' )[ rbinom( opts$samples, p = 0.5, size = 1 ) + 1 ]
ploidy = list(
    'autosomes' = c( 'M' = 'diploid', 'F' = 'diploid' )[ gender ],
    'X' = c( 'M' = 'haploid', 'F' = 'diploid' )[ gender ],
    '0X' = c( 'M' = 'haploid', 'F' = 'diploid' )[ gender ],
    'Y' = c( 'M' = 'haploid', 'F' = 'zeroploid' )[ gender ],
    '0Y' = c( 'M' = 'haploid', 'F' = 'zeroploid' )[ gender ]
)
chromosomes = sample( c( rep( "autosomes", 10 ), "X", "0X", "Y", "0Y" ), size = opts$variants, replace = T )

generate <- function( n, prob, levels ) {
    a = rmultinom( n, 1, p = prob )
    sapply( 1:ncol(a), function(i) { levels[which( a[,i] == 1 )]})
}

G = matrix( '.', nrow = opts$variants, ncol = opts$samples )
for( i in 1:opts$variants ) {
    chromosome = chromosomes[i]
    for( ploid in c( 'haploid', 'diploid' )) {
        w = which( ploidy[[ chromosome ]] == ploid )
        if( length(w) > 0 ) {
            G[i,w] = generate( length(w), base.frequencies[[ploid]][i,], genotypes[[ ploid ]] )
        }
    }
}

stats = data.frame(
    alleleA_count = rowSums( G == '0/1' | G == '0' ) + 2*rowSums( G == '0/0' ),
    alleleB_count = rowSums( G == '0/1' | G == '1' ) + 2*rowSums( G == '1/1' ),
    total = rowSums( G != './.' & G != '.' )
)
stats$alleleA_frequency = stats$alleleA_count / ( stats$alleleA_count + stats$alleleB_count )
stats$alleleB_frequency = stats$alleleB_count / ( stats$alleleA_count + stats$alleleB_count )

############################
# Now prepare and run qctool
############################
vcf.filename = tempfile( tmpdir = "/tmp/", fileext = '.vcf' )
sample.filename = tempfile( tmpdir = "/tmp/", fileext = '.sample' )
snp.stats.filename = tempfile( tmpdir = "/tmp/", fileext = '.tsv' )

VAR = data.frame(
    `#CHROM` = as.character( sample( c( sprintf( "%02d", 1:22 ), sprintf( "%d", 1:22 )), nrow( G ), replace = T ) ),
    `POS` = 1:opts$variants,
    `ID` = sprintf( 'V%d', 1:opts$variants ),
    `REF` = 'A',
    `ALT` = 'G',
    `QUAL` = '.',
    `FILTER` = '.',
    `INFO` = '.',
    `FORMAT` = 'GT',
    check.names = F,
    stringsAsFactors = F
)
VAR[which(chromosomes != "autosomes"),1] = chromosomes[ which( chromosomes != "autosomes" ) ]

write( '##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Type=String,Number=1,Description="genotype">', file = vcf.filename )
write.table( cbind( VAR, G ), file = vcf.filename, sep = '\t', col.names = T, row.names = F, quote = F, append = T )

write( 'ID_1 sex\n0 D', file = sample.filename )
write.table(
    data.frame(
        ID_1 = sprintf( "sample_%d", 1:ncol(G)),
        sex = gender
    ),
    sample.filename, append = T, col.names = F, row.names = F, quote = F
)

cmd = sprintf(
    "%s -g %s -s %s -snp-stats -osnp %s",
    opts$qctool,
    vcf.filename, sample.filename, snp.stats.filename
)

system( cmd )

cat( "I simulated variants on these chromosomes:\n" )
table( VAR$`#CHROM` )
############################
# Get and compare qctool results
############################

compare.floats <- function( a, b, tolerance = 1E-6 ) {
    return( abs( a - b ) < tolerance )
}


results = read.table( snp.stats.filename, sep = '\t', comment = '#', header = T )

columns = c( "alleleA_count", "alleleB_count", "alleleA_frequency", "alleleB_frequency", "A", "B", "AA", "AB", "BB" )
failed.columns = c()
for( column in columns ) {
    tryCatch({
        w = which( is.na( results[, column] ) != is.na( stats[, column ] ) | !compare.floats( results[,column], stats[,column] ) )
        if( length( w ) > 0 ) {
            failed.columns = c( failed.columns, column )
            cat( sprintf( "Expected results for column \"%s\" differ, differences are:\n", column ))
            print( cbind( results[w,1:6], expected = stats[w,column], observed = results[w,column] ) )
        }
    }, error = function(e) {
        failed.columns = c( failed.columns, column )
    })
}
