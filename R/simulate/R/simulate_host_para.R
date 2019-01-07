library( argparse )
parser <- ArgumentParser( description = 'Generate and test snp and samples stats for randomly generated genotype files')
parser$add_argument(
	"--samples",
	type = "integer",
	nargs = 1,
	help = "Number of samples to simulate",
	default = 1000
)
parser$add_argument(
	"--variants",
	type = "integer",
	nargs = 1,
	help = "Number of variants to simulate",
	default = 200
)
parser$add_argument(
	"--interactors",
	type = "integer",
	nargs = 1,
	help = "Number of interactions to include in simulation",
	default = 100
)
parser$add_argument(
	"--fakehomozygotes",
	action = 'store_true',
	default = FALSE,
	help = "Sample 1 allele for each parasite then report as homozygous"
)
parser$add_argument(
	"--tag",
	type = "character",
	nargs = 1,
	help = "Tag to add to filename"
)
parser$add_argument(
	"--covariates",
	type = "integer",
	nargs = 1,
	help = "Number of (std normal) covariates",
	default = 0
)
opts = parser$parse_args()

N = opts$samples
L = opts$variants

host = matrix( NA, nrow = L, ncol = N )
covariates = matrix( rnorm( N*opts$covariates ), nrow = N, ncol = opts$covariates )
covariate.betas = matrix( rnorm( opts$covariates*L, mean = 0, sd = 1 ), byrow = T, nrow = L, ncol = ncol( covariates ))
covariate.host.betas = covariate.betas
for( i in 1:nrow( covariate.betas )) {
	# Simulate confounding by having effect of covariate on host as well
	x = rbeta( 1, shape2 = 1, shape1 = 20 )
	covariate.host.betas[i,] = ( covariate.host.betas[1,] * x )
}
colnames( covariate.betas ) = sprintf( "cov%d.beta", 1:ncol( covariates ))
para = matrix( NA, nrow = L, ncol = N )

# Simulate host genotypes
# Frequency of the b allle
fHost = rbeta( L, shape1 = 1, shape2 = 10 )

# Frequency in host = 0 parasites
fParaH0 = rbeta( L, shape1 = 1, shape2 = 10 )

interactions = list(
	host = 1:L,
	para = 1:L
)

simulate <- function( N, f.baseline, spec, size = 2 ) {
	logodds.baseline = log( f.baseline / ( 1 - f.baseline ))
	x = rep( logodds.baseline, N ) ;
	for( i in 1:length( spec )) {
		x = x + spec[[i]]$value * spec[[i]]$beta
	}
	f = exp( x ) / ( 1 + exp( x ))
	sapply( f, function(p) { rbinom( 1, prob = p, size = size ) } )
}

host.spec = data.frame()
parasite.spec = data.frame()
for( i in 1:L ) {
	host.covariate.spec = lapply( 1:ncol( covariates ), function(j) { list( beta = covariate.host.betas[i,j], value = covariates[,j] ) } )
	host[i,] = simulate( N, fHost[i], host.covariate.spec )
	para.covariate.spec = lapply( 1:ncol( covariates ), function(j) { list( beta = covariate.betas[i,j], value = covariates[,j] ) } )
	if( opts$fakehomozygotes ) {
		para[i,] = simulate(
			N, fParaH0[i], para.covariate.spec, size = 1
		) * 2
	} else {
		para[i,] = simulate(
			N, fParaH0[i], para.covariate.spec, size = 2
		)
	}
	host.spec = rbind(
		host.spec,
		data.frame(
			variant = sprintf( "H%d", i ),
			frequency = fHost[i],
			interactors = NA,
			beta = NA
		)
	)
	parasite.spec = rbind(
		parasite.spec,
		data.frame(
			variant = sprintf( "P%d", i ),
			frequency = fParaH0[i],
			interactors = NA
		)
	)
}
host.spec = cbind( host.spec, covariate.betas )
for( i in 1:opts$interactors ) {
	host.i = interactions$host[i]
	para.i = interactions$para[i]
	beta = rnorm( 1, mean = 0, sd = 1 )
	predictor.spec = c(
		list( list( beta = beta, value = host[host.i,] )),
		lapply( 1:ncol( covariates ), function(j) { list( beta = covariate.betas[para.i,j], value = covariates[,j] ) } )
	)
	if( opts$fakehomozygotes ) {
		para[para.i,] = simulate(
			N, fParaH0[para.i], predictor.spec, size = 1
		) * 2
	} else {
		para[para.i,] = simulate(
			N, fParaH0[para.i], predictor.spec, size = 2
		)
	}
	host.spec[host.i,]$interactors = sprintf( "P%d", para.i )
	host.spec[host.i,]$beta = beta
	parasite.spec[para.i,]$interactors = sprintf( "H%d", host.i )
}

colnames( host ) = sprintf( "sample_%d", 1:N )
colnames( para ) = sprintf( "sample_%d", 1:N )

# Add sporatdic missingness to both host and parasite
para[ sample( 1:(nrow(para) * ncol(para)), nrow(para) * ncol(para) * 0.01 ) ] = "./."
host[ sample( 1:(nrow(host) * ncol(host)), nrow(para) * ncol(para) * 0.01 ) ] = "./."

samples = cbind(
	data.frame(
		ID = c( sprintf( "sample_%d", 1:N ))
	),
	covariates
)
if( ncol( covariates ) > 0 ) {
	colnames( samples ) = c( "ID", sprintf( "cov%d", 1:ncol( covariates )))
}

write.vcf <- function(
	samples,
	variants,
	genotypes,
	filename
) {
	cat( "##fileformat=VCFv4.2\n", file = filename )
	cat( "##FORMAT=<ID=GT,Type=String,Number=1,Description=\"Genotype\">\n", file = filename, append = T )
	cat(
		sprintf(
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n",
			paste( samples, collapse = "\t" )
		),
		file = filename,
		append = T
	)
	mode( genotypes ) = "character"
	genotypes[ which(genotypes == 0) ] = "0/0"
	genotypes[ which(genotypes == 1) ] = "0/1"
	genotypes[ which(genotypes == 2) ] = "1/1"
	write.table(
		cbind(
			data.frame(
				CHROM = variants$chromosome,
				POS = variants$position,
				ID = variants$rsid,
				REF = variants$alleleA,
				ALT = variants$alleleB,
				QUAL = ".",
				FILTER = ".",
				INFO = ".",
				FORMAT = "GT"
			),
			genotypes
		),
		file = filename,
		append = T,
		col.names = F,
		row.names = F,
		quote = F,
		sep = "\t"
	)
}

if( is.null( opts$tag )) {
	tag = ""
} else {
	tag = sprintf( ".%s", opts$tag )
}


write.vcf(
	samples = samples$ID,
	variants = data.frame(
		chromosome = "H1",
		SNPID = sprintf( "H%d", 1:L ),
		rsid = sprintf( "H%d", 1:L ),
		position = 1:L,
		alleleA = "A",
		alleleB = "G"
	),
	host,
	sprintf( "host%s.vcf", tag )
)

write.vcf(
	samples = samples$ID,
	variants = data.frame(
		chromosome = "P1",
		SNPID = sprintf( "P%d", 1:L ),
		rsid = sprintf( "P%d", 1:L ),
		position = 1:L,
		alleleA = "C",
		alleleB = "T"
	),
	para,
	sprintf( "parasite%s.vcf", tag )
)

write( colnames( samples ), file = sprintf( "samples%s.sample", tag ), ncol = ncol( samples ) )
colnames( samples )[1] = "0"
if( ncol( samples ) > 1 ) {
	colnames( samples )[2:ncol( samples )] = "C"
}
write.table(
	samples,
	file = sprintf( "samples%s.sample", tag ),
	row.names = F,
	col.names = T,
	quote = F,
	append = T
)

write.table( host.spec, file = sprintf( 'host_spec%s.txt', tag ), col.names = T, row.names = F, quote = F )
write.table( parasite.spec, file = sprintf( 'parasite_spec%s.txt', tag ), col.names = T, row.names = F, quote = F )
