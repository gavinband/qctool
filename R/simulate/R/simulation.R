library( RSQLite )
qqplot <- function( ps, plot.new = TRUE, limit = 50, CI = 0.95, thin = TRUE, ... ) {
    cat( "Sorting p-values...\n" )
    ps = sort(ps)
    lambda = -log10( median( ps ) )
    lambda = lambda / -log10( 0.5 )

    cat( "Converting to log10 scale...\n" )
    qs = -log10( ps )
    expected = (1:length(ps))/(length(ps)+1)
    expected = -log10( expected )

    cat( "Computing quantiles...\n" )
    lower.quantile = qbeta( 1-(1-CI)/2, (1:length(ps)), length(ps):1 )
    upper.quantile = qbeta( (1-CI)/2, (1:length(ps)), length(ps):1 )
    lower.quantile = -log10( lower.quantile )
    upper.quantile = -log10( upper.quantile )
    # We're not much interested in anything that lies near the diagonal line.  So don't plot them.
    cat( "Thinning...\n" )
    if( thin ) {
        wDiag = which( abs( qs - expected ) < 0.1 & qs < 5 )
        wNotDiag = which( qs >= 5 | abs( qs - expected ) >= 0.1 )
        wDiag = sample( wDiag, min( 10000, length( wDiag ) ))
        o = order( qs[ c( wDiag, wNotDiag )] )
        qs = qs[ c( wDiag, wNotDiag ) ][o]
        expected = expected[ c( wDiag, wNotDiag ) ][o]
        lower.quantile = lower.quantile[ c( wDiag, wNotDiag ) ][o]
        upper.quantile = upper.quantile[ c( wDiag, wNotDiag ) ][o]
    }
    cat( "Plotting...\n" )
    pch = rep( 46, length( qs )) # = '.'
    pch[ which( qs > upper.quantile | qs < lower.quantile )] = 1
    if( plot.new ) {
        plot( pmin( expected, limit ), pmin( qs, limit ), pch = pch, ... )
        points( pmin( expected, limit ), pmin( lower.quantile, limit ), type = 'l', lty = 2 )
        points( pmin( expected, limit ), pmin( upper.quantile, limit ), type = 'l', lty = 2 )
    } else {
        points( pmin( expected, limit ), pmin( qs, limit ), pch = pch, ... )
        points( pmin( expected, limit ), pmin( lower.quantile, limit ), type = 'l', lty = 2 )
        points( pmin( expected, limit ), pmin( upper.quantile, limit ), type = 'l', lty = 2 )
    }
    abline( a = 0, b = 1, col = "red")
    return( lambda )
}

prob.phenotype.given.genotype <- function(
	mu,
	genotypes,
	genotype.levels,
	effect.sizes,
	pcs = NULL,
	pc.effects = NULL
) {
	# If there are N phenotypes, mu has length N and effect.sizes has length 2N
	# First N are additive effects, followed by heterozygote effects.
	
	# If there are d pcs, pc.effects is a dxN matrix with pc effects for each outcome.
	stopifnot( length( effect.sizes ) == 2 * length( mu ))
	result = matrix( nrow = length( genotypes ), ncol = length( mu ) )
	for( i in 1:length( mu )) {
		predictor = (
			mu[i]
			+ ( effect.sizes[i] * genotype.levels[ genotypes + 1 ] )
			+ ( effect.sizes[i + length(mu) ] * ( genotypes == 1 ) )
		)
		if( !is.null( pcs )) {
			predictor = predictor + ( pcs %*% pc.effects )
		}
		result[,i] = exp( predictor )
	}
	for( i in 1:nrow( result )) {
		result[i,] = result[i,] / sum( result[i,] ) ;
	}
	return( result )
}

prob.genotype.given.phenotype <- function( phenotypes, mu, effect.sizes, population.frequencies, pcs = NULL, pc.effects = NULL, verbose = FALSE ) {
#	P( G | D ) = P( D | G ) P( G ), normalised by P(D) = sum of the top half.
	stopifnot( length( effect.sizes ) == 2 * length( mu ) )
	stopifnot( length( population.frequencies ) == 3 )
	N = length( phenotypes )
	# Make an indexing matrix
	phenotypes = matrix( c( 1:N, phenotypes ), ncol = 2 )
	# Make the result
	result = matrix( NA, nrow = N, ncol = 3 )
	# mean genotype
	mean.genotype = population.frequencies[2] + 2 * population.frequencies[3]
	genotype.levels = 0:2 - mean.genotype
	# Make phenotypes into an indexing matrix.
	for( g in 0:2 ) {
		# Compute probability of each phenotype
		x = prob.phenotype.given.genotype( mu, rep( g, N ), genotype.levels, effect.sizes, pcs, pc.effects ) ;
		if( verbose ) {
			print(x)
			print( phenotypes )
		}
		# Compute probability of the given phenotype for each sample
		x = x[ phenotypes ]
		# Compute unnormalised probability of the genotype given phenotype
		result[,g+1] = x * population.frequencies[g+1] ;
	}
#	print( result )
	for( i in 1:nrow( result )) {
		result[i,] = result[i,] / sum( result[i,] ) ;
	}
	return( result )
}

# mu[i] is the log-odds of phenotype i for an average genotype.
# Let's put this pretty low for disease cases.

numberOfSamples = 2000 # number of cases and number of controls
phenotypes = c( "control", "case1", "case2", "case3" )
phenotype = factor( sample( phenotypes, numberOfSamples, replace = T, prob = c( 0.5, 0.5/3, 0.5/3, 0.5/3) ), levels = phenotypes )

numberOfNullSNPs = 1000
numberOfEffectSNPs = 200
numberOfSNPs = numberOfNullSNPs + numberOfEffectSNPs

frequency = rbeta( numberOfSNPs, shape1 = 1.1, shape2 = 5 )
pop.frequencies = matrix(
	c( ( 1 - frequency )^2, 2 * frequency * ( 1 - frequency ), frequency^2 ),
	nrow = numberOfSNPs,
	ncol = 3
)

effect.size.scale = 0.5
effect.sizes = matrix(
	0,
	nrow = numberOfEffectSNPs,
	ncol = 2 * length( phenotypes ),
	dimnames = list(
		c( sprintf( "binomial%d", 1:100 ), sprintf( "multinomial%d", 1:100 ) ),
		c( sprintf( "%s:add", phenotypes ), sprintf( "%s:het", phenotypes ))
	)
)
# 100 SNPs that are case/control
for( i in 1:100 ) {
	w = 2:length(phenotypes)
	effect = rbeta( 1, shape1 = 4, shape2 = 3 ) * effect.size.scale
	effect.sizes[i,w] = effect
	# 
}
# 5 of those have non-additive effects
for( i in 95:100 ) {
	wAdd = 1:length(phenotypes)
	wHet = wAdd + length(phenotypes)
	a = runif(1)
	if( a < 0.3333 ) {
		# dom model
		effect.sizes[i,wHet] = effect.sizes[i,wAdd]
	} else if( a < 0.66666 ) {
		# dom model
		effect.sizes[i,wHet] = -effect.sizes[i,wAdd]
	} else {
		# het model
		effect.sizes[i,wHet] = effect.sizes[i,wAdd]
		effect.sizes[i,wAdd] = 0
	}
}

# 100 SNPs with independent effects
for( i in 101:200 ) {
	w = 2:length(phenotypes)
	effect.sizes[i,w] = rbeta( length(phenotypes)-1, shape1 = 4, shape2 = 3 ) * effect.size.scale
	a = runif( length( phenotypes ) - 1, -1, 1 )
	effect.sizes[i,w] = effect.sizes[i,w] * sign( a )
}
# 5 of those have non-additive effects
for( i in 195:200 ) {
	wAdd = 1:length(phenotypes)
	wHet = wAdd + length(phenotypes)
	a = runif(1)
	if( a < 0.3333 ) {
		# dom model
		effect.sizes[i,wHet] = effect.sizes[i,wAdd]
	} else if( a < 0.66666 ) {
		# dom model
		effect.sizes[i,wHet] = -effect.sizes[i,wAdd]
	} else {
		# het model
		effect.sizes[i,wHet] = effect.sizes[i,wAdd]
		effect.sizes[i,wAdd] = 0
	}
}

# 2000 NULL SNPs
effect.sizes = rbind(
	effect.sizes,
	matrix( 0, nrow = numberOfNullSNPs, ncol = (2 * length( phenotypes )))
)

# simulate two principal components
scale = 1.5
pcs = matrix( NA, nrow = length( phenotype ), ncol = 2 )
w = which( phenotype %in% c( "control" ))
pcs[ w, 1 ] = rbeta( length( w ), shape1 = 2*scale, shape2 = 5*scale ) + rnorm( length(w), mean = 0, sd = 0.1 )
pcs[ w, 2 ] = rbeta( length( w ), shape1 = 5*scale, shape2 = 2*scale ) + rnorm( length(w), mean = 0, sd = 0.1 )
w = which( phenotype %in% c( "case1", "case2", "case3" ))
pcs[ w, 1 ] = rbeta( w, shape1 = 5*scale, shape2 = 3*scale ) + rnorm( length(w), mean = 0, sd = 0.1 )
w = which( phenotype %in% c( "case2" ))
pcs[ w, 2 ] = rbeta( length( w ), shape1 = 4*scale, shape2 = 3*scale ) + rnorm( length(w), mean = 0, sd = 0.1 )
w = which( phenotype %in% c( "case1", "case3" ))
pcs[ w, 2 ] = rbeta( length( w ), shape1 = 2*scale, shape2 = 5*scale ) + rnorm( length(w), mean = 0, sd = 0.1 )
dir.create( "images" )

pdf( file = "images/pcs.pdf", width = 4, height = 4 )
plot( pcs[,1], pcs[,2], col = phenotype, xlab = "PC 1", ylab = "PC 2", xlim = c( -0.2, 1.5 ), ylim = c( -0.2, 1.5 ) )
legend( "topright", col = 1:length( phenotype ), legend = levels( phenotype ), pch = 19, bty = 'n' )
dev.off()

pc.effects = matrix( c( -0.3, 0.2 ), ncol = 1 )

mu = c( 0, -3, -3, -3 )
# Show the phenotype probs for each genotype:
print( prob.phenotype.given.genotype( mu, 0, 0:2, c(0,0,0,0,0,0,0,0) ) ) 
print( prob.genotype.given.phenotype( 1:length(phenotypes), mu, c(0,0,0,0,0,0,0,0), pop.frequencies[1,] ) )
print( prob.genotype.given.phenotype( 1:length(phenotypes), mu, c(0,0,1,0,0,0,0,1), pop.frequencies[1,] ) )

# Now simulate genotypes.
G = matrix( NA, nrow = numberOfSNPs, ncol = numberOfSamples )
for( i in 1:nrow(G) ) {
	cat( i, "of", nrow(G), "...\n" )
	p = prob.genotype.given.phenotype( phenotype, mu, effect.sizes[i,], pop.frequencies[i,], pcs, pc.effects )
	genotypes = sapply( 1:numberOfSamples, function(i) { rmultinom( 1, 1, prob = p[i,] ) } )
	G[i,] = sapply( 1:ncol( genotypes ), function(i) { which( genotypes[,i] == 1 ) } ) - 1
}
sample.frequency = rowSums( G ) / ( 2 * ncol(G) )
plot( sample.frequency, frequency, xlim = c( 0, 1 ), ylim = c( 0, 1 ), col = ( rowSums( effect.sizes ) > 0 ) + 1 )


######################################
# Shall we do a quick check?
# First 100 are case/cnotrol
library( nnet )
g = glm( as.integer( phenotype != "control" ) ~ G[2,] + pcs, family="binomial"); summary(g)$coeff
g = multinom( phenotype ~ as.factor( G[1,] ) + pcs ); summary(g)
g = multinom( phenotype ~ G[101,] ); summary(g)
g = multinom( phenotype ~ G[201,] ); summary(g)


######################################
# Write a sample file and a vcf file
dir.create( "data" )
version = "v3"
output.file = sprintf( "data/test_data_%s.vcf", version )
sample.filename = sprintf( "data/test_data_%s.sample", version )
samples = data.frame(
	ID_1 = sprintf( "sample_%d", 1:ncol( vcf )),
	ID_2 = sprintf( "sample_%d", 1:ncol( vcf )),
	missing = 0,
	case = as.integer( phenotype != 'control' ),
	sub = phenotype,
	pc1 = pcs[,1],
	pc2 = pcs[,1]
)

write( colnames( samples ), sample.filename, ncolumns = 7 )
write( c( 0, 0, 0, 'B', 'D', 'C', 'C' ), sample.filename, ncolumns = 7, append = T )
write.table( samples, file = sample.filename, col.names = F, row.names = F, quote = F, append = T )

# Write a VCF file with the results...
vcf = G
mode( vcf ) = "character"
vcf[ which( G == 0 )] = '0/0'
vcf[ which( G == 1 )] = '0/1'
vcf[ which( G == 2 )] = '1/1'

colnames( vcf ) = sprintf( "sample_%d", 1:ncol( vcf ))

variants = data.frame(
	CHROM = '22',
	POS = 1:numberOfSNPs,
	ID = c( sprintf( "binomial%d", 1:94 ), sprintf( "binomial%d_nonadd", 95:100 ), sprintf( "multinomial%d", 1:94 ), sprintf( "multinomial_nonadd%d", 95:100 ), sprintf( "null%d", 1:numberOfNullSNPs ) ),
	REF = 'A',
	ALT = 'G',
	QUAL = '.',
	FILTER = '.',
	INFO = '.',
	FORMAT = 'GT'
)
write( "##fileformat=VCFv4.1", output.file )
write( "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype calls\">", output.file, append = T )
write.table( cbind( variants, vcf ), file = output.file, quote = F, append = T, col.names = T, row.names = F, sep = "\t" )

######################################
# Run SNPTEST
dir.create( 'scans' )
system(
	paste(
		'snptest_v2.6-dev -data data/test_data_v3.vcf data/test_data_v3.sample -genotype_field GT',
		'-cov_names pc1 pc2',
		'-method newml -frequentist add -pheno case',
		'-o scans/test_data_case_v3.snptest',
		'-max_iterations 100'
	)
)
system(
	paste(
		'snptest_v2.6-dev -data data/test_data_v3.vcf data/test_data_v3.sample -genotype_field GT',
		'-cov_names pc1 pc2',
		'-method newml -frequentist add -pheno case',
		'-full_parameter_estimates -o scans/test_data_case_v3_full.snptest',
		'-max_iterations 100'
	)
)

system(
	paste(
		'snptest_v2.6-dev -data data/test_data_v3.vcf data/test_data_v3.sample -genotype_field GT',
		'-cov_names pc1 pc2',
		'-method newml -frequentist add -pheno sub',
		'-baseline_phenotype control -o scans/test_data_sub_v3.snptest',
		'-max_iterations 100'
	)
)

system(
	paste(
		'snptest_v2.6-dev -data data/test_data_v3.vcf data/test_data_v3.sample -genotype_field GT',
		'-cov_names pc1 pc2',
		'-method newml -frequentist gen -pheno sub ',
		'-baseline_phenotype control -o scans/test_data_sub_gen_v3.snptest',
		'-max_iterations 100'
	)
)

system(
	paste(
		'snptest_v2.6-dev -data data/test_data_v3.vcf data/test_data_v3.sample -genotype_field GT',
		'-cov_names pc1 pc2',
		'-method newml -frequentist gen -pheno sub ',
		'-baseline_phenotype control',
		'-o /tmp/snptest.out',
		'-max_iterations 100',
		'-debug',
		'-snpid binomial1'
	)
)

######################################
# Run Bingwa

dir.create( 'meta' )
file.remove( 'meta/test_data_v3.sqlite' )
system(
	paste(
		'bingwa_v2.0-dev -data scans/test_data_case_v3.snptest',
		'-cohort-names test',
		'-define-sd-set .=0.2,0.75,1',
		'-simple-prior rho=1/sd=[.]',
		'-simple-prior-name fix-add',
		'-o meta/test_data_v3.sqlite',
		'-analysis-name case',
		'-table-name CC',
		'-extra-columns comment'
	)
)

system(
	paste(
		'bingwa_v2.0-dev -data scans/test_data_sub_v3.snptest',
		'-cohort-names test',
		'-define-sd-set .=0.2,0.75,1',
		'-simple-prior rho=1/sd=[.],[.],[.]/cor=1,1,1',
		'-simple-prior rho=1/sd=[.],[.],[.]/cor=0,0,0',
		'-simple-prior-name fix-add',
		'-simple-prior-name ind-add',
		'-o meta/test_data_v3.sqlite',
		'-analysis-name "sub:add"',
		'-table-name subAdd',
		'-extra-columns comment'
	)
)


# add 1 0 1 0 1 0
# het   1 0 1 0 1
# add     1 0 1 0
# het       1 0 1
# add         1 0
# het           1
system(
	paste(
		'bingwa_v2.0-dev -data scans/test_data_sub_gen_v3.snptest',
		'-cohort-names test',
		'-define-sd-set .=0.2,0.75,1',
		'-simple-prior rho=1/sd=[.],0.01,[.],0.01,[.],0.01/cor=0,1,0,1,0,0,1,0,1,0,1,0,0,1,0',
		'-simple-prior-name fix-add',
		'-simple-prior rho=1/sd=[.],0.01,[.],0.01,[.],0.01/cor=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0',
		'-simple-prior-name ind-add',
		'-simple-prior rho=1/sd=[.],[.],[.],[.],[.],[.]/cor=0,1,0,1,0,0,1,0,1,0,1,0,0,1,0',
		'-simple-prior-name fix-gen',
		'-simple-prior rho=1/sd=[.],[.],[.],[.],[.],[.]/cor=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0',
		'-simple-prior-name ind-gen',
		'-o meta/test_data_v3.sqlite',
		'-analysis-name "sub:gen"',
		'-table-name subGen3',
		'-extra-columns comment'
		
	)
)

######################################
# Look at results
db = dbConnect( dbDriver( "SQLite" ), "meta/test_data_v3.sqlite" )
D = dbGetQuery(
	db,
	paste(
		'SELECT CC.rsid AS rsid, CC."test:B_allele_frequency" AS B_allele_frequency,',
		'CC."FixedEffect:pvalue" AS CC_meta_pvalue,',
		'"Bayesian:fix-add/sd=0.2/cor=:bf" AS CC_fix_add_small,',
		'"Bayesian:fix-add/sd=0.75/cor=:bf" AS CC_fix_add_medium,',
		'"Bayesian:fix-add/sd=1/cor=:bf" AS CC_fix_add_large,',
		'CC."test:comment" AS CC_comment,',
		'"Bayesian:fix-add/sd=0.2,0.2,0.2/cor=1,1,1:bf" AS SA_fix_add_small,',
		'"Bayesian:fix-add/sd=0.75,0.75,0.75/cor=1,1,1:bf" AS SA_fix_add_medium,',
		'"Bayesian:fix-add/sd=0.75,0.75,0.75/cor=1,1,1:bf" AS SA_fix_add_large,',
		'"Bayesian:ind-add/sd=0.2,0.2,0.2/cor=0,0,0:bf" AS SA_ind_add_small,',
		'"Bayesian:ind-add/sd=0.75,0.75,0.75/cor=0,0,0:bf" AS SA_ind_add_medium,',
		'"Bayesian:ind-add/sd=1,1,1/cor=0,0,0:bf" AS SA_ind_add_large,',
		'SA."test:comment" AS SA_comment,',
		'"Bayesian:fix-add/sd=0.2,0.01,0.2,0.01,0.2,0.01/cor=0,1,0,1,0,0,1,0,1,0,1,0,0,1,0:bf" AS SG_fix_add_small,',
		'"Bayesian:fix-add/sd=0.75,0.01,0.75,0.01,0.75,0.01/cor=0,1,0,1,0,0,1,0,1,0,1,0,0,1,0:bf" AS SG_fix_add_medium,',
		'"Bayesian:fix-add/sd=1,0.01,1,0.01,1,0.01/cor=0,1,0,1,0,0,1,0,1,0,1,0,0,1,0:bf" AS SG_fix_add_large,',
		'"Bayesian:ind-add/sd=0.2,0.01,0.2,0.01,0.2,0.01/cor=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0:bf" AS SG_ind_add_small,',
		'"Bayesian:ind-add/sd=0.75,0.01,0.75,0.01,0.75,0.01/cor=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0:bf" AS SG_ind_add_medium,',
		'"Bayesian:ind-add/sd=1,0.01,1,0.01,1,0.01/cor=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0:bf" AS SG_ind_add_large,',
		'"Bayesian:fix-gen/sd=0.2,0.2,0.2,0.2,0.2,0.2/cor=0,1,0,1,0,0,1,0,1,0,1,0,0,1,0:bf" AS SG_fix_gen_small,',
		'"Bayesian:fix-gen/sd=0.75,0.75,0.75,0.75,0.75,0.75/cor=0,1,0,1,0,0,1,0,1,0,1,0,0,1,0:bf" AS SG_fix_gen_medium,',
		'"Bayesian:fix-gen/sd=1,1,1,1,1,1/cor=0,1,0,1,0,0,1,0,1,0,1,0,0,1,0:bf" AS SG_fix_gen_large,',
		'"Bayesian:ind-gen/sd=0.2,0.2,0.2,0.2,0.2,0.2/cor=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0:bf" AS SG_ind_gen_small,',
		'"Bayesian:ind-gen/sd=0.75,0.75,0.75,0.75,0.75,0.75/cor=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0:bf" AS SG_ind_gen_medium,',
		'"Bayesian:ind-gen/sd=1,1,1,1,1,1/cor=0,0,0,0,0,0,0,0,0,0,0,0,0,0,0:bf" AS SG_ind_gen_large,',
		'SG."test:comment" AS SG_comment',
		'FROM CCView CC',
		'INNER JOIN subAdd SA ON SA.variant_id == CC.variant_id',
		'INNER JOIN subGen3 SG ON SG.variant_id == CC.variant_id'
	)
)

D$CC_fix_add = ( D$CC_fix_add_small + D$CC_fix_add_medium + D$CC_fix_add_large ) / 3
D$SA_fix_add = ( D$SA_fix_add_small + D$SA_fix_add_medium + D$SA_fix_add_large ) / 3
D$SA_ind_add = ( D$SA_ind_add_small + D$SA_ind_add_medium + D$SA_ind_add_large ) / 3
D$SG_fix_add = ( D$SG_fix_add_small + D$SG_fix_add_medium + D$SG_fix_add_large ) / 3
D$SG_fix_gen = ( D$SG_fix_gen_small + D$SG_fix_gen_medium + D$SG_fix_gen_large ) / 3
D$SG_ind_add = ( D$SG_ind_add_small + D$SG_ind_add_medium + D$SG_ind_add_large ) / 3
D$SG_ind_gen = ( D$SG_ind_gen_small + D$SG_ind_gen_medium + D$SG_ind_gen_large ) / 3

D$CC_mean_bf = D$CC_fix_add
D$SA_mean_bf = 0.8 * D$SA_fix_add + 0.2 * D$SA_ind_add
D$SG_mean_bf = 0.8 * D$SG_fix_add + 0.2 * D$SG_ind_add
D$SG_mean_bf = 0.8^2* D$SG_fix_add + 0.8*0.2 * D$SG_fix_gen + 0.2 * 0.8 * D$SG_ind_add + 0.2^2 * D$SG_ind_gen

D$maf = D$B_allele_frequency
D$maf = pmin( D$maf, 1 - D$maf )
D$maf_cut = cut( D$maf, breaks = c( 0, 0.1, 0.2, 0.3, 0.4, 0.5 ))
shapes = c( 1, 3, 4, 8 )
model.names = c( "null", "additive; case/control", "additive non case/control", "general" )
pch = rep( shapes[1], nrow(D) )
pch[ grep( "binomial", D$rsid ) ] = shapes[2]
pch[ grep( "multinomial", D$rsid ) ] = shapes[3]
pch[ grep( "nonadd", D$rsid ) ] = shapes[4]

pdf( file = "images/qqplot.pdf", width = 4, height = 4, xlab = "Exptected -log10( P-value )", ylab = "Observed -log10( P-value )" )
qqplot( D$CC_meta_pvalue[ grep( "null", D$rsid )], )
dev.off()

#
mar = c( 4, 5, 2, 2 )
width = 6
height = 5
pdf( file = "images/fixed_effect_pvalue_comparison.pdf", width = width, height = height )
par( mar = mar )
plot( log10( D$CC_fix_add ), -log10( D$CC_meta_pvalue ), col = D$maf_cut, pch = pch,
	xlab = "log10( fixed effect Bayes factor )",
	ylab = "-log10( Wald p-value )"
)
legend( "bottomright", col = 1:length( levels( D$maf_cut )), legend = levels( D$maf_cut ), pch = 19 )
legend( "topleft", legend = model.names, pch = shapes )
dev.off()

pdf( file = "images/fixed_effect_pvalue_comparison_zoom.pdf", width = width, height = height )
par( mar = mar )
plot( log10( D$CC_fix_add ), -log10( D$CC_meta_pvalue ), col = D$maf_cut, pch = pch, xlim = c( -1, 2 ), ylim = c( 0, 3 ),
	xlab = "log10( case/control Bayes factor )",
	ylab = "-log10( case/control Wald p-value )"
)
legend( "bottomright", col = 1:length( levels( D$maf_cut )), legend = levels( D$maf_cut ), pch = 19 )
legend( "topleft", legend = model.names, pch = shapes )
dev.off()

pdf( file = "images/cc_fix_add_vs_sa_fix_add.pdf", width = width, height = height )
par( mar = mar )
plot( log10( D$CC_fix_add ), log10( D$SA_fix_add ), col = D$maf_cut, pch = pch,
	xlab = "log10( additive c/c BF )\n(logit model)",
	ylab = "log10( additive c/c BF )\n(multinomial model)"
)
abline( a = 0, b = 1, col = 'red' )
legend( "bottomright", col = 1:length( levels( D$maf_cut )), legend = levels( D$maf_cut ), pch = 19 )
legend( "topleft", legend = model.names, pch = shapes )
dev.off()

w = 1:nrow(D); #which( D$SG_comment == 'NA' )
pdf( file = "images/cc_fix_add_vs_sg_fix_add.pdf", width = width, height = height )
par( mar = mar )
plot( log10( D$CC_fix_add[w] ), log10( D$SG_fix_add[w] ), col = D$maf_cut[w], pch = pch[w],
xlab = "log10( additive c/c BF )\n(logit model)",
ylab = "log10( additive c/c BF )\n(general multinomial model)"
)
abline( a = 0, b = 1, col = 'red' )
legend( "bottomright", col = 1:length( levels( D$maf_cut )), legend = levels( D$maf_cut ), pch = 19 )
legend( "topleft", legend = model.names, pch = shapes )
dev.off()

w = 1:nrow(D); #which( D$SG_comment == 'NA' )
pdf( file = "images/cc_fix_add_vs_sg_fix_gen.pdf", width = width, height = height )
par( mar = mar )
plot( log10( D$CC_fix_add[w] ), log10( D$SG_fix_gen[w] ), col = D$maf_cut[w], pch = pch[w],
xlab = "log10( additive c/c BF )\n(logit model)",
ylab = "log10( general c/c BF )\n(general multinomial model)"
)
abline( a = 0, b = 1, col = 'red' )
legend( "bottomright", col = 1:length( levels( D$maf_cut )), legend = levels( D$maf_cut ), pch = 19 )
legend( "topleft", legend = model.names, pch = shapes )
dev.off()

w = 1:nrow(D); #which( D$SG_comment == 'NA' )
pdf( file = "images/cc_fix_add_vs_sg_mean_bf.pdf", width = width, height = height )
par( mar = mar )
plot( log10( D$CC_fix_add[w] ), log10( D$SG_mean_bf[w] ), col = D$maf_cut[w], pch = pch[w],
xlab = "log10( additive c/c BF )\n(logit model)",
ylab = "log10( model-averaged BF )\n(general multinomial model)"
)
abline( a = 0, b = 1, col = 'red' )
legend( "bottomright", col = 1:length( levels( D$maf_cut )), legend = levels( D$maf_cut ), pch = 19 )
legend( "topright", legend = model.names, pch = shapes )
dev.off()

###############################################################
# ROC curve
w = which( !is.na( D$SG_mean_bf ) & D$SG_mean_bf != 0 )
D2 = D[w, ]
ROC.data = data.frame()
D2 = D2[ order( D2$CC_mean_bf, decreasing = T ), ]
ROC.data = rbind(
	ROC.data,
	data.frame(
		method = "additive, case/control",
		true.positive = cumsum( !1:nrow(D2) %in% grep( "null", D2$rsid )),
		false.positive = cumsum( 1:nrow(D2) %in% grep( "null", D2$rsid ))
	)
)
D2 = D2[ order( D2$SA_mean_bf, decreasing = T ), ]
ROC.data = rbind(
	ROC.data,
	data.frame(
		method = "additive, multinomial mean BF",
		true.positive = cumsum( !1:nrow(D2) %in% grep( "null", D2$rsid )),
		false.positive = cumsum( 1:nrow(D2) %in% grep( "null", D2$rsid ))
	)
)
D2 = D2[ order( D2$SG_mean_bf, decreasing = T ), ]
ROC.data = rbind(
	ROC.data,
	data.frame(
		method = "general, multinomial mean BF",
		true.positive = cumsum( !1:nrow(D2) %in% grep( "null", D2$rsid )),
		false.positive = cumsum( 1:nrow(D2) %in% grep( "null", D2$rsid ))
	)
)

pdf( file = "images/power.pdf", width = 8, height = 5 )
( ggplot( data = ROC.data ) + geom_line( aes( x = false.positive, y = true.positive, col = method ))
+ xlab( "Number of false positives" )
+ ylab( "Number of true positives" )
+ theme_bw( 16 )
)
dev.off()



###############################################################
# plot effect sizes
CC = dbGetQuery( db, "SELECT * FROM CCView" )
SA = dbGetQuery( db, "SELECT * FROM SubAddView" )
count = 1
for( term in c( "add", "het" )) {
	count = count + 1
	for( case in c( "case1", "case2", "case3" )) {
		effect.size.name = sprintf( "%s:%s", case, term )
		estimate.name = sprintf( "test:beta_.*:%s/sub=%s", term, case )
		if( count == 2 ) {
			plot( E[, grep( estimate.name, colnames( E ) ) ], effect.sizes[, effect.size.name ], xlim = c( -0.5, 0.5 ), ylim = c( -0.5, 0.5 )  )
		} else {
			points( E[, grep( estimate.name, colnames( E ) ) ], effect.sizes[, effect.size.name ] )
		}
		count = count + 1
	}
}

