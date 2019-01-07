logistic <- function(x) {
	exp(x) / (1+exp(x))
}

dlogf <- function( x, nu1, nu2 ) {
	a = nu1/2
	b = nu2/2
	l = logistic(x)
	(l^a * (1 - l)^b)/beta(a,b)
}

plogf <- function( x, nu1, nu2 ) {
	a = nu1/2
	b = nu2/2
	pbeta( logistic(x), shape1 = a+1, shape2 = b+1 )
}

dlogf <- function( x, nu1, nu2 ) {
	a = nu1/2
	b = nu2/2
	exp(-b * x ) / (( 1 + exp(-x))^(a+b) * beta(a,b))
}

dnmixture <- function( x, sds = c( 0.05, 0.1, 0.2, 0.4 ), alphas = rep( 0.25, 4 ) ) {
	result = numeric(length(x))
	for( i in 1:length( sds )) {
		result = result + alphas[i] * dnorm( x, sd = sds[i] )
	}
	return( result )
}

hptest = "../../../build/release/hptest_v2.1-dev" 
system(
	sprintf( "%s -predictor host.vcf -outcome parasite.vcf -s samples.sample -o /tmp/output.txt", hptest )
)
X = read.table( "/tmp/output.txt", hea=T, as.is=T, comment = '#', check.names = F )
host.spec = read.table( "host_spec.txt", hea=T, as.is= T)
interactors = host.spec[ !is.na( host.spec$interactors ), ]
interactors$index = match(
	paste( interactors$variant, interactors$interactors, sep = ':'),
	paste( X$`predictor:rsid`, X$`outcome:rsid`, sep = ':' )
)

X$col = 'black'
X$col[ interactors$index ] = 'red'
X$true_beta = 0
X$true_beta[ interactors$index ] = interactors$beta
X$fit_beta = X$`model_1:beta_1:add/gp=1`
X$fit_beta_lower = X$fit_beta - 1.96 * X$`model_1:se_1`
X$fit_beta_upper = X$fit_beta + 1.96 * X$`model_1:se_1`

X = X[ order( X$`model_1:log10_bf`, decreasing = T ), ]

par( mfrow = c( 2, 3 ))
plot( X$true_beta, X$fit_beta, pch = 19, cex = 0.75, col = X$col )
abline( a = 0, b = 1, col = 'red' )
segments( x0 = X$true_beta, x1 = X$true_beta, y0 = X$fit_beta_lower, y1 = X$fit_beta_upper, col = 'grey' )
legend(
	"bottomright",
	legend = c( "beta = 0", "|beta| != 0" ),
	pch = 19,
	col = c( "black", "red" )
)

## ROC curve
plot(
	cumsum( X$true_beta == 0 ) / nrow(X), cumsum( X$true_beta != 0 ) / sum( X$true_beta != 0.0 ),
	type = 'l',
	xlab = "False positive rate",
	ylab = "True positive rate"
)
points(
	cumsum( X$true_beta == 0 ) / nrow(X), cumsum( abs( X$true_beta ) > 0.2 ) / sum( abs( X$true_beta ) > 0.2 ),
	type = 'l',
	col = 'grey'
)
points(
	cumsum( X$true_beta == 0 ) / nrow(X), cumsum( abs( X$true_beta ) > 0.5 ) / sum( abs( X$true_beta ) > 0.5 ),
	type = 'l',
	col = 'blue'
)
legend(
	"bottomright",
	legend = c( "beta!=0", "|beta| > 0.2", "|beta| > 0.5" ),
	pch = 19,
	col = c( "black", "grey", "blue" )
)
abline( a = 0, b = 1, col = 'red' )

### Compute rates
quantiles = c( 0.01, 0.025, 0.05 )
rates = data.frame(
	quantile = quantiles,
	rate = sapply(
		quantiles,
		function(q) {
			length( which( ( abs( X$fit_beta - X$true_beta ) / X$`model_1:se_1` ) > qnorm( 1 - (q/2) ))) / nrow(X)
		}
	)
)

x = seq( from = -5, to = 5, by = 0.01 )
plot( x, dnmixture( x ), type = 'l', lty = 2, col = 'black' )
points( x, dnmixture( x, sds = c( 0.2, 0.4, 0.6, 0.8 ), alphas = rep( 0.25, 4 ) ), type = 'l', lty = 2, col = 'grey' )
points( x, dnmixture( x, sds = c( 0.2 ), alphas = 1 ), type = 'l', lty = 2, col = 'green' )

points( x, dlogf( x, 16, 16 ), type = 'l', col = 'black' )
points( x, dlogf( x, 32, 32 ), type = 'l', col = 'red' )
points( x, dlogf( x, 64, 64 ), type = 'l', col = 'blue' )

plot( x, dlogf( x, 4, 4 ), type = 'l', col = 'blue' )
points( x, dlogf( x, 1, 1 ), type = 'l', col = 'black' )
points( x, dlogf( x, 2, 2 ), type = 'l', col = 'red' )

hist( X$`model_1:iterations`, 100 )
