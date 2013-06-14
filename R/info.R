info.strange <- function( genotypes ) {
    I = sum( rowSums( genotypes ))
    N = nrow( genotypes )
    
    return( I / N ) ;
}

info.original <- function( genotypes ) {
    I = sum( rowSums( genotypes ))
    if( I == 0 ) {
        return( 0 ) ;
    }
    N = nrow( genotypes )
    f = genotypes[, 2] + 4 * genotypes[, 3]
    e = genotypes[, 2] + 2 * genotypes[, 3]
    theta = sum(e) / (2*I)
    ratio = sum( f-e^2 ) / ( 2*N*theta*(1-theta) )
    return( 1 - ratio )
}

info.original.terms <- function( genotypes ) {
    I = sum( rowSums( genotypes ))
    if( I == 0 ) {
        return( 0 ) ;
    }
    N = nrow( genotypes )
    f = genotypes[, 2] + 4 * genotypes[, 3]
    e = genotypes[, 2] + 2 * genotypes[, 3]
    theta = sum(e) / (2*I)
    ratio = sum( f-e^2 ) / ( 2*N*theta*(1-theta) )
    return( list( info = 1 - ratio, variance.terms = sum( f-e^2 ), denominator = ( 2*N*theta*(1-theta) ) ))
}

info.variant1 <- function( genotypes ) {
    I = sum( rowSums( genotypes ))
    if( I == 0 ) {
        return( 0 ) ;
    }
    N = nrow( genotypes )
    f = genotypes[, 2] + 4 * genotypes[, 3]
    e = genotypes[, 2] + 2 * genotypes[, 3]
    adj1 = ( 1.0 - rowSums( genotypes ) ) * e
    adj2 = ( 1.0 - rowSums( genotypes ) )^2

    theta = sum(e) / (2*I)
    ratio = sum( f-e^2 ) / ( 2*N*theta*(1-theta) )

    numerator = sum( f - e^2 ) + (( N - I ) * ( 2 * theta * ( 1 + theta ))) - 4 * theta * sum( adj1 ) - 4 * theta^2 * sum( adj2 )
    denominator = 2*N*theta*(1-theta)
    return( 1 - ( numerator / denominator ))
}

info.variant1.terms <- function( genotypes ) {
    I = sum( rowSums( genotypes ))
    if( I == 0 ) {
        return( 0 ) ;
    }
    N = nrow( genotypes )
    f = genotypes[, 2] + 4 * genotypes[, 3]
    e = genotypes[, 2] + 2 * genotypes[, 3]
    adj1 = ( 1.0 - rowSums( genotypes ) ) * e
    adj2 = ( 1.0 - rowSums( genotypes ) )^2

    theta = sum(e) / (2*I)
    ratio = sum( f-e^2 ) / ( 2*N*theta*(1-theta) )

    numerator1 = sum( f - e^2 )
    numerator2 = (( N - I ) * ( 2 * theta * ( 1 + theta ))) - 4 * theta * sum( adj1 ) - 4 * theta^2 * sum( adj2 )
    numerator = numerator1 + numerator2
    denominator = 2*N*theta*(1-theta)

    return( list( info = 1 - ( numerator / denominator ), variance.terms = numerator1, covariance.terms = numerator2, denominator = denominator ))
}

info.variant2 <- function( g ) {
    I = sum( g[, 1:3] )
    if( I == 0 ) {
        return( 0 ) ;
    }
    N = nrow( g )
    e = g[,2] + 2*g[,3]
    theta.hat = sum(e) / ( 2 * I )

    vF0 = sum(g[,1]*(1-g[,1]))
    vF1 = sum(g[,2]*(1-g[,2]))
    vF2 = sum(g[,3]*(1-g[,3]))

    c01 = -sum(g[,2]*g[,1])
    c02 = -sum(g[,1]*g[,3])
    c12 = -sum(g[,2]*g[,3])


    numerator1 = ( 4 * (theta.hat^2) * vF0 ) + ((1-(2*theta.hat))^2 * vF1) + (4*(1-theta.hat)^2 * vF2)

    numerator2 = - (4*theta.hat*(1-(2*theta.hat))*c01) + (4*(1-theta.hat)*(1-(2*theta.hat))*c12) - (8*theta.hat*(1-theta.hat)*c02)
    numerator = numerator1 + numerator2
    denominator = 2 * I * theta.hat * ( 1 - theta.hat )

    #cat( "theta_mle =", theta.hat, "v0 =", vF0, "v1 =", vF1, "v2 =", vF2, "c01 =", c01, "c02 =", c02, "c12 =", c12, ".\n" )
    #cat( "numerator1 =", numerator, "numerator2 =", numerator2, "numerator =", numerator, ".\n" )
    return( 1.0 - ( numerator / denominator ) )
}

info.variant2.terms <- function( g ) {
    I = sum( g[, 1:3] )
    if( I == 0 ) {
        return( 0 ) ;
    }
    N = nrow( g )
    e = g[,2] + 2*g[,3]
    theta.hat = sum(e) / ( 2 * I )

    vF0 = sum(g[,1]*(1-g[,1]))
    vF1 = sum(g[,2]*(1-g[,2]))
    vF2 = sum(g[,3]*(1-g[,3]))

    c01 = -sum(g[,2]*g[,1])
    c02 = -sum(g[,1]*g[,3])
    c12 = -sum(g[,2]*g[,3])


    numerator1 = ( 4 * (theta.hat^2) * vF0 ) + ((1-(2*theta.hat))^2 * vF1) + (4*(1-theta.hat)^2 * vF2)

    numerator2 = - (4*theta.hat*(1-(2*theta.hat))*c01) + (4*(1-theta.hat)*(1-(2*theta.hat))*c12) - (8*theta.hat*(1-theta.hat)*c02)
    numerator = numerator1 + numerator2
    denominator = 2 * I * theta.hat * ( 1 - theta.hat )
    return( list( info = ( 1 - ( numerator / denominator )), variance.terms = numerator1, covariance.terms = numerator2, denominator = denominator ))
}

X = read.table( "PS_01_illumina.top.snp-stats" )

genotypes = read.table( "PS_01_illumina.top.gen", hea=F, as.is = T )
colnames( genotypes )[1:5] = c( "SNPID", "rsid", "pos", "allele1", "allele2" )

D = data.frame(
    rsid = genotypes$rsid,
    pos = genotypes$pos,
    AA = NA,
    AB = NA,
    BB = NA,
    missing = NA,
    original = NA,
    variant1 = NA,
    variant2 = NA,
    strange = NA
)

for( i in 1:10) {
    G = matrix( data = as.numeric( genotypes[i, 6:ncol( genotypes) ]), ncol = 3, byrow = T )
    D$AA[i] = sum( G[, 1] )
    D$AB[i] = sum( G[, 2] )
    D$BB[i] = sum( G[, 3] )
    D$missing[i] = ( nrow(G) - sum( G[, 1:3] )) / nrow( G )
    D$original[i] = info.original( G )
    D$variant1[i] = info.variant1( G )
    D$variant2[i] = info.variant2( G )
    D$strange[i] = info.strange( G )
}

original.terms = sapply( c( 61, 540 ), function( i ) {
    G = matrix( data = as.numeric( genotypes[i, 6:ncol( genotypes) ]), ncol = 3, byrow = T )
    return( info.original.terms( G ) )
})

variant1.terms = sapply( c( 61, 540 ), function( i ) {
    G = matrix( data = as.numeric( genotypes[i, 6:ncol( genotypes) ]), ncol = 3, byrow = T )
    return( info.variant1.terms( G ) )
})

variant2.terms = sapply( c( 61, 540 ), function( i ) {
    G = matrix( data = as.numeric( genotypes[i, 6:ncol( genotypes) ]), ncol = 3, byrow = T )
    return( info.variant2.terms( G ) )
})

original.terms
variant1.terms
variant2.terms
