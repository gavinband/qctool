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
    if( theta == 0 || theta == 1 ) {
        return( 1 )
    }
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
    if( theta == 0 || theta == 1 ) {
        return( 1 )
    }
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
    if( theta.hat == 0 || theta.hat == 1 ) {
        return( 1 )
    }

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

load.snp.stats <- function( filename, info_filename = NULL ) {
    X = read.table( filename )
    X$I = rowSums( X[, c( "AA", "AB", "BB" )])

    if( !is.null( info_filename )) {
        info = read.table( info_filename, as.is = T, header = F )
        colnames( info )[c(1:3,7)] = c( "SNPID", "rsid", "position", "type" )
        stopifnot( X$rsid == info$rsid )
        X$type = info$type
        X$type = factor( X$type, levels = c( 1, 2 ))
        levels( X$type ) = c( "imputed", "genotyped" )
    } else {
        X$type = factor( "genotyped", levels = c( "imputed", "genotyped" ))
    }

    X = melt( X, measure.vars = c("information", "jonathans_information" ))
    X$variable = as.character( X$variable )
    X$variable[ which( X$variable == "information" )] = "info.variant1"
    X$variable[ which( X$variable == "jonathans_information" )] = "info.variant2"
    return( X )
}

library( ggplot2 )

genotyped = load.snp.stats( "PS_01_illumina.snp-stats" )
imputed = load.snp.stats( "PS.illumina.22.80.500.imputed.gz.snp-stats", "PS.illumina.22.80.500.info" )

make.plot <- function( data ) {
    #X11( type = "Xlib" )
    p = ggplot( data = data )
    p = p + layer( mapping = aes( x = old_information, y = value, colour = MAF, shape = type, solid = type ), geom = "point" )
    p = p + scale_shape_manual( value =  c( 20, 3 ) )
    #p = p + scale_colour_manual( values = c( "red", "lightblue" ) )
    #p = p + stat_abline( xintercept = 0, slope = 1, colour = "red" )
    p = p + facet_wrap( ~ variable )
    p = p + xlab( "Original information measure" )
    p = p + ylab( "New information measure" )
    p = p + scale_colour_continuous( high = "black", low = "red" )
    return( p )
}
p = make.plot( imputed )
p = p + opts( title = "(PS chr22 imputed / genotyped data) info variant versus original" )
#print(p)
print(p + scale_x_continuous( limits = c( 0.9, 1.0 ) ) + scale_y_continuous( limits = c( 0.9, 1.0 ) ) )
ggsave( p, file = "PS.illumina.22.80.500.imputed_info_comparison.png", dpi = 72, width = 16, height = 8 )
ggsave( p + scale_x_continuous( limits = c( 0.9, 1.0 ) ) + scale_y_continuous( limits = c( 0.9, 1.0 ) ), file = "PS.illumina.22.80.500.imputed_info_comparison_zoom.png", dpi = 72, width = 16, height = 8 )

q = make.plot( genotyped )
q = q + opts( title = "(PS chr1 genotyped data) info variant versus original" )
print( q + scale_x_continuous( limits = c( 0.9, 1.0 ) ) + scale_y_continuous( limits = c( 0.9, 1.0 ) ) )
ggsave( q + scale_x_continuous( limits = c( 0.0, 1.0 ) ) + scale_y_continuous( limits = c( 0.0, 1.0 ) ), file = "PS_01_info_comparison.png", dpi = 72, width = 16, height = 8 )
ggsave( q + scale_x_continuous( limits = c( 0.9, 1.0 ) ) + scale_y_continuous( limits = c( 0.9, 1.0 ) ), file = "PS_01_info_comparison_zoom.png", dpi = 72, width = 16, height = 8 )

genotyped$missing_category = cut( genotyped$I, breaks = c( 0, 1000, 2000, 2500, 2600, 2610, 2620, 2622 ), labels = c( "I >= 0", "I >= 1000", "I >= 2000", "I >= 2500", "I >= 2600", "I >= 2610", "I >= 2620" ))
genotyped$MAF_category = cut( genotyped$MAF, breaks = c( -1, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5 ), labels = c( "MAF < 0.01", "MAF < 0.1", "MAF < 0.2", "MAF < 0.3", "MAF < 0.4", "MAF < 0.5" ))

p = ggplot( data = genotyped[ which( !is.na( genotyped$MAF )), ] )
p = p + layer( mapping = aes( x = old_information, y = value, colour = -log10( missing ), shape = MAF_category ), geom = "point" )
#p = p + scale_colour_manual( values = c( "red", "lightblue" ) )
#p = p + stat_abline( xintercept = 0, slope = 1, colour = "red" )
p = p + facet_wrap( ~ variable )
p = p + xlab( "Original information measure" )
p = p + ylab( "New information measure" )
p = p + scale_colour_continuous( high = "black", low = "red" )
print(p + scale_x_continuous( limits = c( 0.0, 1.0 ) ) + scale_y_continuous( limits = c( 0.0, 1.0 ) ))
ggsave( p + scale_x_continuous( limits = c( 0.0, 1.0 ) ) + scale_y_continuous( limits = c( 0.0, 1.0 ) ), file = "PS_01_info_comparison_coloured_by_missingness.png", dpi = 72, width = 16, height = 16 )
ggsave( p + scale_x_continuous( limits = c( 0.9, 1.0 ) ) + scale_y_continuous( limits = c( 0.9, 1.0 ) ), file = "PS_01_info_comparison_coloured_by_missingness_zoom.png", dpi = 72, width = 16, height = 16 )

p = ggplot( data = genotyped[ which( !is.na( genotyped$MAF )), ] )
p = p + layer( mapping = aes( x = old_information, y = value, colour = -log10( missing ) ), geom = "point" )
#p = p + scale_shape_manual( value =  c( 20, 3 ) )
#p = p + scale_colour_manual( values = c( "red", "lightblue" ) )
#p = p + stat_abline( xintercept = 0, slope = 1, colour = "red" )
p = p + facet_grid( MAF_category ~ variable )
p = p + xlab( "Original information measure" )
p = p + ylab( "New information measure" )
p = p + scale_colour_continuous( high = "black", low = "red" )
p = p + opts( title = "(PS chr1 genotyped data) Info variants vs original,facetted by MAF")
print(p + scale_x_continuous( limits = c( 0.0, 1.0 ) ) + scale_y_continuous( limits = c( 0.0, 1.0 ) ) )
ggsave( p + scale_x_continuous( limits = c( 0.0, 1.0 ) ) + scale_y_continuous( limits = c( 0.0, 1.0 ) ), file = "PS_01_info_by_MAF.png", dpi = 72, width = 16, height = 24 )
ggsave( p + scale_x_continuous( limits = c( 0.9, 1.0 ) ) + scale_y_continuous( limits = c( 0.9, 1.0 ) ), file = "PS_01_info_by_MAF_zoom.png", dpi = 72, width = 16, height = 24 )

snp.stats = read.table( "PS_01_illumina.snp-stats" )
snp.stats$MAF_category = cut( snp.stats$MAF, breaks = c( -1, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5 ), labels = c( "MAF < 0.01", "MAF < 0.1", "MAF < 0.2", "MAF < 0.3", "MAF < 0.4", "MAF < 0.5" ))
snp.stats$I = rowSums( snp.stats[, c( "AA", "AB", "BB" )] )
snp.stats$I_category = cut( snp.stats$I, breaks = c( 0, 1000, 2000, 2500, 2600, 2610, 2620, 2625 ), labels = c( "I < 1000", "I < 2000", "I < 2500", "I < 2600", "I < 2610", "I < 2620", "I <= 2622" ))

p = ggplot( data = snp.stats[ which( !is.na( snp.stats$MAF )), ] )
p = p + layer( mapping = aes( x = jonathans_information, y = information, colour = -log10( missing ), shape = MAF_category ), geom = "point" )
#p = p + scale_shape_manual( value =  c( 20, 3 ) )
#p = p + scale_colour_manual( values = c( "red", "lightblue" ) )
#p = p + stat_abline( xintercept = 0, slope = 1, colour = "red" )
p = p + xlab( "info variant2" )
p = p + ylab( "info variant1" )
#p = p + facet_wrap( ~MAF_category )

p = p + scale_colour_continuous( high = "black", low = "red" )
p = p + opts( title = "(PS chr01 genotyped data).  Info variant 1 vs variant 2.")
print(p)
ggsave( p, file = "PS_01_info_variant1_vs_variant2.png", dpi = 72, width = 8, height = 8 )

which( genotyped$rsid %in% imputed$rsid )

##########################


genotypes = read.table( gzfile( "PS_01_illumina.top.gen.gz" ), hea=F, as.is = T )
colnames( genotypes )[1:5] = c( "SNPID", "rsid", "pos", "allele1", "allele2" )

genotypes = read.table( "PS_master_01_illumina.top.gen", hea=F, as.is = T )
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
    strange = NA,
    original.flipped = NA,
    variant1.flipped = NA,
    variant2.flipped = NA,
    strange.flipped = NA
)

for( i in 1:100) {
    if( i %% 10 == 0 ) {
        cat( i, 'of', nrow(D), '...' )
    }
    G = matrix( data = as.numeric( genotypes[i, 6:ncol( genotypes) ]), ncol = 3, byrow = T )
    G.flipped = G[, c(3,2,1)]

    D$AA[i] = sum( G[, 1] )
    D$AB[i] = sum( G[, 2] )
    D$BB[i] = sum( G[, 3] )
    D$missing[i] = ( nrow(G) - sum( G[, 1:3] )) / nrow( G )
    D$original[i] = info.original( G )
    D$variant1[i] = info.variant1( G )
    D$variant2[i] = info.variant2( G )
    D$strange[i] = info.strange( G )

    D$original.flipped[i] = info.original( G.flipped )
    D$variant1.flipped[i] = info.variant1( G.flipped )
    D$variant2.flipped[i] = info.variant2( G.flipped )
    D$strange.flipped[i] = info.strange( G.flipped )
}
cat( i, 'of', nrow(D), '\n' )

which( abs( D$original - D$original.flipped ) > 0.0001 )
which( abs( D$variant1 - D$variant1.flipped ) > 0.0001 )
which( abs( D$variant2 - D$variant2.flipped ) > 0.0001 )

original.terms = sapply( 1:10, function( i ) {
    G = matrix( data = as.numeric( genotypes[i, 6:ncol( genotypes) ]), ncol = 3, byrow = T )
    return( info.original.terms( G ) )
})

variant1.terms = sapply( 1:10, function( i ) {
    G = matrix( data = as.numeric( genotypes[i, 6:ncol( genotypes) ]), ncol = 3, byrow = T )
    return( info.variant1.terms( G ) )
})

variant2.terms = sapply( 1:10, function( i ) {
    G = matrix( data = as.numeric( genotypes[i, 6:ncol( genotypes) ]), ncol = 3, byrow = T )
    return( info.variant2.terms( G ) )
})

original.terms
variant1.terms
variant2.terms