
// simulate 100 individuals with these frequencies

genotypes1 = data.frame(
    chromosome = "01",
    SNPID = sprintf( "SNP%d", 1:100 ),
    rsid = sprintf( "rs%d", 1:100 ),
    position = 1:100,
    allele1 = "A",
    allele1 = "G"
)

genotypes2 = data.frame(
    chromosome = "01",
    SNPID = sprintf( "SNP%d", 1:100 ),
    rsid = sprintf( "rs%d", 1:100 ),
    position = 1:100,
    allele1 = "A",
    allele1 = "G"
)

V = data.frame(
    index = 1:100,
    v0 = NA,
    v1 = NA,
    v2 = NA,
    v3 = NA
)

genotypes1[,7:306] = NA
genotypes2[,7:306] = NA
for( i in 1:nrow( genotypes1 )) {
    v = c( sort( runif( 3, min = 0, max = 1 ) ), 1 )
    V[i,2:5] = v
    for( j in 1:100 ) {
        genotypes = get_genotypes( v ) ;
        genotypes1[i, (7 + ( 3 * (j-1) )):(9 + ( 3 * (j-1) ))] = genotypes$first_snp
        genotypes2[i, (7 + ( 3 * (j-1) )):(9 + ( 3 * (j-1) ))] = genotypes$second_snp
    }
}

V$pi00 = V$v0
V$pi01 = V$v1 - V$v0
V$pi10 = V$v2 - V$v1
V$pi11 = V$v3 - V$v2

get_probs <- function( g ) {
    A = c( 0, 0, 0 )
    result[g+1] = 1
    return( result ) ;
}

get_genotypes <- function( v ) {
    p = runif( 2, min = 0, max = 1 ) ;

    first_snp = 0
    second_snp = 0

    for( i in 1:2 ) {
        if( p[i] < v[1] ) {
        
        } else if( p[i] < v[2] ) {
            second_snp = second_snp + 1
        } else if( p[i] < v[3] ) {
            first_snp = first_snp + 1
        } else {
            first_snp = first_snp + 1
            second_snp = second_snp + 1
        }
    }
    
    A = c( 0, 0, 0 )
    B = c( 0, 0, 0 )
    A[ first_snp + 1 ] = 1 ;
    B[ second_snp + 1 ] = 1 ;
    
    return( list( first_snp = A, second_snp = B ) ) ;
}
