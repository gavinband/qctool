data = matrix( NA, nrow = 10, ncol = 3 )
data[,1] = runif(10)
data[,2] = runif(10) * 1-data[,1]
data[,3] = 1-(data[,1]+data[,2])

generate.data <- function( sexes ) {
  N = length( sexes )
  result = matrix( NA, nrow = N, ncol = 3 )
  result[,1] = runif( N )

  wM = which( sexes == 'M' )
  if( length(wM) > 0 ) {
    result[wM,2] = 0
    result[wM,3] = 1 - result[wM,1]
  }

  wF = which( sexes == 'F' )
  if( length( wF ) > 0 ) {
      result[wF,2] = (1-result[wF,1]) * runif( length( wF ))
      result[wF,3] = 1 - result[wF,1] - result[wF,2]
  }
  return( result )
}

freq <- function( data, gender ) {
  wM = which( gender == 'M' )
  wF = which( gender == 'F' )
  alleleBCount = 0
  totalCount = 0
  if( length( wM ) > 0 ) {
    # Males coded like homozygotes
    alleleBCount = alleleBCount + sum( data[wM,3] )
    totalCount = totalCount + sum( data[wM,] )
  }
  if( length( wF ) > 0 ) {
    alleleBCount = alleleBCount + sum( data[wF,2] + 2*data[wF,3] )
    totalCount = totalCount + 2 * sum( data[wF,] )
  }
  return( alleleBCount / totalCount )
}

info <- function( data, gender ) {
  wM = which( gender == 'M' )
  wF = which( gender == 'F' )
  f = freq( data, gender )
  E = rep( 0, length( gender ))
  D = rep( 0, length( gender ))
  if( length( wM ) > 0 ) {
    mean.genotype = data[wM,3]
    mean.squared.genotype = data[wM,3]
    E[wM] = mean.squared.genotype - mean.genotype^2
    D[wM] = f*(1-f)
  }
  if( length( wF ) > 0 ) {
    mean.genotype = data[wF,2] + data[wF,3]*2
    mean.squared.genotype = (data[wF,2] * 1 + data[wF,3]*4)
    E[wF] = mean.squared.genotype - mean.genotype^2
    D[wF] = 2*f*(1-f)
  }
  wNotNA = which( !is.na( gender ))
  return( list(
    freq = f,
    E = E,
    D = D,
    info = 1 - (sum(E[wNotNA]/D[wNotNA]) / sum(data[wNotNA,]))
  ))
}

impute.info <- function( data ) {
  N = nrow( data )
  f = freq( data, rep('F', N) )
  E = rep( 0, N)
    mean.genotype = data[,2] + data[,3]*2
    mean.squared.genotype = (data[,2] * 1 + data[,3]*4)
    E = mean.squared.genotype - mean.genotype^2
    D = 2*f*(1-f)
  return( list(
    freq = f,
    E = E,
    D = D,
    info = 1 - (sum(E/D) / sum(data))
  ))
}

do.qctool <- function( data, sex, chromosome = 'X', qctool = "qctool_v2.0-beta6" ) {
  tmp1 = tempfile( fileext = '.gen' )
  tmp2 = tempfile( fileext = '.sample' )
  write.table(
    cbind(
      data.frame( chromosome = chromosome, SNPID = 'SNP1', rsid = 'rs1', position = 1, alleleA = 'A', alleleB = 'B' ),
      t( as.numeric( t(data)) )
    ),
    file = tmp1,
    col.names = F,
    row.names = F,
    quote = F
  )
  
  write( "ID_1 sex\n0 D", tmp2 )
  write.table(
    data.frame(
      ID = 1:nrow(data),
      sex = sex
    ),
    col.names = F,
    row.names = F, 
    quote = F,
    append = T,
    file = tmp2
  )
  tmp3 = tempfile( fileext = '.stats' )
  cmd = sprintf( '/Users/gav/Projects/Software/qctool/build/release/%s -g %s -s %s -snp-stats -osnp %s', qctool, tmp1, tmp2, tmp3 )
  system( cmd )
  result = read.table( tmp3, head = T, as.is=T )
  return( result )
}

