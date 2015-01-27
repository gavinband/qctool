\name{bgen.read_snp_probability_data}
\alias{bgen.read_snp_probability_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bgen.read_snp_probability_data(header, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{
%%     ~~Describe \code{header} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (header, data) 
{
    con = header$connection
    endian = header$endian
    if (header$flags == 0) {
        data$genotypes = readBin(con, integer(), size = 2, n = (data$number_of_samples * 
            3), endian = endian)
    }
    else {
        library(Rcompression)
        compressed_size = readBin(con, integer(), size = 4, n = 1, 
            endian = endian)
        compressed_genotypes = readBin(con, "raw", n = compressed_size, 
            endian = endian)
        uncompressed_genotypes = uncompress(content = compressed_genotypes, 
            size = data$number_of_samples * 6, asText = FALSE)
        data$genotypes = readBin(uncompressed_genotypes, integer(), 
            size = 2, n = data$number_of_samples * 3, signed = TRUE, 
            endian = endian)
        data$genotypes = matrix(data = data$genotypes/10000, 
            ncol = 3, nrow = data$number_of_samples, byrow = T)
    }
    data$ploidy = 2
    return(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
