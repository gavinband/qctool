\name{bgen.read_header}
\alias{bgen.read_header}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bgen.read_header(con, endian = "little")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{con}{
%%     ~~Describe \code{con} here~~
}
  \item{endian}{
%%     ~~Describe \code{endian} here~~
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
function (con, endian = "little") 
{
    if (class(con) == "character") {
        con = file(con, open = "rb")
    }
    data = as.list(readBin(con, integer(), n = 5, size = 4, endian = endian))
    names(data) = c("offset", "header_length", "number_of_snps", 
        "number_of_samples", "reserved")
    if (data$header_length > 20) {
        data[["free_data"]] = readChar(con, nchars = (data$header_length - 
            20))
    }
    else {
        data[["free_data"]] = ""
    }
    data[["flags"]] = readBin(con, integer(), n = 1, size = 4, 
        endian = endian)
    data[["endian"]] = endian
    data[["type"]] = "bgen"
    data$connection = con
    return(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
