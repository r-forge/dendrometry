\name{Dendrometry-package}
\alias{Dendrometry-package}
\alias{Dendrometry}
\docType{package}

\title{
Tree Dendrometry Processing with Barr & Stroud Dendrometers
}

\description{
This package will process dendrometry files that are in STX format for
trees that have been dendrometered with the Barr and Stroud FP-12 or
FP-15. It will parse the STX files and determine overall tree and
segment information. It will not do any sampling analysis, use STX for
that. See the routines for more information.
}

\details{
\tabular{ll}{
Package: \tab Dendrometry\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2011-11-02\cr
License: \tab GPL (>= 3)\cr
LazyLoad: \tab yes\cr
}

The \code{\link{SplineVolume}} function is an interface to some
fortran code that reads the STX files and dis-assembles them into
something that can be used easily in \R.
}

\author{
Jeffrey Gove

Maintainer: jgove@fs.fed.us
}

\references{
STX, STXMOD source code.
}

\keyword{ package }

\seealso{
 \code{\link{SplineVolume}} 
}

\examples{
#
# see routines for examples
#
}
