\name{hardcopyLattice}
\alias{hardcopyLattice}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Prints a trellis/lattice object to a file device }

\description{
  This routine takes a lattice/trellis object and prints it to one of
  the following devices: jpg, tiff, eps, or pdf.
}

\usage{
hardcopyLattice(plt = NULL,
                fileName = '',
                filePath = file.path(getwd(),'figures'),
                eps2pdf = TRUE,
                family = 'Helvetica',
                pointsize.text = 16,
                runQuiet = FALSE,
                columnsPerRow = 2,
                pdf.height = 8,
                pdf.width = 8
               )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plt}{ a list of lattice plot object of class 'trellis' OR a single
              object of class trellis that will be wrapped into a list for
              processing }
  \item{fileName}{ character file name to send the output to; the
	extension determines the type of device to open }
  \item{filePath}{ the complete path to fileName }
  \item{eps2pdf}{ \code{TRUE}: convert eps to pdf; \code{FALSE}: leave
	as eps. \emph{This option requires the \code{epstopdf} utility to be
	installed on your system!} }
  \item{family}{ font family }
  \item{pointsize.text}{ overall pointsize, not real discerning at this
	point (pun intended) }
  \item{runQuiet}{ \code{FALSE}: let user know the file is being written;
    \code{TRUE}: run silently}
  \item{columnsPerRow}{ number of columns for each row in output matrix }
  \item{pdf.height}{ in inches; note that this can affect the way panels are
                    presented. For example, if numPlts=4 and pdf.height=8,
                    with 2 time periods (2 panels), they are side-by-side;
                    but if pdf.height=10, they are one atop the other in each
                    of the 4 split displays.}
  \item{pdf.width}{ width of a pdf output; making this larger that pdf.height
                    is one way to get "landscape" in pdf output }
}

\details{
  This routine does the same thing for each device type (except when
  \code{eps2pdf} is in effect). After rescaling the font in the object,
  it simply opens an output device of the correct type and prints the
  lattice object to that device, then closes it.

  It now allows one to pass a list of trellis plot objects for printing
  on a single page as arranged via \code{columnsPerRow}.
}

\value{
  Nothing is returned.
}
%\references{ ~put references to the literature/web site here ~ }

\author{ Jeffrey Gove }
%\note{ ~~further notes~~ 
%}

\seealso{ \code{\link{print.trellis}}, \code{\link{device}}, \code{\link{lattice}}  }

\examples{
\dontrun{
  Depth <- equal.count(quakes$depth, number=8, overlap=.1)
  j = xyplot(lat ~ long | Depth, data = quakes)
  hardcopyLattice(j, 'j.pdf')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ device }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
