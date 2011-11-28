\name{SplineVolume}
\alias{SplineVolume}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Translate an \acronym{STX} dendrometry file and calculate some volumes }

\description{ This routine can be used to read a file in \acronym{STX} dendrometry
  format and process the Barr and Stroud FP-12 or FP-15 measurements
  into corresponding segment diameters and lengths. More details are
  given below and in the Fortran and \R source code.
}

\usage{
SplineVolume(fileName = "odr_wp_1977.stx",
             filePath = getwd(),
             dataPath = "data",
             reportPath = "reports",
             splineMethods = c("natural","monoH"),
             ... )
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileName}{ An input file name with extension corresponding to an
	\acronym{STX} file on disk. }

  \item{filePath}{ The full file path to the current \R
    workspace, the default is \code{getwd()}. }

  \item{dataPath}{ The relative path from \code{filePath} to where the
	\acronym{STX} dendrometry data file resides; i.e, if the current \R workspace is
	in /Rwork, and the \acronym{STX} file is in /Rwork/data, then the default will
	work fine. The file should have the 10 control records (cards)
	before the first tree record as in a normal \acronym{STX} file. }

  \item{reportPath}{ The relative subdirectory where the report files
	should be written. }

  \item{splineMethods}{ One or more of the methods available in the
    \code{\link{spline}} and \code{\link{splinefun}} routines in the
    \code{utils} package. The defaults are the \sQuote{natural} and
    \sQuote{monoH.FC} spline methods.}

  \item{\dots}{ Just gobbled for now.}

}

\details{

This routine has its roots in some old Fortran code used to calculate
cubic, Hermite and Akima spline volumes (hence the name) in addition to
conic volumes for segments for dendrometry data. In the old Fortran
code, input data were in a form that Lew Grosenbaugh's \acronym{STX}
program would use. This is still the case, but the Fortran code has been
adapted to be used with the \R \code{SplineVolume} interface to read and
process the dendrometry data such that it is returned in data frames for
further use in \R. In addition, three files are written as output from a
run of \code{SplineVolume} (Fortran code) to the \code{reportPath}
subdirectory: (1) a tree file with one record per tree, (2) a segment
file with one record per measurement (dendrometer reading) for each
tree, and (3) a report file with a formatted tabular report.

The Akima and Hermite splines are not calculated in the current Fortran
version because the original code used licensed IMSL routines at the
time. They can be fitted and volumes derived from the results using the
appropriate \R functions based on the results returned from this
routine. A Hermite version (see \code{splineMethods} argument) is fitted
in \code{SplineVolume} and Akima can be fitted using the \pkg{akima}
package. Similarly, I have removed the \sQuote{natural cubic spline}
fits from the Fortran because of the source license (Numerical Recipes)
restrictions. A similar natural cubic spline (see \code{splineMethods}
argument) is fitted in the \R code.

This function is a wrapper for the Fortran code, which was origianlly
designed to be a stand-alone program.  At present, parsing the
\acronym{STX} file to determine the total number of segments is somewhat
of a pain. This is done tree-by-tree as the dendrometer measurements are
processed in the code. Because each tree can have variable numbers of
measurements, and because arrays are allocated and deallocated on a
tree-by-tree basis in the Fortran, when combined with \R's
\code{character} variable limitations with Fortran, it was too
problematic to pass back the actual arrays to \R. Rather the results are
written out to the files mentioned above in the Fortran, and then
\code{read.table} is used to recover these results in the \R code. This
does result in some loss of precision, but we are only talking about
diameters and lengths here, nothing that really requires full precision
so it should be enough until someone asks for more. More information on
the process is detailed in the \R and Fortran source files in the
package.

Note that the dendrometry measurements in the \acronym{STX} input file
are assumed to be in the form of direct measurements in the first two
segments as illustrated on page 31 in Space (1974); e.g., the first two
readings are\ldots

\describe{

\item{Reading 1}{ \acronym{FGRAD} = Stump diameter, and \acronym{SINELV} =
  height from stump to dbh, signalled with a negative (-999)
  \acronym{TGRAD}.}

\item{Reading 2}{ \acronym{FGRAD} = \acronym{DBH}, and \acronym{SINELV} =
  measured height from \acronym{DBH} to the first Barr & Stroud
  measurement, often something convenient like 4 feet above
  \acronym{DBH} (\acronym{TGRAD} = 0.0).}  }

The remainder of the measurements are assumed to have been taken with
the Barr & Stroud beyond that. The program should work if all
measurements have been taken using a Barr & Stroud (the routine that
processes measurements was taken directly from \acronym{STXMOD} and
this option is handled in the code, but has not been tested). It would be
simple enough to add other alternatives, like all direct measurements,
but it has not been done as all of our dendrometry has been taken with
the above model.

The code can not be guaranteed, unfortunately. A lot of it was borrowed
from \acronym{STX} (actually \acronym{STXMOD}), and that is not easy or
simple. There could be errors. Everything that has been tested so far
looks fine and the original versions were tesed against \acronym{STX}
for conversions of \acronym{B&S} readings.

Finally, please note that this is not a substitute for \acronym{STX} in
the sense of \acronym{3P} sampling as it is only used to process
measurements. 
  
}

\value{
  \describe{
	A list is returned invisibly as\ldots
	
  \item{\bold{z}}{A list returned from the Fortran program with self evident
	components as used in the call to the Fortran code.}
  
  \item{\bold{trees}}{ A data frame with the tree information, one record per
	tree.}

  \itemize{
	\item \emph{TreeNo}: The tree number	

	\item \emph{Spp}: Species

	\item \emph{dbhob}: Diameter at breast height outside bark.

	\item \emph{dbThick}: Double bark thinkness.

	\item \emph{isForked}: \code{TRUE}: the tree is forked;
	\code{FALSE}: not forked.

	\item \emph{nReadings}: Number of dendrometer + direct readings for
	the tree.

	\item \emph{totLen}: The total length of the tree including forks lengths.
	
	\item \emph{lenFlag}: \code{TRUE}: There was a problem with negative
	lengths in one of the segments measured, this is probably due to an
	error in the dendrometry reading; \code{FALSE}:
	everything is fine.
	
	\item \emph{conicVol}: Tree volumes using the formula for a
	conic frustum.

	\item \emph{natVol, monVol}: The tree volumes from the spline
	fits. The names of these fields depend on what is chosen as the
	splines to be fitted, the defaults are given above in the
	\code{splineMethods} argument. So be aware that the name of this
	column may be different for different data sets. The rule is to use
	the first 3 characters from the names and \code{splineMethods} and
	append \dQuote{Vol} for the name.

  }

  
  \item{\bold{segments}}{ A data frame with the segment information, one record
    per dendrometer reading (so \dQuote{readings} would probably have been
	a more appropriate name) within each tree. The columns in the data
	frame are the following\ldots }

  \itemize{
	\item \emph{TreeNo}: The tree number	

	\item \emph{Reading}: The dendrometer reading number; note that
	  begin at the top reading and decrease to the bottom reading within
	  a tree.

	\item \emph{Spp}: Species
	  
	\item \emph{dob}: Diameter outside bark
	
	\item \emph{dib}: Diameter inside bark--may be the same as
	  \code{dob} if no bark thinkness was recorded.
	
	\item \emph{prod}: The produce code for the segment as described
	below.
	
	\item \emph{tgrad}: The \sQuote{tgrad} readings

	\item \emph{fgrad}: The \sQuote{fgrad} readings

	\item \emph{sinelv}: The \sQuote{sinelv} readings

	\item \emph{instRange}: The distance from the instrument to the
	  reading on the tree.

	\item \emph{hgtHoriz}: The height of the readings above horizontal
	as determined by the instrument height.
	
	\item \emph{segLen}: The segment length between two dendrometer
	readings.
	
	\item \emph{height}: The cumulative height from the previous field;
	note that this will not equal the real height of the tree in the
	case of forked trees.

	\item \emph{conicVol}: Segment volumes using the formula for a
	conic frustum.

	\item \emph{natVol, monVol}: The segment volumes from the spline
	fits. The names of these fields depend on what is chosen as the
	splines to be fitted, the defaults are given above in the
	\code{splineMethods} argument. So be aware that the name of this
	column may be different for different data sets. The rule is to use
	the first 3 characters from the names and

  } %itemize
  
}

In the \code{segment} data frame, please note that the records
correspond to dendrometer readings from the bottom to top dendrometer
reading as recorded in the field. If the first reading is \code{i=1},
then the segment volumes, lengths, etc., are such that \code{i=2}
contains the segment information for readings \code{c(1,2)}; obviously
then reading \code{i+1} contains the segment information for readings
\code{c(i,i+1)}.

A final word concerning the spline volumes. The splines can be very
(too) flexible and give fits that are not always good depending on the
tree shape. One can see this by plotting the results in
\code{\link{taperPlot}}. So please be aware that the volume can be off
and that the conic versions are perhaps the best all around.

}

\references{
  
\describe{
%
\item{}{Grosenbaugh, L. R. 1974. \acronym{STX} 3-3-73: Tree Content and
  Value Estimation Using Various Sample Designs, Dendrometry Methods,
  and V-S-L Conversion Coefficients. U.S. Forest Service Research Paper
  SE-117. }
%
\item{}{Rennie, J. C. 1977. \acronym{STX}: A Computer System for
  Processing Timber Inventory Data. Unpublished manuscript.}
%
\item{}{Space, J. C. 1974. 3-P Forest Inventory: Design, Procedures,
  Data Processing. U.S. Forest Service State and Private Forestry,
  Southeastern Area, Atlanta, GA. }
%
 } 
}

\author{
Jeffrey Gove
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

%\examples{
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
