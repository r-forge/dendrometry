\name{splineSegVol}
\alias{splineSegVol}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Fit a spline and calculate segment volumes }

\description{ This function is used to fit a spline from one of the
  basic methods offered through the \code{splinefun} function in the
  \code{utils} package. It then calculates segment volumes from the
  spline fit. }

\usage{
splineSegVol(segments, method = "mono", ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{segments}{ \emph{One} tree's worth of segment (really readings)
	information from a \code{segments} data frame generated from
	\code{\link{SplineVolume}}. Note that this should be a data frame in
	the same format.}

  \item{method}{ One of the methods allowed in \code{\link{spline}} and
	\code{\link{splinefun}}.}

  \item{\dots}{ Any other arguments to be passed to \code{\link{splinefun}}.}
}

\details{ The routine was not really written to be used outside of
  \code{SplineVolume}, but it can be, and it is a simple way to add
  other spline fits to the data frame returned from that routine. 
}

\value{
The routine returns a vector of segment volumes. It has length one less
than the segment data frame (number of dendrometer readings) passed to
it. The order is from base segment to top segment, so it must be
reversed and padded with a zero element if one desires to add it to the
input data frame.
}

\author{
Jeffrey Gove
}


%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
#
# if treesegs is a data frame from SplineVolume...
#
\dontrun{
segs = treesegs[treesegs$TreeNo == 1,] #assume TreeNo 1 exists
vol = splineSegVol(segs, 'mono')
segs[,'monoVol'] = c(rev(vol), 0.0) #pad initial reading
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
