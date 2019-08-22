\name{stx2sTs}
\alias{stx2sTs}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{

  Convert a List of Dendrometered Trees to a \dQuote{\pkg{sampSurf}}
  \code{standingTrees} Object

}

\description{

  This routine will convert the full output from a run of
   \code{\link{SplineVolume}} to a list of \pkg{sampSurf}
   \dQuote{\code{\linkS4class{standingTree}}} objects and finally to a
   \dQuote{\code{\linkS4class{standingTrees}}} container object. Note
   that this automatically generates tree coordinates within a
   \dQuote{\code{\linkS4class{bufferedTract}}} buffer.

}

\usage{
stx2sTs(sv.list,
        buffTr,
        units.in = "English",
        units.out = "metric",
        prefix = "odr77:tree",
        inhibitDist = 3,
        runQuiet = TRUE,
        ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sv.list}{ An object generated from \code{SplineVolume}. }
  
  \item{buffTr}{ A \pkg{sampSurf}
	\dQuote{\code{\linkS4class{bufferedTract}}} object. Tree coordinates
	will be drawn from the interior of the buffer.
}

  \item{units.in}{ The units that the dendrometry was take in, this
                should normally be \dQuote{English}; however, it must
                conform to the two potential choices in \pkg{sampSurf}
                of either \dQuote{English} or \dQuote{metric}.  }

\item{units.out}{ The units for output, which must conform to the two
                 potential choices in \pkg{sampSurf} of either
                 \dQuote{English} or \dQuote{metric}.  }
			   
  \item{prefix}{ A character identifier prefix for the tree id used in
	the \code{standingTree} object. It is important that this have some
	unique identification tied to the full set of trees dendrometered at
	a site as passed in \code{sv.list}. The reason is because most sites
	will probably begin tree dendrometry enumeration with tree number
	one, and continue from there; and \code{STX} (and hence
	\code{SplineVolume}) does not allow for much of an option in tree
	numbering beyond this. Importantly, if trees from a number of plots
	or sites are to be concatenated into one set of
	\dQuote{\code{\linkS4class{standingTrees}}} objects, one will be
	unable to create such an object without each tree having a unique
	identifier. This prefix allows such an ID by site for
	differentiation and future reference. }

  \item{inhibitDist}{ The spatial process inhibition distance in
	\code{units.out} units. (See below for more details.)
}

  \item{runQuiet}{ \code{TRUE}: no feedback; \code{FALSE}: a little
	information on the tree volume comparisons. }

  \item{\dots}{ Gobbled. }

} %\arguments

\details{

   This routine uses the \code{\link{stx2sT}} function to actually
   create individual \dQuote{\code{standingTree}} objects from the
   entire population of dendrometered trees in the \code{sv.list}
   argument. Individual tree coordinates are assigned for each tree and
   will fall within the \dQuote{\code{bufferedTract}} buffer's bounding
   box. The actual spatial coordinates are drawn from an inhibition
   process using the \code{\link{SSI}} function in the \pkg{spatial}
   package.

   The reason for passing back a \code{\link{list}} of
   \dQuote{\code{standingTree}} objects as well as the
   \dQuote{\code{standingTrees}} container object is that it facilitates
   joining one or more populations of \dQuote{\code{standingTree}}
   objects by simply concatenating the associated lists. One can not do
   this with \dQuote{\code{standingTrees}} container objects without
   running into problems with the bounding box, &., which will have to
   be recreated.

}

\value{ A \code{list} invisibly with\ldots

  \item{strees.list}{ A list of the individual
	\dQuote{\code{standingTree}} objects. }

  \item{strees}{ The same trees in a \dQuote{\code{standingTrees}}
	container object.}

  \item{prefix}{ The tree ID prefix used for the group of trees passed
	in \code{sv.list} as described above in argument \code{prefix}. }
  
}

\references{

\dQuote{The Stem Class} vignette in th the \pkg{sampSurf} package.

}

\author{ Jeffrey H. Gove }

%% \note{
%% %%  ~~further notes~~
%% }


\seealso{
  
The \pkg{sampSurf} "\code{\linkS4class{standingTrees}}" class, the
\code{\link{standingTrees}} generic and associated methods.

In this package, this routine directly calls \code{\link{stx2sT}}, which
will build an individual "\code{\linkS4class{standingTree}}" object for
each tree in the set of dendrometered trees processed by
\code{SplineVolume}, consisting of the full set passed to this routine
in the \code{sv.list} argument above.

}

\examples{
#
# first process the raw data, then create the standingTrees
# container of standingTree objects, assuming that a 'metric'
# bufferedTract object exists in the "btr" object below...
#
\dontrun{
hatch16 = SplineVolume('Hatch Plot 16.stx')
hp16.sts = stx2sTs(hatch16, btr, units.out = 'metric',
                   prefix = 'hp16:tree', runQuiet = FALSE)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
