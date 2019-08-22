\name{stx2sT}
\alias{stx2sT}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ Convert a Dendrometered Tree to a \dQuote{\pkg{sampSurf}}
\code{standingTree} Object }

\description{

  This routine will convert a tree that has been processed through
  \code{\link{SplineVolume}} into a \pkg{sampSurf}
  \dQuote{\code{\linkS4class{standingTree}}} object for use in
  simulations.
  
}

\usage{
stx2sT(sv.list,
       treeNo = 1,
       centerOffset = c(x = 0, y = 0),
       units.in = "English",
       units.out = "metric",
       prefix = "odr77:tree",
       runQuiet = TRUE,
       ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  
  \item{sv.list}{ An object generated from \code{SplineVolume}.
}

  \item{treeNo}{ The tree number in the \code{sv.list} to process.
}

  \item{centerOffset}{ This is used to locate the
	\dQuote{\code{\linkS4class{standingTree}}} object.  This must be
	specified here if coordinates are desired in order for all the other
	spatial points to be set set up correctly relative to the tree's
	center. See \code{\link{standingTree-methods}} for details. }

  \item{units.in}{ The units that the dendrometry was take in, this
                should normally be \dQuote{English}; however, it must
                conform to the two potential choices in \pkg{sampSurf}
                of either \dQuote{English} or \dQuote{metric}.  }
			  
  \item{units.out}{ The units for output, which must conform to the two
                 potential choices in \pkg{sampSurf} of either
                 \dQuote{English} or \dQuote{metric}.
}

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
	differentiation and future reference.  }

  \item{runQuiet}{ \code{TRUE}: no feedback; \code{FALSE}: a little
	information on the tree volume comparisons.
}

  \item{\dots}{ Gobbled. }

} %\argument section

\details{

 This routine uses the taper measurements taken by dendrometry to create
 a \pkg{sampSurf} "\code{\linkS4class{standingTree}}" object. It does
 this using the "standingTree" method that takes a signature first
 argument of a \code{\link{data.frame}}, representing the taper. The
 object constructor will then create the "standingTree's" taper slot
 from this information, fit a spline to it to compute volume, and use
 the \code{centerOffset} argument to assign a spatial location to the
 tree's pith. As noted above, it is important to assign the spatial
 location at the creation of the "standingTree" object so that all
 associated spatial polygon information (like the tree's taper polygon)
 will be accurately registered to the tree's location in the tract. For
 more information, please see the \pkg{sampSurf} documentation as noted
 below.
  
}

\value{ A \code{list} invisibly with\ldots

  \item{segs}{ The segment data frame for the chosen tree.}

  \item{tree}{ The same for the tree-level data frame (one record).}

  \item{stree}{ The "standingTree" object.}

  \item{useDBT}{ Whether double bark thickness was used (\code{TRUE}) or
	not (\code{FALSE}).}

  \item{volumes.are}{ Redundant with \code{useDBT}; whether volumes are
	inside or outside bark.}

  \item{aveVol}{ The average of the three methods of volume
  determination in SplineVolume. This can be compared against the volume
  from the spline fit used in creation of the "standingTree" object. }
  
}

\references{

\dQuote{The Stem Class} vignette in th the \pkg{sampSurf} package.

}

\author{
Jeffrey H. Gove
}

%% \note{
%% %%  ~~further notes~~
%%}

\seealso{
The \pkg{sampSurf} "\code{\linkS4class{standingTree}}" class, the
\code{\link{standingTree}} generic and associated methods.

In this package, this routine is called directly by
\code{\link{stx2sTs}}, which will build a
"\code{\linkS4class{standingTrees}}" container object from a set of
dendrometered trees processed by \code{SplineVolume}; these could, for
example, consist of the full set passed to this routine in the
\code{sv.list} argument above.

}

\examples{
#
# first process the raw data, then create a single
# standingTree object from tree #11...
#
\dontrun{
hatch16 = SplineVolume('Hatch Plot 16.stx')
h16.t11 = stx2sT(hatch16, 11, units.out = 'English', runQuiet = FALSE)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
