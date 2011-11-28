taperPlot = function(ds,
                     TreeNos = NA,
                     diamKind = c('dib','dob'),
                     diamType = c('diameter','radius','xsection'),
                     standUp = FALSE,
                     ylab = NULL,
                     xlab = NULL,
                     rugReadings = c('SC','PC'), 
                     method = 'mono',
                     maxTrees = 50,
                     fileName = '',
                     pointsize.text = 12,
                     ... 
                    )
{
#---------------------------------------------------------------------------
#
#   This routine will create one or more lattice conditional graphs for
#   profiles of individual trees in the Dendrometry package. Note that you
#   can select several different ways to display these based on measurements
#   (e.g., dib vs dob, etc.), which are controlled by the arguments as listed
#   below.
#
#   Arguments...
#     ds = a list that must be generated from SplineVolume
#     TreeNos = the tree number to be plotted; if NA ot NULL, all trees will
#               be plotted
#     diamKind = either inside or outside bark
#     diamType = here we can choose to plot diameters, radius (just diameter/2)
#                or cross-sectional area
#     standUp = TRUE: the tree is displayed vertically; FALSE: it is displayed
#               horizontally in height
#     ylab = NULL or NA for default, otherwise a desired label
#     xlab = see ylab
#     rugReadings = a character vector of product codes where you would like
#                   rug marks--default are for base of live crown on our data;
#                   NULL or NA for no rug marks
#     method = one of the legal spline methods or NA/NULL for no spline plotted
#     maxTrees = a resettable maximum, this is just there in case there the
#                data have lots of trees and would create an absurd plot, you
#                can change it if you want to allow the display of more trees
#     fileName = '' for no hardcopy, otherwise a filename compatible with
#                hardcopyLattice
#     pointsize.text = see hardcopyLattice for details
#     ... = passed to xyplot
#
#   Returns...
#     A valid lattice object containing the plot
#    
#
#Author...									Date: 15-Nov-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   quick check...
#
    if(!is(ds, 'SplineVolume'))
      stop('You must supply a list from SplineVolume!')
    diamKind = match.arg(diamKind)
    diamType = match.arg(diamType)
    kfac = pi/(4*144)                  #all English units
    
#
#   "unlist"...
#
    segments = ds$segments
    trees = ds$trees

#
#   check that tree numbers are valid, and there are not too many...
#
    trees = trees[!trees$isForked,]               #discard forked trees
    if(is.null(TreeNos) || is.na(TreeNos))
      TreeNos = trees$TreeNo                      #not necessarily consecutive
    else {
      tns = intersect(TreeNos, trees$TreeNo)
      if(length(tns) == 0)
        stop('No tree numbers matching your request in argument \"TreeNos"')
      else
        TreeNos = tns
    }
    TreeNos = unique(TreeNos)

    numTrees = length(TreeNos)
    if(numTrees > maxTrees)
      stop(paste('You requested',numTrees,'trees for plotting, but maxTrees is',maxTrees))


#
#   determine the cumulative heights and diameters required for each tree and
#   create the new data frame for plotting...
#    
    TreeNo = NULL                                 #tree number
    height = NULL                                 #cumulative heights for all trees
    diam = NULL                                    #diameter
    prod = NULL
    for(i in TreeNos) {                           #one tree at a time
      segs = segments[segments$TreeNo == i,]
      TreeNo = c(TreeNo, segs$TreeNo)
      height = c(height, rev(segs$height))
      dx = rev(segs[,diamKind])
      if(diamType == 'radius')
        dx = dx/2
      else if(diamType == 'xsection')
        dx = dx*dx*kfac
      diam = c(diam, dx)
      prod = c(prod, rev(as.character(segs$prod)))      #from factor to character for plotting
    }
    TreeNo = ordered(TreeNo)
    prod = ifelse(is.na(prod), '  ', prod)
    df = data.frame(TreeNo, height, diam)          #all but prod or use stringsAsFactors=FALSE
    df$prod = prod                                 #keep it character for plotting below

#
#   the panel function for plotting...
#
    thePanel = function(x,y, ..., subscripts) {
      panel.xyplot(x,y, ...)
      if(standUp) {
        hgt = y
        dia = x
      }
      else {
        hgt = x
        dia = y
      }
      if(!is.null(method) && !is.na(method)) {            #plot the splines if desired
        sp.ht = splinefun(hgt, dia, method=method)
        h = seq(min(hgt), max(hgt), length=201)
        if(standUp)
          panel.lines(sp.ht(h), h, lty='dashed', col='grey50')
        else
          panel.lines(h, sp.ht(h), lty='dashed', col='grey50')
      }
      pcodes = df$prod[subscripts]
      pdx = ifelse(pcodes %in% rugReadings, hgt, NA) #this will catch NA & NULL rugReadings
      if(standUp)
        panel.rug(y=pdx)
      else
        panel.rug(pdx)
    }

#
#   and plot it...
#
    if(is.null(xlab) || is.na(xlab)) {
      if(standUp)
        xlab = diamType
      else
        xlab = 'height'
    }
    if(is.null(ylab) || is.na(ylab)) {
      if(standUp)
        ylab = 'height'
      else
        ylab = diamType
    }
    if(standUp)
      modelForm = formula(height ~ diam | TreeNo)
    else
      modelForm = formula(diam ~ height | TreeNo)
    plt = xyplot(modelForm, data = df,
                 panel = thePanel,
                 xlab = xlab, ylab = ylab,
                 ...
                )
    print(plt)
    
#
#   hardcopy?...
#    
    if(nchar(fileName) > 0) 
      hardcopyLattice(plt, fileName=fileName,
                      #columnsPerRow = columnsPerRow,
                          pointsize.text = pointsize.text#,
                      #    pdf.height=pdf.height,
                      #    pdf.width=pdf.width
                      )

    return(invisible(plt))
}   #taperPlot

      
