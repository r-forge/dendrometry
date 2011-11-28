SplineVolume = function(fileName='odr_wp_1977.stx', filePath = getwd(), dataPath = 'data',
                        reportPath = 'reports', splineMethods = c('natural','monoH'),
                        ...
                                    )
{
#---------------------------------------------------------------------------
#
#   This routine is the interface to the fortran program SplineVolume, which
#   does all of the work reading and translating LRG's STX file format. Please
#   see the comments in the fortran code for details. Much of what I did comes
#   directly from code borrowed from STX, that I augmented.
#
#   Note below that I pass a few arrays to be filled directly by the fortan
#   code. This approach was not used though in the end since the R-fortran
#   interface does not allow passing character vectors, which are required
#   for the segment product codes. So the simplest thing to do was simply
#   to read the output files from the fortran back in here to create the
#   basic data frame. We could set up an index array matching the codes in
#   the fortran (like a factor) and decode them back here in the R code,
#   but there's really no sense to that.
#
#   Note also that you can choose to estimate spline volumes for any one of the
#   methods in R's spline and splinefun routines in the utils package. These
#   will be add in a single column to the data frame, but not to the files
#   residing on disk. See the code below for an example of how you might add
#   other spline volumes from this package.
#
#   Arguments...
#     fileName = the STX "cards" file name
#     filePath = the full file path to the current (R session) working directory
#     dataPath = the subfolder, if any, for the STX input file
#     reportPath = the subfolder for the output from SplineVolume
#     splineMethods = one or more of the legal methods in the spline and
#                     splinefun functions
#     ... = gobbled
#
#   Returns...
#     A named list (invisibly--this can be quite large!) with...
#        z = a list of the arguments passed to SplineVolume (fortran)
#        trees = the tree data frame
#        segments = the segments data frame
#
#   Please check the results carefully. STX format can allow for some odd
#   tree measurments, and I can not guarantee that I have accounted for
#   every possibility. The *.rpt file in the reportPath subdirectory should
#   allow for easy checking of the results.
#
#Author...									Date: 2-Nov-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   check file names...
#
    fn = unlist(strsplit(fileName, '\\.'))[1]
    fn.stx = file.path(filePath,dataPath,fileName)
    if(!file.exists(fn.stx))
      stop(paste('There is no file named:',fn.stx))
    if(nchar(fn.stx) > 255)
      stop(paste('data file name is longer than max of 255 bytes =',nchar(fn.stx)))
    n.stx = nchar(fn.stx)
    fn.dvf = file.path(filePath,reportPath,paste(fn,'tre',sep='.'))  #tree file
    fn.dvl = file.path(filePath,reportPath,paste(fn,'seg',sep='.'))  #segments file
    fn.out = file.path(filePath,reportPath,paste(fn,'rpt',sep='.'))  #report
    n.fn = nchar(fn.out)

#
#   this has been added so we can find out how many trees and total segments
#   are in the data file, allowing us to pre-dimension the arrays that will
#   be passed back from SplineVolume...
#
    nTrees = 0
    nSegs = 0
    nReadings = 0
    recs = .Fortran('enumerateSTX',
                    nTrees = as.integer(nTrees),
                    nSegs = as.integer(nSegs),
                    nReadings = as.integer(nReadings),
                    stxFile = as.character(fn.stx),
                    n.stx = as.integer(n.stx)
                   )

#
#   okay, now set up the arrays that will be passed to and filled by SplineVolume
#
    nSegs = recs$nSegs
    nTrees = recs$nTrees
    nReadings = recs$nReadings
    tgrads = rep(0, nReadings)
    fgrads = rep(0, nReadings)
    sinelv = rep(0, nReadings)

   
    
#
#   get the tree and segment information from the STX file...
#
    z = .Fortran('SplineVolume',
                 nTrees = as.integer(nTrees),
                 nSegs = as.integer(nSegs),
                 nReadings = as.integer(nReadings),
                 stxFile = as.character(fn.stx),
                 outFile = as.character(fn.out),
                 dvfFile = as.character(fn.dvf),
                 dvlFile = as.character(fn.dvl),
                 n.stx = as.integer(n.stx),
                 n.fn = as.integer(n.fn),
                 tgrads = as.double(tgrads),
                 fgrads = as.double(fgrads),
                 sinelv = as.double(sinelv)
                )

#
#   now read in the output files just created...
#
    trees = read.table(fn.dvf, TRUE)
    segments = read.table(fn.dvl, TRUE)


#
#   now calculate spline volumes from the spline package...
#    
    nTrees = z$nTrees
    TreeNo = trees$TreeNo
    for(sm in splineMethods) {
      treeSegVols = NULL                            #segment volumes for all trees
      totVol = rep(NA,nTrees)                       #total volume for each tree
      for(i in seq_len(nTrees)) {                   #one tree at a time
        tn = TreeNo[i]                              #number can be anything
        segs = segments[segments$TreeNo == tn,]
        if(trees[i,'isForked']) {                   #don't try to spline forked trees
          nsegs = nrow(segs) - 1
          tsegs = rep(0, nsegs)
        }
        else 
          tsegs = splineSegVol(segs, method=sm)
        treeSegVols = c( treeSegVols, c(rev(tsegs), 0.0) )  #add the zero placeholder for first measurement
        totVol[i] = sum(tsegs)
      }
      vname = paste(substr(sm,1,3),'Vol',sep='')    #new column name for both data frames
      segments[,vname] = round(treeSegVols, 3)      #no need for more digits
      trees[,vname] = round(totVol,3)               #same
    } #spline methods

    rl = list(z=z, trees=trees, segments=segments)
    class(rl) = c('SplineVolume',class(rl))
    return( invisible(rl) )
}   #SplineVolume

