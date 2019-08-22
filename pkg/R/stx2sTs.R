stx2sTs = function(sv.list,
                   buffTr,
                   units.in = 'English',
                   units.out = 'metric',
                   prefix = 'odr77:tree',
                   inhibitDist = 3,
                   runQuiet = TRUE,
                   ...
                  )
{
#---------------------------------------------------------------------------
#
#   This routine will convert the full output from a run of SplineVolume()
#   to a list of standingTree() objects and finally to a standingTrees()
#   object. Note that this is a prototype, and generates tree coordinates
#   within a bufferedTract buffer.
#
#   The spatial coordinates are drawn within the tract using an inhibition
#   process; the "spatial" package code is used here because it is faster than
#   the corresponding routine in the spatstats package.
#
#   Arguments...
#     sv.list = a list from SplineVolume()
#     buffTr = a bufferedTract object
#     units.in = the units that the dendrometry was take in, this should
#                always be "English"
#     units.out = the units for output, which must conform to the two
#                 potential choices in sampSurf
#     prefix = an identifier prefix for the tree id
#     inhibitDist = inhibition distance in units.out units
#     runQuiet = TRUE: no feedback; FALSE: a little info
#     ... = gobbled
#
#   Returns...
#     a list invisibly with...
#       -- a list of standingTree objects
#       -- the same trees in a standingTrees container object
#       -- the tree ID prefix used for the group of trees
#
#   The reason for passing back a list too is that it facilitates joining one
#   or more populations of standingTree object by simply concatenating
#   the associated lists. One can not do this with standingTrees objects
#   without running into problems with the bounding box, etc., which will
#   have to be recreated.
#
#Author...									Date: 20-May-2019
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   a few checks...
#
    if(!is(sv.list, 'SplineVolume'))
      stop('sv.list must be of class "SplineVolume"!')
      
    if(!require(sampSurf, quietly=TRUE))
      stop('You must have package sampSurf installed!')
      
    if(units.out != buffTr@units)
      stop('***>Units for output (units.out) do not match those for bufferedTract!')
      
    if(!require(spatial, quietly = TRUE))
      stop('You must have package spatial installed!')
      
      
    trees = sv.list$trees
    n.trees = nrow(trees)
    treeNums = trees$TreeNo
    treeIDs = paste(prefix, treeNums, sep='.')
    
#
#   generate spatial locations within the tract bbox...
#
    bbox = buffTr@bufferRect
    ppregion(bbox[1,1], bbox[1,2], bbox[2,1], bbox[2,2]) #set the domain for SSI()
    xy = SSI(n.trees, inhibitDist)
    xx = xy$x
    yy = xy$y

    
#
#   loop through all trees in the list...
#
    strees.list = vector('list', n.trees)
    names(strees.list) = treeIDs
    if(!runQuiet)
      cat('\nProcessing tree number...')
    for(i in seq_len(n.trees)) {
      if(!runQuiet)
        cat(treeNums[i],', ', sep='')
      treeNo = treeNums[i]                         #not always consecutive, could be missing numbers
      centerOffset = c(x = xx[i], y = yy[i])
      zz = stx2sT(sv.list, treeNo, centerOffset, units.in, units.out, prefix, TRUE) #quietly always
      strees.list[[i]] = zz$stree
    }
    
    if(!runQuiet) {
#      cat('\nGenerate a standingTrees object from the return list...')
      cat('\n')
    }
    
    strees = standingTrees(strees.list)
    
    return(invisible(list(strees.list = strees.list,
                          strees = strees,
                          prefix = prefix
                         )
                    )
          )

}   #stx2sTs


