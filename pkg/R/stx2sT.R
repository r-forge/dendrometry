stx2sT = function(sv.list,
                  treeNo = 1,
                  centerOffset = c(x = 0.0, y = 0.0),
                  units.in = 'English',
                  units.out = 'metric',
                  prefix = 'odr77:tree',
                  runQuiet = TRUE,
                  ...
                 )
{
#---------------------------------------------------------------------------
#
#   This routine will convert a tree that has been processed through
#   SplineVolume and convert it into a sampSurf::standingTree object.
#   Here we use the full segment information from dendrometry for the taper
#   information and create the tree based on this.
#
#   Note: The dendrometry is normally only 6-8 segments. This is rather a
#         small number for subsequent spline interpolation. Also, the
#         dendrometry assumes a stump height, which is not assumed in the
#         taper data, and for some of the trees in the Old Reservoir, for
#         example, this can be 1.5ft, which is a substantial loss of 
#         information. We could "back-out" a stump diameter from a taper 
#         equation, but we'd need to fit the parameter(s). So, here we 
#         use the taper data for the object, but artificially set dbh
#         based on the observed value; though this is at odds with the 
#         taper relationship! I see no other way at the present and it
#         should not hurt anything in sampSurf.
#
#   See the companion program stx2sTs() for making a standingTrees
#   container object for multiple trees.
#
#   Arguments...
#     sv.list = a list directly from SplineVolume.
#     treeNo = the current tree number to extract from sv.list
#     centerOffset = see methods?standingTree; this must be specified here
#                    in order for all the other spatial points to be set
#                    set up correctly relative to the tree's center
#     units.in = the units that the dendrometry was take in, this should
#                always be "English"
#     units.out = the units for output, which must conform to the two
#                 potential choices in sampSurf
#     prefix = an identifier prefix for the tree id
#     runQuiet = TRUE: no feedback; FALSE: a little info
#     ... = gobbled
#
#   Returns...
#     a list invisibly with...
#       -- the segment data frame for the chosen tree
#       -- same for the tree-level data frame (one record)
#       -- the "standingTree" object
#       -- whether double bark thickness was used
#       -- volume i.b. or o.b.; redundant with above
#       -- average of the three methods of volume determination in SplineVolume
#
#Author...									Date: 16-May-2019
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
      
    if(length(treeNo) > 1)
      stop('Only one tree may be processed at a time!')

    ssUnits = unlist(sampSurf::.StemEnv$msrUnits) #adds the unwanted list names
    #must strip names or identical() tests in sampSurf .validObjects will fail...
    names(ssUnits) = NULL   #a very subtle possible problem averted
    units.in = match.arg(units.in, ssUnits)
    units.out = match.arg(units.out, rev(ssUnits))

#
#   determine the correct input and output unit conversions...
#
    in2cm = 2.54
    cm2in = 1/in2cm
    ft2m = 0.3048
    m2ft = 1/ft2m
    cuft2cum = 0.028316847
    cum2cuft = 1/cuft2cum
    if(units.in == 'English' && units.out == 'metric') {
      u.mult = sampSurf::.StemEnv$in2ft * ft2m                  #converts inches to meters
      u.vol = cuft2cum
      u.dbh = in2cm
      u.hgt = ft2m
    }
    else if(units.in == 'metric' && units.out == 'English') {
      u.mult = sampSurf::.StemEnv$cm2m * m2ft                   #converts cm to feet
      u.vol = cum2cuft
      u.dbh = cm2in
      u.hgt = m2ft
    }
    else if(units.in == 'English' && units.out == 'English') {
      u.mult = sampSurf::.StemEnv$in2ft
      u.vol = 1.0
      u.dbh = 1.0
      u.hgt = 1.0
    }
    else if(units.in == 'metric' && units.out == 'metric') {
      u.mult = sampSurf::.StemEnv$cm2m
      u.vol = 1.0
      u.dbh = 1.0
      u.hgt = 1.0
    }
    else
      stop('Something is wrong with "units.*"!')

#
#   extract the tree's aggregate info...
#
    trees.all = sv.list$trees
    tdx = which(trees.all$TreeNo %in% treeNo)
    if(length(tdx) == 0)
      stop(paste('No tree number',treeNo,'was found in the input list!'))
    tree = trees.all[tdx, ]

#
#   get segments for this tree...
#
    treeID = paste(prefix, treeNo, sep='.')
    segs.all = sv.list$segments
    seg.names = colnames(segs.all)
    sdx = which(segs.all$TreeNo %in% treeNo)
    if(length(sdx) == 0)
      stop(paste('Requested "treeNo"', treeNo, 'is not in the segment data frame!'))
    segs = segs.all[sdx, ]
    usedDBT = !identical(segs$dob, segs$dib)
    if(!runQuiet && usedDBT) 
      cat('\nBark deductions present: STX volumes are inside bark.')

    segs$dob = segs$dob * u.mult
    segs$dib = segs$dib * u.mult
    segs$height = segs$height * u.hgt
    
#
#   now create the standingTree object using the taper method...
#
    odx = order(segs$height, decreasing = FALSE)                #base at the top!
    if(usedDBT) {                         #only matters if volumes were i.b., could use dib always
      taper = segs[odx, c('dib', 'height')]
      dbh.u = with(tree, dbhob - dbThick) * u.mult              #note dbh here is in ft or m
    }
    else {
      taper = segs[odx, c('dob', 'height')]
      dbh.u = tree$dbhob * u.mult                               #note dbh here is in ft or m
    }
    names(taper) = c('diameter', 'height')
    stree = standingTree(taper,
                         centerOffset = centerOffset,
                         species = as.character(tree$Spp), 
                         treeID = treeID,
                         units = units.out
                        )
    #should not normally do the following, but the high stump throws this off...
    #note that some trees can have dbh > butt--measurement error?...
    stree@dbh = min(taper[1, 'diameter'], dbh.u)               #both must be in ft or m here
    stree@ba = .StemEnv$baFactor[units.out]*stree@dbh*stree@dbh    #must adjust this too!!
    #this check can be removed eventually if desired, though it does no hurt anything...
    if(stree@dbh >= taper[1, 'diameter'] || stree@dbh <= taper[nrow(taper), 'diameter'])
      cat('\n\n Bad dbh tree', treeID, '\n')
    
    if(usedDBT)
      volumes.are = 'inside bark'
    else
      volumes.are = 'outside bark'

#
#   just for kicks, average volume over all three methods actually compares to that from
#   standingTree()...
#
    aveVol = with(tree, mean(c(conicVol, natVol, monVol)))
    if(!runQuiet) {
      cat('\nTree volumes...')
      cat('\n  Spline volume from standingTree =', stree@treeVol)
      cat('\n  Average tree volume (3 methods below) =',  aveVol)
      cat('\n  SplineVolume conicVol =', tree$conicVol)
      cat('\n  SplineVolume natVol =', tree$natVol)
      cat('\n  SplineVolume monVol =', tree$monVol)
      cat("\nTree's pith location centerOffset: x =", centerOffset['x'], 'y =', centerOffset['y'])
      cat('\n')
    }
    
    return(invisible(list(segs = segs,
                          tree = tree,
                          stree = stree,
                          sv.name = deparse(substitute(sv.list)),
                          usedDBT = usedDBT,
                          volumes.are = volumes.are,
                          aveVol = aveVol
                         )
                    )  
          )
    
}   #stx2sT

