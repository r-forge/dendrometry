splineSegVol = function(segments, method='mono',  ... )
{
#---------------------------------------------------------------------------
#
#   This function will simply fit a spline function using any of the
#   available methods in the spline & splinefun routines. It will then
#   simply calculate the volumes per segment by the method requested.
#
#   Arguments...
#     segments = a set of tree segments from SplineVolume for ONE tree only
#     method = one of the legal spline methods in the stats package.
#     ... = to be passed on to the spline routines as desired
#
#   Returns...
#     the segments data frame with the new spline volumes calculated
#
#Author...									Date: 8-Nov-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   just extract the diameters and heights...
#   
    height = rev(segments$height)
    diameter = rev(segments$dib)

#
#   spline function to be integrated (cross sectional area for volume)...
#    
    taper.spline = splinefun(height, diameter, method=method, ...)
    v.spline = function(hgt) { 
      diam = taper.spline(hgt)
      csa = pi*diam*diam/(4*144)
      return(csa)
    } #v.spline

#
#   calculate volume for each segment, pad first with zero (should probably be NA)...
#
    nSegs = nrow(segments) - 1
    segVol = rep(NA, nSegs)
    for(i in seq_len(nSegs)) {
      segVol[i] = integrate(v.spline, height[i], height[i+1])$value       #integrates for each segment
    }
    #vname = paste(method,'Vol',sep='')
    #segments[,vname] = c(rev(segVol), 0.0)                #zero placeholder for first measurement

    return(segVol)
} #splineSegVol

