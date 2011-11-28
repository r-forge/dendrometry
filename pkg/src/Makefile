#----------------------------------------------------------------
#
#***** Please Note...
#      Building directly through R is now the preferred method
#      because things can change over the version of R in RFLAGs,
#      and also Rflags are different between Ubuntu and Fedora.
#      Finally, architecture could change (32-bit vs 64-bit) so
#      it is better to just use R's built-in flags for any given
#      Linux version and architecture.
#***** JHG, 5-Mar-2010
#
#
#
#   This file is meant to be used by R to compile the Dendrometry
#   fortran code. As noted above, it relies on R to supply all of the 
#   correct flags for building. Thus, please do not invokee this from
#   the command line to try to build the so; you can, however, invoke
#   it for cleaning from the appropriate ../src/ directory...
#
#    "make -fMakefile Dendrometry.so" PLEASE DON'T DO THIS!!!!
#    "make -fMakefile clean" to clean before building often helps.
#
#Author...									Date: 3-Nov-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	270 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#  
#----------------------------------------------------------------
#
#  compile the modules first, separately, then put it together...
#

Dendrometry.so :
	R CMD COMPILE ioMod.f90
	R CMD COMPILE TreeMod.f90
	R CMD COMPILE GrosenbaughMod.f90
	R CMD COMPILE FPMod.f90
	R CMD SHLIB -o Dendrometry.so SplineVolume.f90 ioMod.f90 TreeMod.f90 GrosenbaughMod.f90 FPMod.f90 \
enumerateSTX.f90


#
# clean safely...
#
.PHONY : clean
clean :
	rm -f *.mod *.o *.so
