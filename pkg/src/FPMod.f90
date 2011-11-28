!-----------------------------------------------------------------
!	Routine: FPMod.f90  
!
!	Module that contains tree oriented routines for the Barr and 
!	Stroud Model FP-12 and FP-15 dendrometers.
!
!	In Grosenbaugh-ese: D=diameter; R=range; E=elevation/height
!						above horizontal.
!
!	***Always call FPSetup() onec/first to set constants.
!
!Author...									Date: 8-Dec-97
!	Jeffrey H. Gove
!	USDA Forest Service
!	Northeastern Forest Experiment Station
!	P.O. Box 640
!	Durham, NH 03824
!	jhgove@christa.unh.edu
!	phone: 603-868-7667	fax: 603-868-7604
!-----------------------------------------------------------------
!
module FPMod



contains !********************************************************

include 'FP_ProcessTree.f90'		!process measurements for B&Stroud


end Module FPMod
