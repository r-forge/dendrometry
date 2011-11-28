subroutine testDendro(fileName, nfn, stxFile, nstx, gg2  &
                     )
!
!Author...									Date: 2-Nov-2011
!	Jeffrey H. Gove
!	USDA Forest Service
!	Northern Research Station
!	271 Mast Road
!	Durham, NH 03824
!	jhgove@unh.edu
!	phone: 603-868-7667	fax: 603-868-7604
!-----------------------------------------------------------------------------
!
use GrosenbaughMod

implicit none

real*8, intent(inout) :: gg2

integer, intent(in) :: nstx, nfn

character(nstx), intent(inout) :: stxFile
character(nfn), intent(in) :: fileName 
!
!======================================================================*
!------------- C o m m o n   A r e a s --------------------------------*
!======================================================================*
!
Real*8   gg,cg,w,z,q,b,u,g
common/lrg/gg,cg,w,z,q,b,u,g


!
!set Barr & Stroud constants...
!
call FPSetup()
print "('gg =',f12.4)", gg

gg2=gg



stxFile = fileName//'.stx'

return
end subroutine testDendro
