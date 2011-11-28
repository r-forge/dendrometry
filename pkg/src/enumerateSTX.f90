subroutine enumerateSTX(nTrees, nTotSegs, nTotReadings, inFile, nfin)
!======================================================================*
!--------------- P r o g r a m   D e s c r i p t i o n ----------------*
!======================================================================*
!
!     {enumerateSTX.f90}
!
!======================================================================*
!--------------- V a r i a b l e   D e c l a r a t i o n  -------------*
!======================================================================*
!

use ioMod
use TreeMod


IMPLICIT NONE

integer, intent(in) :: nfin

integer, intent(inout) ::   nTrees        !total number of trees processed
integer, intent(inout) ::   nTotSegs      !total number of segments in nTrees
integer, intent(inout) ::   nTotReadings  !total dendrometer readings in nTrees

character(nfin), intent(in) ::    inFile

!
!======================================================================*
!-------------------  B e g i n   t h e   R o u t i n e  --------------*
!======================================================================*
!


Finished = .false.
nCtlRecs = 10                   !# of STX control cards/records

!
!open files...
!
open(unit=lunIn, file=inFile, status='old')              !.stx

!
!read in control cards...
!
call GetCtl()

!
!go tree-by-tree through the input (STX) file...
!
nTrees = 0
nTotSegs = 0
nTotReadings = 0
do
    call GetTree()
    if(Finished) exit                                    !no more trees left

!
!   clean up and get ready for the next tree...
!
    nTrees = nTrees + 1
    nTotSegs = nTotSegs + nSegs
    nTotReadings = nTotReadings + nReadings
    
    if(associated(TGRADS)) deallocate(TGRADS)
    if(associated(FGRADS)) deallocate(FGRADS)
    if(associated(SINELV)) deallocate(SINELV)
    if(associated(GAMATH)) deallocate(GAMATH)
    if(associated(dob)) deallocate(dob)
    if(associated(dib)) deallocate(dib)
    if(associated(SegLength)) deallocate(SegLength)
    if(associated(instRange)) deallocate(instRange)
    if(associated(HtAboveHorizontal)) deallocate(HtAboveHorizontal)
    if(associated(xsArea)) deallocate(xsArea)
    if(associated(NatCSVol)) deallocate(NatCSVol)
    if(associated(AkimaVol)) deallocate(AkimaVol)
    if(associated(HermiteVol)) deallocate(HermiteVol)
    if(associated(ConicVol)) deallocate(ConicVol)
    if(associated(ForkedAt)) deallocate(ForkedAt)
enddo

close(lunIn)
return
end subroutine enumerateSTX

