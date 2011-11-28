subroutine GetCtl()
!-----------------------------------------------------------------
!    Routine: GetCtl.f90  
!
!    This routine simple dummy reads any control cards that may be
!    at the beginning of an STX input file.
!
!Author...                                    Date: 4-Dec-97
!    Jeffrey H. Gove
!    USDA Forest Service
!    Northeastern Forest Experiment Station
!    P.O. Box 640
!    Durham, NH 03824
!    jhgove@christa.unh.edu
!    phone: 603-868-7667    fax: 603-868-7604
!-----------------------------------------------------------------
!
use ioMod

implicit none
integer i

!
!    get all of the control records out of the way...
!
if(nCtlRecs.gt.0) then
    if(nCtlRecs.gt.10) then
        write(kbout,*)'Too many control records specified (10 max)!'
        stop    'please change your input'
    else
        do i=1,nCtlRecs
            if(i.eq.1) then
                read(lunIn,'(t5,a60)')title
            else
                read(lunIn,*)
            endif
        enddo
    endif
endif

end subroutine GetCtl
