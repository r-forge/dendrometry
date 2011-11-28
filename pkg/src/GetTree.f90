subroutine GetTree()
!-----------------------------------------------------------------
!    Routine: GetTree.f90 
!
!    This routine reads in one tree's worth of information from the
!    combination of tree plus accompanying dendrometer cards in STX
!    file format.
!
!    Notice that allocatable (pointer) arrays are used so that they
!    are custom fits for each tree depending upon the number of 
!    segments in the tree. ***But, this means that you must remember
!    to deallocate this space before calling this routine from the
!    second tree through the last.*** 
!
!   Note that there is nothing in this routine that is tied to the
!   Barr and Stroud other than the names of the variables. Therefore, 
!   it can be used to input information taken by other measuring 
!   devices, or by direct measurement. Subsequent routines can
!   do the device-specific processing. This is similar to STXMOD.
!   (JHG June 1998)
!
!Author...                                    Date: 26-Nov-97
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
!locals...
!LRG's variables first in (ALL CAPS), paired with mine if changed...
character*1    CERT,JIM,TERM
character*4 BETATH
integer        KREENO,JQ,KPI,LST,MBK,MUL,JAM,PLOTNO 
real        UMAXL,UDORT,XTRA,XTRB,DEDT,GRTH

character*2 gam(4)                !used to input eacg tree record only
integer nRecs                    !# of tree records
integer    idxFork                    !fork index keeper-tracker-of-er
integer i,j,ns                    !misc. counters
integer    nleft,ntop,ngood,ngroup    !for segment arithmetic
real tg(4),fg(4),selv(4)        !used to input each tree record only
!
!-----------------------------------------------------------------
!
!take care of some initializing business first...
!
nReadings = 0            !global, so zero it first thing!
isForked = .false.      !assume the best to begin--the eternal optimist

!
!    first read in the tree card...
!
read(lunIn, 9)KREENO,JQ,KPI,LST,CERT,BETATH,DBH,JIM,METH,MBK,MUL,    &
              JAM,BKA,BKB,UMAXL,UDORT,XTRA,XTRB,TERM,DEDT,PLOTNO
9 FORMAT(I4,I1,I4,I1,A1,A4,1X,F5.1,A1,3I1,I1,2F4.1,1X,F3.0,F3.3,    &
         2F15.0, A1, F4.0, I4)

!
!    and see if we have any trees left to process...
!  
IF (KREENO .GE. 9999) then
    Finished = .true.
    return
endif                                   
TreeNo = KREENO         
Spp = BETATH 
         
                                                   
!
!    now read all of the dendrometer cards for this tree to figure out
!   how big the arrays must be to hold the measurements...                                                                       
!
nRecs = 0
nInterrupts = 0
do 
    READ(lunIn,10)KREENO,JQ,(tg(i),fg(i),selv(i), gam(i), i=1,4),    &
                  TERM, GRTH, PLOTNO                         
10    FORMAT(I4, I1, 6X, 4(2F4.1, F5.4, A2), A1, F4.0, I4 )
    nRecs = nRecs + 1
    do i=1,4
        if((tg(i).ne.0.0 .or. fg(i).ne.0.0) .and. selv(i).ne.0.0)    & 
            nReadings = nReadings + 1
    enddo
    if(TERM.eq.'+') nInterrupts = nInterrupts + 1
    if(TERM.eq.'*') exit
enddo

!
!we want to include the measurement point below the fork in ForkedAt too,
!so add one for this...
!
if(nInterrupts.gt.0) then
    nForkedPoints = nInterrupts + 1
else
    nForkedPoints = 0
endif

!
!    backup so we can read them again; then allocate tree arrays...
!                
do i=1,nRecs
    backspace lunIn
enddo
allocate(tgrads(nReadings),fgrads(nReadings),sinelv(nReadings),    &
         gamath(nReadings))

if(nForkedPoints.gt.1) then
    allocate(ForkedAt(nForkedPoints))       !only if needed
    ForkedAt = 0                            !all null to begin with
    isForked = .true.                       !yes, a forked tree!
endif
!
!    now read all of the dendrometer cards...                                                                       
!
idxFork = 1                                 !pointer to ForkedAt() cells
ns = 1                                      !reading # to start at
do j=1,nRecs
    READ(lunIn,10)KREENO,JQ,                            &
                 (tg(i),fg(i),selv(i),gam(i), i=1,4),    &
                  TERM, GRTH, PLOTNO
!
!    number of measurements on this record...
!
    ngroup = 0
    do i=1,4
        if((tg(i).ne.0.0 .or. fg(i).ne.0.0) .and. selv(i).ne.0.0)    & 
            ngroup = ngroup + 1
    enddo
!
!    figure out how to put them in the arrays...
!
    nleft = nReadings - (ns    - 1)            !readings available to assign
    if(nleft.ge.ngroup) then                !in groups
        ntop = ns + (ngroup - 1)
        ngood = ngroup
    else
        ntop = nReadings
        ngood = nleft
    endif                            
    TGRADS(ns:ntop) = tg(1:ngood)
    FGRADS(ns:ntop) = fg(1:ngood)
    SINELV(ns:ntop) = selv(1:ngood)
    GAMATH(ns:ntop) = gam(1:ngood)
!
!    if at a fork, save the reading number for later...
!
    if(TERM.eq.'+') then
        if(idxFork.eq.1) then           !save reading # just below the fork
            ForkedAt(idxFork) = ntop
            idxFork = idxFork + 1
        endif
        ForkedAt(idxFork) = ntop + 1    
        idxFork = idxFork + 1
    endif
    ns = ntop + 1                       !reset base reading pointer
enddo                
nSegs = nReadings - 1                   !one less segments than readings

return
end subroutine GetTree
