!-----------------------------------------------------------------
!    Routine: TreeMod.f90  
!
!    Module that contains tree oriented routines for SplineVolume
!    (with the exception of GetCtl).
!
!    **The ForkedAt array always stores (# of forks+1) measurements...
!        (1)        [reading number just below the fork]
!        (2)        [reading number for 1st msrment in first/smallest fork]
!        (3)        [reading number for 1st msrment in second/tallest fork]
!       etc...
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
module TreeMod

!
!All of the following variables are global to whatever routines use
!this module--no need to put them in commons.
!
!LRG's variables are first (ALL CAPS), paired with mine if changed;
!note that I've made things double precision...
!
logical*1                 LengthFlag       !flags a negative section length in tree
logical*1                 isForked         !flags a forked tree
character*2, pointer ::   GAMATH(:)        !product class/code
character*4               Spp              !species
character*60              title            !job title from STX control card #1 (ALFATH)
integer                   nTreeErrors      !# of trees processed with errors found
integer                   TreeNo           !tree number
integer                   nReadings        !# of measurements=#segments+1
integer                   nSegs            !# of segments=nReadings-1
integer                   nInterrupts      !# of fork or new setup of dendrometer
integer                   nForkedPoints    !dimension of ForkedAt()
integer                   METH             !tree measurement instrument/method
integer, pointer ::       ForkedAt(:)      !**reading numbers for forks if isForked
                                           !see comments above
real*8                    DBH,BKA,BKB      !dbh & bark thicknesses
real*8                    dbThickness      !double bark thickness
real*8                    dbhob            !dbh outside bark==DBH
real*8                    dbhib            !dbh inside bark
real*8                    totLength        !total length
real*8                    totNatCSVol      !total tree natural cubic spline vol c.f.
real*8                    totAkimaVol      !total tree Akima volume in c.f.
real*8                    totHermiteVol    !total tree Hermite volume in c.f.
real*8                    totConicVol      !total tree conic section volume in c.f.
real*8, pointer ::        TGRADS(:),FGRADS(:),SINELV(:)
real*8, pointer ::        SegLength(:), dib(:), dob(:)
real*8, pointer ::        CumLength(:)     !cumulative length==height if not forked
real*8, pointer ::        NatCSVol(:)      !natural cubic spline volume in c.f.
real*8, pointer ::        AkimaVol(:)      !Akima cubic spline volume in c.f.
real*8, pointer ::        HermiteVol(:)    !Hermite cubic spline volume in c.f.
real*8, pointer ::        ConicVol(:)      !conic section volume in c.f.
real*8, pointer ::        instRange(:)     !instrument range to tree
real*8, pointer ::        HtAboveHorizontal(:)    !says it all!
real*8, pointer ::        xsArea(:)        !stem cross-sectional area at msrmnts



contains !***************************************************************


include 'GetCtl.f90'                 !reads STXMOD job control records
include 'GetTree.f90'                !reads tree & associated dendro records


end Module TreeMod
