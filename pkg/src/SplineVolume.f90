subroutine splineVolume(nTrees, nTotSegs, nTotReadings, &
                        inFile, outFile, dvfFile, detailFile, nfin, nfrpt, &
                        Rtgrads, Rfgrads, Rsinelv)
!======================================================================*
!--------------- P r o g r a m   D e s c r i p t i o n ----------------*
!======================================================================*
!
!     {SplineVolume.f90}
!
!     ***Please read this update 8-Nov-2011...
!
!        Now I have added some of what is mentioned below with regard to
!        reading in the number of segments in another routine (enumerateSTX.f90)
!        and pre-dimensioning the arrays. Note that tgrads, fgrads, and sinelv
!        are now all passed back this way. 
!  
!        Now the problem is that the character vector for product for each
!        segment would need to be flattened into one long character variable
!        which can only have length 255. This is a huge limitation, as soon
!        as we get more than 255/2 segments (products are 2 bytes) we will
!        overflow the variable. R can't handle fortran character vectors 
!        right now. We could set up an index array with each of the possible
!        product codes like a factor in R, and then decode when we send it back
!        to R, but this seems to be needlessly complicated when everything is 
!        written to the files and can simply be read into R.
!
!        So, I have kept the above changes, but have added all of the B&S
!        dendrometry readings and computed distances and heights to the
!        detailed file, which can now simply be read into the R calling
!        routine.
!
!        I think this is adequate. I also took some of the tree-level
!        variables out of the detail segment file to make more room
!        for the above w/o the clutter of extra repetative columns. You
!        can merge them all in R easily enough. 
!
!        As a result, SplineVolume differs considerably from the old original
!        code w/r to what is saved to files.
!
!     ***Please read this update 4-Nov-2011...
!
!        This version has been adapted and modified to be called from
!        within R. Akima and Hermite splines have been dropped because
!        I used code from the imsl library that we had a license for
!        back when this was written. Probably I should drop the cubic
!        splines too as I do not think we are supposed to distribute
!        the NR sorce code, though it may be okay. Akima and Hermite
!        are both available in R, and there will be an R routine to
!        calculate them after this has been run. One can also compare 
!        the NR cubic spline results with those in R if desired. The conic
!        sections are probably as good as anything if they are kept 
!        relatively short.
!
!        The fortran code uses allocatable arrays extensively because
!        we don't know how many trees or readings within trees there
!        will be ahead of time; and in raw stx format, this is a pain
!        to figure out. So each tree was read, processed and written
!        within one allocation/deallocation cycle and does no lend
!        itself, in the end, to returning this information as arguments
!        to this function. But everything is written to files that have been
!        adapted to be easily read with R's read.table. So the R interface
!        routine does just that after calling this fortran code. Of course, 
!        the drawback to this approach is the loss of precision, but it should
!        be okay for these types of measurements.
!
!        One last thought here is that we could write a separate little routine
!        that called getTree to determine the total number of segments (and trees)
!        that could be callable from R so that we could know the sizes of the
!        arrays required beforehand if need be. But let's leave that for another
!        time if someone really needs full precision! (Done, please see addendum 
!        above and enumerateSTX.f90.)
!
!        The stx/dendrometer stuff is not simple to read in, and I can not 
!        guarantee the correctness of all possible options for trees, so
!        please scrutinize the output closely. Notably, forked trees are
!        always problematic. So use this cautiously please. JHG.
!
!     **************
!
!     This is a simple front-end which calls some routines "lifted"
!     from STXMOD used to calculate tree diameter, range, and height
!     above horizontal for any section using a Barr and Stroud Fp-12
!     or Fp-15. It also calculates the following for each section...
!       --length in feet
!       --splined cf volume using natural cubic splines (Numerical Recipes)
!       --splined cf volume using Akima splines (IMSL)
!       --splined cf volume using Hermite splines (HValentine--old IMSL)
!       --cf volume based on conic sections (frustrum of a cone and cone)
!     totals are given for the tree for these latter quantities. 
!
!     Three output files are created...
!       1. The report file with tables that is best viewed in Word
!       2. A data file with the tree indentifiers and totals for S-Plus
!       3. A data file with detailed section info and Akima coefficients
!
!     This program is very limited in its scope compared to STX; that is, 
!     very few of the control options found in STX are useful here. For
!     example, only one type of bark option is used and it is always used
!     and that is option #2. If you have bark thicknesses, the dibs reported
!     will be calculated via this formula. If no thicknesses are taken, then
!     the DBT==0 and you get outside bark diameter (***>this is in contrast
!     with STXMOD which gives you 0.9*dob as a default if no DBT is taken!).
!
!     About the only error that is checked for now on individual trees is
!     for sections with negative lengths. These can occur for any number
!     of reasons, most obviously because of decreasing SINELV readings; but
!     oddities in TGRADS/FGRADS can also throw this off. Other checks might
!     also be added in the future.
!
!     Another place where it would be useful to modify this routine is 
!     for the other instruments. It would be simple to allow direct
!     measurements of course, and the tele-relascope would also be of
!     interest since we own one! Right now, METH<=1 is all that is allowed;
!     i.e., the B&S FP-12 or FP-15.
!
!     Also, it would be nice to use the log length info on the last control
!     card and partition up the tree into merchantable units as STXMOD does,
!     but based on spline interpolation here.
!
!   **On trees where the dendrometer has been moved (only once!) there will
!     be only ONE "+" TERM card. The program currently looks at this as a
!     forked tree, but needs to be able to distinguish this special case/use
!     of the "+"; why some other character was not originally used for 
!     new setup is beyond me!**
!
!----Author...
!     Jeffrey H. Gove             (603) 868-7667
!     USDA Forest Service
!     Forestry Sciences Lab       26-Dec-97
!     P.O. Box 640
!     Durham, N.H. 03824          jhgove@christa.unh.edu
!
!======================================================================*
!-------------- V a r i a b l e    D i c t i o n a r y ----------------*
!======================================================================*
!
!-->Variables...
!
!     Name                   Description
!   --------          ---------------------------------------------
!     tgrads            tgrad reading from Barr and Stroud
!     fgrads            fgrad reading from Barr and Stroud
!     sinelv            sine of elevation reading from B&S
!     meth              from tree card column #23, instrument type
!     dbhob             diameter at breast height outside bark
!     dbhib             diameter at breast height inside bark
!     dob               diameter calculated from dendrometer measurements
!      dib               inside bark diameter from dob and Grosenbaugh's
!                       bark model #2
!     Range             distance from instrument to tree
!     HtFromHoriz       height of dob from horizontal
!
!======================================================================*
!--------------- V a r i a b l e   D e c l a r a t i o n  -------------*
!======================================================================*
!
!program SplineVolume

use ioMod
use TreeMod
use GrosenbaughMod
use FPMod


IMPLICIT NONE

integer, intent(in) :: nfin, nfrpt

character(nfin), intent(in) ::    inFile
character(nfrpt), intent(in) ::    outFile,dvfFile,detailFile


integer, intent(inout) ::   nTrees        !total number of trees processed
integer, intent(inout) ::   nTotSegs      !total number of segments in nTrees
integer, intent(inout) ::   nTotReadings  !total dendrometer readings in nTrees

real*8, intent(inout) :: Rtgrads(nTotReadings), Rfgrads(nTotReadings), Rsinelv(nTotReadings)


integer    nLinesOut           !# of lines written for the current page
integer    nHeaderLines        !# of lines/tree for header & summary
integer    nTitleLines         !# of lines for page title

integer    begSeg, endSeg      !offsets for storing segment info 


!
!======================================================================*
!------------- C o m m o n   A r e a s --------------------------------*!
!======================================================================*
!
Real*8   gg,cg,w,z,q,b,u,g
common/lrg/gg,cg,w,z,q,b,u,g

real*8    kFactor, pi
common/jhgconstants/kFactor

!
!======================================================================*
!----------- V a r i a b l e   I n i t i a l i z a t i o n  -----------*
!======================================================================*
!
! see modules...

!
!======================================================================*
!-------------------  B e g i n   t h e   R o u t i n e  --------------*
!======================================================================*
!
pi = 2.0d0*dacos(0.0d0)
kFactor = pi/(4.0d0*144.0d0)        !basal area conversion constant--English



Finished = .false.
nCtlRecs = 10                   !# of STX control cards/records
nLinesPP = 50                   !set page length
nHeaderLines = 8                !set table header/summary length
nTitleLines = 2                 !set page title line size

write(kbout,5)
5   format(t20,"Processing of Barr and Stroud Dendrometry data"/    &
           t30,"in STX-formatted data files"/)

!
!open files...
!
open(unit=lunIn, file=inFile, status='old')              !.stx
open(unit=lunOut, file=outFile, status='replace')        !.out 
open(unit=lunDVF, file=dvfFile, status='replace')        !.dvf
open(unit=lunDetail, file=detailFile, status='replace')  !.dvl
write(kbout,16)trim(adjustl(outFile)),trim(adjustl(dvfFile)),trim(adjustl(detailFile))
16 format('Output report will be written to file...'/a   &
          /'Dendrometry summary volumes will be written to file...'/a   &
          /'Detail section information will be written to file...'/a/)



!
!set Barr & Stroud constants...
!
call FPSetup()


!
!read in control cards...
!
call GetCtl()
write(lunOut,"(t20,a60/)")title
nLinesOut = nTitleLines                           !for title
write(kbout,"('--->',a60)")title


write(lunDVF, 20)
20 format('TreeNo ','Spp ','dbhob ', 'dbThick ','isForked ','nReadings ','totLen ', &
               !'ncsVol ',
               'LenFlag ','conicVol ')

write(lunDetail, 22)
22 format('TreeNo ','Reading ','Spp ', &  !'dbhob ', 'dbThick ','isForked ',
          'dob ','dib ','prod ','tgrad ','fgrad ','sinelv ','instRange ','hgtHoriz ', &          
          'segLen ', 'height ','conicVol ') !'ncsVol ',)

!
!go tree-by-tree through the input (STX) file...
!
nTrees = 0
nTotSegs = 0
nTreeErrors = 0
nTotReadings = 0
do
    call GetTree()
    if(Finished) exit                                    !no more trees left

!
!   check to see if a new page is required...
!
    if((nLinesOut+nReadings+nHeaderLines).gt.nLinesPP) then
        write(lunOut,"(a1)")achar(12)
        write(lunOut,"(t20,a60/)")title
        nLinesOut = nTitleLines                         !for title
    endif
    nLinesOut = nLinesOut + nReadings + nHeaderLines    !add another table

!
!    calculate double bark thickness and dbh inside bark (from TreeMod)...
!
    dbThickness = BKA + BKB
    dbhob = DBH
    dbhib = dbhob - dbThickness

!
!   print the report while actually processing the tree...
!
    write(lunOut,50)TreeNo,Spp,dbhob,dbThickness,isForked
50  format('Tree number: ',i5,5x,'Species: ',a4,5x,                 &
           'Dbhob = ',f5.1,5x,'D.B.T. = ',f5.1,5x,                  &
           'Forked = ',L1/90('=')) !111('='))
    if(meth.le.1) then
        call FP_ProcessTree()
    else
        write(lunOut,75)
75      format('****>Invalid measurement instrument--this tree'     &
               ' will not be processed'/102('-')/)
        cycle
    endif
    write(lunOut,100)nReadings,totLength,                          &
                     !totNatCSVol,totAkimaVol,totHermiteVol,   &        
                     totConicVol
100 format(90('-')/'Number of Readings = ',i5,t64,'Totals:',t72,   &
           f5.1,t80,f5.1/)
           !f5.1,t80,f5.1,t88,f5.1,t96,f5.1,t105,f5.1/)

!
!   now the summary data for further processing elsewhere...
!
    write(lunDVF,200)TreeNo,Spp,dbhob,dbThickness,isForked,nReadings,   &
                     totLength, & !totNatCSVol,totAkimaVol,totHermiteVol,   &
                     LengthFlag, totConicVol
200 format(i5,1x,a4,2(1x,f5.1),1x,L1,1x,i5, 1x,f8.3, 1x,L1, 1x, f8.3)


!
!   store each tree in the segments arrays-- **Note that they are stored
!   *top down* as in the outpur files...
!
    begSeg = nTotReadings + 1
    endSeg = nTotReadings + nReadings
    Rtgrads(begSeg:endSeg) = TGRADS(nReadings:1:-1)
    Rfgrads(begSeg:endSeg) = FGRADS(nReadings:1:-1)
    Rsinelv(begSeg:endSeg) = SINELV(nReadings:1:-1)

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

write(kbout,"('Number of dendrometer readings process = ',i5)")nTotReadings
write(kbout,"('Number of segments process = ',i5)")nTotSegs
write(kbout,"('Number of trees process = ',i5)")nTrees
write(kbout,"('Number of trees with errors encountered = ',i5)")nTreeErrors


close(lunIn)
close(lunOut)
close(lunDVF)
close(lunDetail)
!stop    'Clean completion'
return
end subroutine splineVolume
