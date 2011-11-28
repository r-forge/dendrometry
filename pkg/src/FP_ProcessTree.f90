subroutine FP_ProcessTree()
!-----------------------------------------------------------------
!    Routine: FP_ProcessTree.f90 
!
!    This routine is used to process the desired measurements for
!    an individual tree. It takes the raw Barr & Stroud dendrometer
!    measurements and converts them to diameter, range and height
!    from horizontal. It then calculates dib's and section lengths.
!    Finally, cross sectional area (i.b.) is calculated for each
!    measurement and a natural cubic spline is fitted to the tree.
!    Using the coefficients of the spline, cubic foot volumes
!    are subsequently calculated.
!
!   This routine handles forked trees (repositioning of the instrument)
!   as long as the forks are all in the same group; that is, no nesting
!   of forks within forks is allowed currently. This could be added by
!   making ForkedAt, nForkedPoints a record structure and applying
!   the same set of code to each individual fork set. But not now!
!   Follow the STXMOD convention for handling volume and length of the
!   crotch section between the last stem reading before the fork and
!   the first fork readings if you want to include this volume and length.
!   If you do not want to include it you just do not repeat the lower
!   measurement on one of the forks. JHG (June 1998)
!
!   Note: It can happen that both the Akima and Hermite (one in the same!)
!         splines can both produce negative section volumes at the tip
!         of the tree (the last section) depending on where the measurements
!         were taken. This does not happen too often, but has occurred,
!         and the negative volumes are not very large (normally on the 
!         order of 0.2cuft). I have now put in the code a flag for this;
!         any negative section volume is set to 0.0 and a "-" flag is 
!         printed on the tree report to the right of that section volume
!         to let you know that this volume has been changed. (JHG July 1998)
!
!
!Author...                                    Date: 8-Dec-97
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
use TreeMod
use GrosenbaughMod
!use NRSplineMod
!use msimsl
!use HermiteSplineMod

implicit none
!locals...
logical*1   ForkMatch           !used to flag beginning of forked sections
character*1    ForkMark            !to mark forks on printouts
!character*1, pointer::        &
!            NegativeNCSV(:),  & !flags a negative section volume for NatCSVol
!            NegativeAV(:),    & !flags a negative section volume for AkimaVol
!            NegativeHV(:)       !flags a negative section volume for HermiteVol
character*2, pointer :: prod(:)    !product class/code w/o blanks for *.DVL/S+
integer    i,j
integer iSegBeg, iSegEnd        !beginning & ending measurement # for forked
                                !section processing of diameter, range, etc.
!
!    natural cubic spline-related variables (see Numerical Recipes)...
!
!logical*1 doIntegrate            !T: volumes calculated   
!real    x                        !a point for spline interpolation
!real    y                        !and its interpolated value
!real    yint                     !cf volume of the section containing x
!real, pointer::scndDeriv(:)        !solved spline coefficients    (NR)
!real, pointer::xa(:),ya(:)        !data input for spline (NR)

!
!    Akima spline-related variables (see IMSL)...
!
!real*8 a,b
!real*8, pointer::xdata(:),fdata(:),break(:),cscoef(:,:)

!
!    Hermite spline-related variables...
!
!real*8, pointer::b1(:),b2(:),b3(:)            !spline coefficients

!
!-----------------------------------------------------------------
!

real*8    kFactor    
common/jhgconstants/kFactor

!-----------------------------------------------------------------
!
!    allocate the storage required (from TreeMod)...
!
allocate(dob(nReadings),dib(nReadings),SegLength(nReadings),    &
         instRange(nReadings),HtAboveHorizontal(nReadings),     &
         xsArea(nReadings),                                     &
         NatCSVol(nReadings),AkimaVol(nReadings),HermiteVol(nReadings),    &
         ConicVol(nReadings),CumLength(nReadings) )

!
!   allocate local storage...
!
!allocate(NegativeNCSV(nReadings),NegativeAV(nReadings),NegativeHV(nReadings))

!
!    initialize these arrays (etc.) for the current tree...
!
dob = 0.0d0
dib = 0.0d0
SegLength = 0.0d0
CumLength = Seglength
instRange = 0.0d0
HtAboveHorizontal = 0.0d0
xsArea =  0.0d0
NatCSVol = 0.0d0
totNatCSVol = 0.0d0
AkimaVol = 0.0d0
totAkimaVol = 0.0d0
HermiteVol = 0.0d0
totHermiteVol = 0.0d0
ConicVol = 0.0d0
totConicVol = 0.0d0
totLength = 0.0d0
!NegativeNCSV = ' '
!NegativeAV = NegativeNCSV
!NegativeHV = NegativeNCSV

!
!    get the converted measurements from the Barr & Stroud...
!
if(isForked) then                           !either fork or moved instrument
    do i=1,nForkedPoints                    !do each component seperately   
        if(i.eq.1) then                     !main stem before first fork
            iSegBeg = 1
            iSegEnd = ForkedAt(i)                
        elseif(i.eq.nForkedPoints) then     !last fork to tip
            iSegBeg = ForkedAt(i)
            iSegEnd = nReadings                
        else                                !intermediate fork change
            iSegBeg = ForkedAt(i)
            iSegEnd = ForkedAt(i+1)-1            
        endif
        call jhgSBRD(iSegBeg,iSegEnd)        !convert segments
    enddo
else                                        !no forks--do entire tree
    call jhgSBRD(1,nReadings)                
endif

!
!    get dib and section lengths; if you want to include the volume of the
!   crotch at a fork (and its length) then you must repeat the measurement
!   of the last reading below the crotch on one of the forks--this is
!   standard STXMOD convention...
!
LengthFlag = .false.        !assume the best--all positive lengths        
do i=1,nReadings
    dib(i) = BarkModelTwo(dbhib,dbhob,dob(i))
    if(i.gt.1) then
!
!       check to see if we are at one of the beginning measurements of a
!       fork; if so, then the segment just below gets zero length, etc...
!
        if(isForked) then           !check further on forked trees
            ForkMatch = .false.
            do j=2,nForkedPoints    !skip last msrmnt before fork on main trunk
                if(i.eq.ForkedAt(j)) then
                    ForkMatch = .true.
                    exit
                endif
            enddo
            if(ForkMatch) then      !got a match, no length or volume here
                SegLength(i) = 0.0d0
             else                    !business as usual otherwise
                SegLength(i) = HtAboveHorizontal(i) - HtAboveHorizontal(i-1)
            endif
         else                        !no forks--it's easy
            SegLength(i) = HtAboveHorizontal(i) - HtAboveHorizontal(i-1)
        endif
        CumLength(i) = CumLength(i-1) + SegLength(i)
    endif
!
!   check for negative section length which will mess things up...
!
    if(SegLength(i).lt.0.0) then
        write(lunOut,25)i-1, SegLength(i)
25      format('****>Negative section length encountered on section #', &
               i2,'; Section length = ',f10.5/                          &
               '     Processing of tree completed.')
        nTreeErrors = nTreeErrors + 1
        LengthFlag = .true. 
        return
    endif
enddo
totLength = sum(SegLength)

!
!    calculate cross sectional area for each diameter inside bark...
!
xsArea =  kFactor*dib*dib

!
!    now calculate the spline coefficients based on cross sectional area;
!   obviously, splines can not be calculated for forked trees (unless, of
!   course we somehow do it piecemeal!)...
!
!if(.not.isForked) then
!
!   first the natural cubic spline from Numerical Recipes...
!
  !  allocate(xa(nReadings),ya(nReadings),scndDeriv(nReadings))
  !  xa = sngl(CumLength)                !x input points
  !  ya = sngl(xsArea)                    !y input points
  !  call spline(xa,ya,nReadings,0.0,0.0,scndDeriv)  !get second derivatives

!
!    and, using the midway point in each section, numerically integrate
!    the section for volume in cf...
!
  !  doIntegrate = .true.
  !  do i=2,nReadings
  !      x = sngl((CumLength(i) + CumLength(i-1))/2.0d0)
  !      call splintGQ(xa,ya,scndDeriv,nReadings,x,y,yint,doIntegrate)
  !      NatCSVol(i) = dble(yint)
  !      if(NatCSVol(i).lt.0.0d0) then
  !          NatCSVol(i) = 0.0d0
  !          NegativeNCSV(i) = '-'       !flag this as having been negative
  !      endif
  !  enddo
!    deallocate(xa,ya,scndDeriv)
  !  totNatCSVol = sum(NatCSVol)



!
!    IMSL spline routine: Akima method; also Harry's Hermite method...
!
  !  allocate(xdata(nReadings),fdata(nReadings),break(nReadings),    &
  !          cscoef(4,nReadings))                        !for Akima
  !  fdata = xsArea
  !  xdata = CumLength
  !  call dcsakm(nReadings,xdata,fdata,break,cscoef)        !IMSL Akima
  !  allocate(b1(nReadings),b2(nReadings),b3(nReadings))    !for Hermite
  !  call hermid(nReadings,xdata,fdata,b1,b2,b3)            !Hermite
  !  do i=2,nReadings
  !      a = CumLength(i-1)
  !      b = CumLength(i)
  !      AkimaVol(i) = dcsitg(a,b,nReadings,break,cscoef)
  !      if(AkimaVol(i).lt.0.0d0) then
  !          AkimaVol(i) = 0.0d0
  !          NegativeAV(i) = '-'       !flag this as having been negative
  !      endif
!
!        HermiteVol(i) = (b-a)*(xsArea(i-1) + (b-a)*(b1(i-1)/2.0d0 +        &
!                        (b-a)*(b2(i-1)/3.0d0 + (b-a)*b3(i-1)/4.0d0)))
    !    HermiteVol(i) = IntegrateHermite(i,a,b,b1,b2,b3,xsArea(i-1))
    !    if(HermiteVol(i).lt.0.0d0) then
    !        HermiteVol(i) = 0.0d0
    !        NegativeHV(i) = '-'       !flag this as having been negative
    !    endif
    !enddo
!    deallocate(xdata,fdata,break) !,cscoef)
!    deallocate(b1,b2,b3)
   ! totAkimaVol = sum(AkimaVol)
   ! totHermiteVol = sum(HermiteVol)
!endif

!
!    compute cubic foot volume via conic sections...
!
do i=2,nReadings
    if(tgrads(i).le.0.0d0 .and. fgrads(i).le.0.0d0 .or. i.eq.nReadings) then
        ConicVol(i) = SegLength(i)*xsArea(i-1)/3.0d0    !cone
    else                                                !frustrum of cone
        ConicVol(i) = SegLength(i)*(dib(i-1)*dib(i-1) +             &
                      dib(i-1)*dib(i) +  dib(i)*dib(i))*kFactor/3.0d0
    endif
enddo
totConicVol = sum(ConicVol)


!
!    output the results...
!
write(lunOut,50)
50 format(t2,6('-'),'Instrument Readings',5('-'),t40,           &                 
          'Instr',t48,'HtAbv',t63,'Total',t71,'Segmnt',         &
          !t80,7('-'),'Cubic Foot Volume',7('-')/                &
          t80, 'Conic'/                                         &
          t5,'tgrad',t15,'fgrad',t25,'sinelv',t34,'prod',       &
          t40,'Range',t48,'Horiz',t57,'dib',t63,'Length',       &
          t71,'Length',                                         &
          !t80,'NatCS',t88,'Akima', t95,'Hermite',               &
          !t105,'Conic'/111('-'))
          t80,'Volume'/90('-'))


allocate(prod(nReadings))
do i=nReadings,1,-1
!
!   first, set the mark character on the report for a forked point...
!
    if(isForked) then
        ForkMatch = .false.
        do j=2,nForkedPoints    !skip last msrmnt before fork on main trunk
            if(i.eq.ForkedAt(j)) then
                ForkMatch = .true.
                exit
            endif
        enddo
        if(ForkMatch) then
            ForkMark = '+'                      !mark the fork
        else
            ForkMark = ' '                      !no fork here
        endif
    else
        ForkMark = ' '                          !normal measurement
    endif
!
!   then print the report...
!
    write(lunOut,100)TGRADS(i),FGRADS(i),SINELV(i),GAMATH(i),ForkMark,    &
                     instRange(i),HtAboveHorizontal(i),dib(i),            &
                     CumLength(i),SegLength(i),                           &
                     !NatCSVol(i),NegativeNCSV(i),                         &
                     !AkimaVol(i),NegativeAV(i),                           &
                     !HermiteVol(i),NegativeHV(i),                         &
                     ConicVol(i)
100 format(t5,f5.1,t15,f5.1,t25,f7.4,t35,a2,a1,t40,f5.1,t48,f5.1,t55,   &
           f5.1,t64,f5.1,t72,f5.1,                                      &
           !t80,f5.1,a1,t88,f5.1,a1,t96,f5.1,a1,                        &
           !t105,f5.1)
           t80,f5.1)

!
!   and save the detailed section information to its file...
!
    prod(i) = GAMATH(i)                     !assign product for output
    if(prod(i).eq."  ") prod(i) = "NA"      !no blanks for R input
    write(lunDetail,120)TreeNo,i,Spp, &    !dbhob,dbThickness,isForked,        &
                        dob(i),dib(i),prod(i),        &
                        !NatCSVol(i),
                        TGRADS(i),FGRADS(i),SINELV(i),       &
                        instRange(i),HtAboveHorizontal(i),   &
                        SegLength(i),CumLength(i),ConicVol(i)
                     
120 format(i4,1x,i2,1x,a4,1x,        & !f4.1,1x,f4.1,1x,l1,1x,
           2(f6.3,1x),a2, 2(1x,f7.3), 1x,f7.5, 5(1x,f7.3))

enddo
deallocate(prod)




!
!   deallocate local storage; spline functions may not be allocated for
!   forked trees so check before deallocating...
!
!deallocate(NegativeNCSV,NegativeAV,NegativeHV)      !(-) volume flags
!Natural cubic spline...
!if(associated(xa)) deallocate(xa)                   
!if(associated(ya)) deallocate(ya)                   
!if(associated(scndDeriv)) deallocate(scndDeriv) 
!Akima cubic spline...
!if(associated(xdata)) deallocate(xdata)
!if(associated(fdata)) deallocate(fdata)
!if(associated(break)) deallocate(break)
!if(associated(cscoef)) deallocate(cscoef)
!Hermite cubic spline...
!if(associated(b1)) deallocate(b1)
!if(associated(b2)) deallocate(b2)
!if(associated(b3)) deallocate(b3)


return
end subroutine FP_ProcessTree
