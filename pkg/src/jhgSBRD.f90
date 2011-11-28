!-----------------------------------------------------------------
!    Routine: jhgSBRD.for 
!
!    This is simply Grosenbaugh's modified SBRD lifted from
!    STXMOD.
!
!    All caps are from LRG; lower are jhg modifications to sbrd.
!
!    I believe that the comments in "***" boxes are due either
!    to Jim Space or John Rennie.
!
!
!Author (Modifications)...                    Date: 9-Dec-97
!    Jeffrey H. Gove
!    USDA Forest Service
!    Northeastern Forest Experiment Station
!    P.O. Box 640
!    Durham, NH 03824
!    jhgove@christa.unh.edu
!    phone: 603-868-7667    fax: 603-868-7604
!-----------------------------------------------------------------
!
!
! **********************************************************************
! *  SUBROUTINE FOR CONVERTING SHORT-BASE-RANGEFINDER READINGS TO      *
! *  DIAMETER, ELEVATION, AND RANGE IN ST22.                           *
! *  UPDATED TO HANDLE LEANING TREES - SEE 1980 FOREST SCIENCE 26:207. *
! **********************************************************************
!
      SUBROUTINE jhgSBRD(bot,top)
    use TreeMod

    implicit none
!    arguments...
    integer, intent(in):: bot,top
!    locals...
      LOGICAL*1 NOTOTH
    integer l,l1,lll,m,m1,mmm,i
    real*8, pointer::X(:),Y(:),SK(:) 
    real*8, pointer::d(:),r(:),e(:)
    real*8    vl,tn,pze,tr,tm,rc    
!
!
    Real*8   gg,cg,w,z,q,b,u,g
    common/lrg/gg,cg,w,z,q,b,u,g
!
!
!    assign my variables to lrg's and allocate storage...
!
    m = top                    !last reading on stem section
    l  = bot                !first reading on stem section
    allocate(x(nReadings),y(nReadings),sk(nReadings), d(nReadings),r(nReadings),e(nReadings))

!
!    next assignments are from st22...
!
      m1 = m-1
    l1 = l + 1            
!
!    initialize...
!
      VL=1.0d0
      TN=0.0d0
      PZE=0.0d0
      NOTOTH=.TRUE.
      LLL=L
      MMM=M
!
! **********************************************************************
! *  IF LAST FGRADS IS NEGATIVE, LAST SINELV IS TO TREE TOP, OR        *
! *  PROJECTION DESIRED.  IF LAST TGRADS IS NEGATIVE, LAST SINELV      *
! *  IS TO TREE TOP.                                                   *
! **********************************************************************
!
      IF (FGRADS(M) .GE. PZE) GO TO 1        !projection to tip
      MMM=M1
      D(M)=0.0
      IF (TGRADS(M) .GE. PZE) GO TO 1        
      NOTOTH=.FALSE.                        !last sinelv is to tip
!
! **********************************************************************
! *  IF FIRST TGRADS IS NEGATIVE, FIRST 2 TRIOS ARE DIRECT             *
! *  MEASUREMENTS.                                                     *
! **********************************************************************
!
    1 IF (TGRADS(L) .GE. PZE) GO TO 2
      LLL=L1+1                            !skip direct measurements
    2 IF (METH .EQ. 0) GO TO 4
!
! **********************************************************************
! *  ADJUST SINELV FOR B & S MODEL 15                                  *
! **********************************************************************
!
      do I=LLL,MMM
        SINELV(I)=SINELV(I)-1.0
      enddo
!
! **********************************************************************
! *  CALCULATE RANGE (R(I)), HEIGHT FROM HORIZONTAL (E(I)), AND        *
! *  DIAMETER (D(I)), FOR EACH MEASUREMENT POINT (NEXT 17 STATEMENTS). *
! **********************************************************************
!
    4 do I=LLL,MMM
        SK(I) = 1.0d0/(1.0d0 - SINELV(I)**2)
        X(I) = Q*DCOS (TGRADS(I)*CG*.9D0)
        x(I) = G*X(I)*DSQRT(GG-GG**2*X(I)**2) - X(I)*DSQRT(GG-G**2*GG**2*X(I)**2) !xx
        x(I) = (X(I)*DSQRT(1.D0-X(I)**2)+W)/(-X(I)**2+Z)    !t
        x(I) = (2.D0*x(I))/(1.D0-x(I)**2)                    !ta
        x(I) = (DSQRT(SK(I)*x(I)**2+1.D0)-1.D0)/x(I)        !tt
        Y(I) = Q*DCOS(FGRADS(I)*CG*.9D0)
        y(I) = G*Y(I)*DSQRT(GG-GG**2*Y(I)**2) - Y(I)*DSQRT(GG-G**2*GG**2*Y(I)**2) !yy
        y(I) = (y(I)*DSQRT(1.D0-Y(I)**2)+W)/(-Y(I)**2+Z)    !f
        y(I) = (2.D0*y(I))/(1.D0-y(I)**2)                    !fa
        y(I) = (DSQRT(SK(I)*y(I)**2+1.D0)-1.D0)/y(I)        !ff
        R(I) = B*(SK(I)-x(I)*y(I))/(24.0d0*x(I))
        E(I) = R(I)*SINELV(I)
        D(I) = B*DABS((x(I)-y(I))/x(I))
    enddo

      IF (LLL .EQ. MMM)  GO TO  8
      VL = R(MMM)*DSQRT(1.0d0 - SINELV(MMM)**2) - R(LLL)* DSQRT(1.0d0 - SINELV(LLL)**2)
      TN = VL/(R(MMM)*SINELV(MMM)-R(LLL)*SINELV(LLL))
      VL = DSQRT(1.0d0 + TN*TN)

      do I=LLL,MMM
        E(I)=E(I)*VL
      enddo

    8 IF (TGRADS(L) .GE. PZE) GO TO 9        !skip direct if N/A
!
! **********************************************************************
! *  INCORPORATE FIRST 2 SETS OF MEASUREMENTS WHERE DIRECT (NEXT 6     *
! *  STATEMENTS).                                                      *
! **********************************************************************
!
      D(L1) = FGRADS(L1)                !fgrads==dbh here
      R(L1) = 0.001d0                    !unknown range
      E(L1) = E(LLL) - SINELV(L1)        !sinelv==direct length here
      D(L) = FGRADS(L)                !fgrads==dbh here
      R(L) = 0.001d0                    !unknown range
      E(L) = E(L1) - SINELV(L)        !sinelv==direct length here
!
! **********************************************************************
! *  CALCULATE DISTANCE FROM HORIZONTAL TO TREE TIP (NEXT 5 STATEMENTS)*
! **********************************************************************
!
    9 IF (NOTOTH) GO TO 10
      IF (METH .NE. 0) SINELV(M) = SINELV(M) - 1.0d0
      TR = SINELV(M1)/DSQRT(1.0d0 - SINELV(M1)**2)
      TM = SINELV(M)/DSQRT(1.0d0 - SINELV(M)**2)
      RC = R(M1)*SINELV(M1)/TR
      E(M) = E(M1)+RC*VL*(TM-TR)/(1.0d0 - TM*TN)
      R(M) = 0.001d0
      D(M) = 0.1d0

!
!    done, assign/save the computed quantities, etc...
!
   10 continue
    HtAboveHorizontal(bot:top) = E(bot:top)
    instRange(bot:top) = r(bot:top)
    dob(bot:top) = d(bot:top)
       deallocate(x,y,sk,d,r,e)
      RETURN
      END subroutine jhgsbrd
