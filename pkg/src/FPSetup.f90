      SUBROUTINE fpsetup
!======================================================================*
!--------------- P r o g r a m   D e s c r i p t i o n ----------------*
!======================================================================*
!
!     {FPSetup.for}
!
!     This routine simply initializes the constants for the Barr and
!     Stroud dendrometer models FP-12 and FP-15. It was "lifted" out
!     of the STXMOD program (after some bit of searching!).
!
! ***>This routine should be called before any calls to FPDRE are
!     made to calculate diameters, ranges or heights.
!
!----Author...
!     Jeffrey H. Gove             (603) 868-5692 or 834-5798 (FTS)
!     USDA Forest Service
!     Forestry Sciences Lab       14-March-91
!     P.O. Box 640
!     Durham, N.H. 03824
!
!======================================================================*
!-------------- V a r i a b l e    D i c t i o n a r y ----------------*
!======================================================================*
!
!-->Variables...
!
!     Name                   Description
!   --------          ---------------------------------------------
!     B                 8" Barr and Stroud instrument base
!     half              1/2 !!
!
!     **all other variables are from Grosenbaugh--good luck!!
!
!======================================================================*
!------- S u b r o u t i n e s   &   F u n c t i o n s   U s e d ------*
!======================================================================*
!
!     Name                   Description/Library
!   --------          ---------------------------------------------
!     dsin              FORTRAN-intrinsic
!     dcos              FORTRAN-intrinsic
!     dsqrt             FORTRAN-intrinsic
!
!======================================================================*
!--------------- V a r i a b l e   D e c l a r a t i o n  -------------*
!======================================================================*
!
      IMPLICIT NONE
      Real*8   gg,cg,w,z,q,b,u,g,half
!
!======================================================================*
!------------- C o m m o n   A r e a s --------------------------------*
!======================================================================*
!
!     all machine constants passed in the following common...
!
      common/lrg/gg,cg,w,z,q,b,u,g
!
!======================================================================*
!----------- V a r i a b l e   I n i t i a l i z a t i o n  -----------*
!======================================================================*
!
      DATA  Q/.1964673d-1/,U/-1.1905d0/,G/1.5658d0/,half/0.50d0/
      data  b/8.0d0/
!
!======================================================================*
!-------------------  B e g i n   t h e   R o u t i n e  --------------*
!======================================================================*
!
      CG=.1745329251994329D-01
      W=half*DSIN(U*CG)
      Z=half*DCOS(U*CG)+half
      GG=1.0d0/(1.0d0+G**2-2.0d0*DSQRT(G**2-(G**2)*(Q**2)))
      RETURN
      END subroutine FPSetup
