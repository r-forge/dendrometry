    real*8 function BarkModelTwo(dbhib,dbhob,dob)
!-----------------------------------------------------------------
!    Routine: BarkModelTwo.for 
!
!    This is simply Grosenbaugh's second bark model lifted from
!    STXMOD.
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

    implicit none
    real*8, intent(in)::dbhib                !dbh inside bark
    real*8, intent(in)::dbhob                !dbh outside bark
    real*8, intent(in)::dob                    !diameter o.b. at point



    BarkModelTwo = dob*(1.0d0 - (1.0d0 - dbhib/dbhob)*(1.0d0/(2.0d0 - dob/dbhob)))

    end function BarkModelTwo
