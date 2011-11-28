!-----------------------------------------------------------------
!
!Author...                                    Date: 2-Nov-2011
!    Jeffrey H. Gove
!    USDA Forest Service
!    Northern Research Station
!    270 Mast Road
!    Durham, NH 03824
!    jhgove@christa.unh.edu
!    phone: 603-868-7667    fax: 603-868-7604
!-----------------------------------------------------------------
!
module ioMod

!
!input/output-related variables...
!
logical*1       Finished        !EOF 9999 reached for input file?
integer    ::   lunIn=20        !input file LUN for STX data
integer    ::   lunOut=22       !output file LUN for report
integer    ::   lunDVF=24       !output file LUN for dendrometry vol
integer    ::   lunDetail=26    !output file LUN for details
integer    ::   kbin=5          !keyboard in lun
integer    ::   kbout=6         !keyboard out lun
integer         nCtlRecs        !!# of control records: 0<=x<=10
integer         nLinesPP        !# of output records per page

contains

end Module ioMod
