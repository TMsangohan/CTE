! AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ***************************************************************************
! VERSION 1.0 : 
!   AUTHOR    : TOM MERTENS 
!   DATE      : 10/06/2015
!   COPYRIGHT : CERN
!
!   DESCRIPTION : 
!       CODE SPLIT IN SEPERATE FILES THAT WERE CONVERTED TO MODULES
!       CODE ADAPTED TO FORTRAN 95 STANDARD
!       VARIABLES ARE NOW DEFINED WITH A FIXED KIND IN ORDER TO ALLOW FOR A SUBROUTINE TO CHECK IF ON 
!       A DIFFERENT PLATFORM THE RANGES AND PRECISIONS CAN BE OBTAINED = FAILSAFE FOR NUMERICAL INCONSISTENCY
! ***************************************************************************

PROGRAM collider_time_evolution
  USE COMMONDATA
  USE RAN
  USE modgetinput
  IMPLICIT NONE 

! ------------------------------------------------------------------------------------------------------------------
! ******************* DESCRIPTION VARIABLES*************************************************************************
! ------------------------------------------------------------------------------------------------------------------
!
! *** INTEGER KIND I4B *********************************************************************************************
! nb                : CONSTANT DEFINED IN MODULECOMMON nb = 6000000
! kkk               : DUMMY LOOP VARIABLE
! np1               :
! np2               :
! np10              :
! np20              :
! nLostLum1         : NUMBER OF PARTICLES LOST IN BEAM 1 DUE TO LUMINOSITY DURING THE ACTIVE SIMULATION TURN
! nLostLum2         : NUBMER OF PARTICLES LOST IN BEAM 2 DUE TO LUMINOSITY DURING THE ACTIVE SIMULATION TURN
! nLostLumSum1      : TOTAL NUMBER OF PARTICLES LOST IN BEAM 1 DUE TO LUMINOSITY UP UNTIL THE ACTIVE SIMULATION TURN
! nLostLumSum2      : TOTAL NUMBER OF PARTICLES LOST IN BEAM 2 DUE TO LUMINOSITY UP UNTIL THE ACTIVE SIMULATION TURN
! nLostMom1         : 
! nLostBeta1        :
! nLostMomSum1      :
! nLostBetaSum1     :
! nLostMom2         :
! nLostBeta2        :
! nLostMomSum2      :
! nLostBetaSum2     :
! np                :
! iwrite            :
! kcheck            :
! nLostDebunch1     : NUMBER OF PARTICLES LOST IN BEAM 1 DUE TO DEBUNCHING DURING THE ACTIVE SIMULATION TURN
! nLostDebunch2     : NUMBER OF PARTICLES LOST IN BEAM 2 DUE TO DEBUNCHING DURING THE ACTIVE SIMULATION TURN
! nLostDebunchSum1  : TOTAL NUMBER OF PARTICLES LOST IN BEAM 1 DUE TO DEBUNCHING UP UNTIL THE ACTIVE SIMULATION TURN
! nLostDebunchSum2  : TOTAL NUMBER OF PARTICLES LOST IN BEAM 2 DUE TO DEBUNCHING UP UNTIL THE ACTIVE SIMULATION TURN
!
! *** REAL KIND SP *************************************************************************************************
! fmix1             :
! fmix2             :
! lumi              : LUMINOSITY
! redfac            : LUMINOSITY GEOMETRIC REDUCTION FACTOR
! fdum              :
! pnumber1          :
! emixy1            :
! tauhat1           :
! rmsDelta1         :
! rmsBunchLen1      :
! emix1             :
! emiy1             :
! ex1               : HORIZONTAL EMITTANCE BEAM 1 DYNAMICALLY UPDATED DURING RUN
! ex2               : HORIZONTAL EMITTANCE BEAM 2 DYNAMICALLY UPDATED DURING RUN
! ey1               : VERTICAL EMITTANCE BEAM 1 DYNAMICALLY UPDATED DURING RUN
! ey2               : VERTICAL EMITTANCE BEAM 2 DYNAMICALLY UPDATED DURING RUN
!
! *** REAL ARRAY KIND SP ******************************************************************************************
! x1                : ARRAY CONTAINING HORIZONTAL POSITIONS OF PARTICLES IN BEAM 1 - ? DYNAMICALLY UPDATED ?
! px1               : ARRAY CONTAINING HORIZONTAL MOMENTA OF PARTICLES IN BEAM 1 - ? DYNAMICALLY UPDATED ?
! t1                :
! pt1               :
! y1                : ARRAY CONTAINING VERTICAL POSITIONS OF PARTICLES IN BEAM 1 - ? DYNAMICALLY UPDATED ?
! py1               : ARRAY CONTAINING VERTICAL MOMENTA OF PARTICLES IN BEAM 1 - ? DYNAMICALLY UPDATED ?
! xpavg1            :
! ykeep1            :
! ykeep2            :
! hglassfacOld      : THE HOURGLASS FACTOR IF BETALEVELLING IS USED NEEDING A NEW VARIABLE 10 = nIPS WITH LEVELLING
! ------------------------------------------------------------------------------------------------------------------
! ******************* DESCRIPTION VARIABLES **** END ***************************************************************
! ------------------------------------------------------------------------------------------------------------------

  INTEGER(I1B) :: kkk

  INTEGER(I4B) :: np1,np2,np10,np20,nLostLum1,nLostLum2,nLostLumSum1,nLostLumSum2
  INTEGER(I4B) :: nLostMom1,nLostBeta1,nLostMomSum1,nLostBetaSum1,nLostMom2,nLostBeta2,nLostMomSum2,nLostBetaSum2,np
  INTEGER(I4B) :: iwrite,kcheck,nLostDebunch1,nLostDebunch2,nLostDebunchSum1,nLostDebunchSum2

  REAL(SP)     :: fmix1,fmix2,lumi,redfac,fdum
  REAL(SP)     :: pnumber1,emixy1,tauhat1,rmsDelta1,rmsBunchLen1,emix1,emiy1  
  REAL(SP)     :: ex1,ex2,ey1,ey2 
  REAL(SP)     :: pnumber2,emixy2,tauhat2,rmsDelta2,rmsBunchLen2,emix2,emiy2

  REAL(SP), DIMENSION(1:nb) :: x1,px1,t1,pt1,y1,py1,xpavg1  
  REAL(SP), DIMENSION(1:nb) :: avgline1,avglinex1,avgliney1

  REAL(SP), DIMENSION(1:nb) :: x2,px2,t2,pt2,y2,py2,xpavg2  
  REAL(SP), DIMENSION(1:nb) :: avgline2,avglinex2,avgliney2

  REAL(SP), DIMENSION(7,nb) :: ykeep1
  REAL(SP), DIMENSION(7,nb) :: ykeep2

  REAL(SP), DIMENSION(10)   :: hglassfacOld

! -------------------------------------------------------------------------------------------------------------
! ******************* INITIALIZATION **************************************************************************
! -------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
! INITIALIZATION OF BETA-LEVELLING ARRAY VARIABLE
! -----------------------------------------------   

  DO kkk=1, 10
    hglassfacOld(kkk) = 1
  ENDDO

! ----------------------------------
! INITIALIZING SOME OF THE VARIABLES 
! ----------------------------------
  nLostLum1         = 0
  nLostLumSum1      = 0
  nLostLum2         = 0
  nLostLumSum2      = 0
  nLostDebunch1     = 0
  nLostDebunch2     = 0
  nLostDebunchSum1  = 0
  nLostDebunchSum2  = 0
  nLostMom1         = 0
  nLostMomSum1      = 0
  nLostBeta1        = 0
  nLostBetaSum1     = 0
  nLostMom2         = 0
  nLostMomSum2      = 0
  nLostBeta2        = 0
  nLostBetaSum2     = 0

! -------------------------------------------------------------------------------------------------------------
! ******************* INITIALIZATION **** END *****************************************************************
! -------------------------------------------------------------------------------------------------------------

! -------------------------------------------------------------------------------------------------------------
!	******************* WRITING STARTUP INFORMATION TO THE SCREEN PART 1	***************************************
! -------------------------------------------------------------------------------------------------------------

WRITE(6,*)
WRITE(6,*)
WRITE(6,*) '*************************************'
WRITE(6,*) '*      COLLIDER TIME EVOLUTION      *'
WRITE(6,*) '*************************************'
WRITE(6,*)
WRITE(6,*) 'processing input file ...'

! ------------------------------------------------------------
! READING IN SIMULATION INPUT
! ------------------------------------------------------------
call getinput(pnumber1,pnumber2,emix1,emiy1,emix2,emiy2,tauhat1,tauhat2,rmsDelta1,rmsBunchLen1,rmsDelta2,rmsBunchLen2)

WRITE(6,*)
WRITE(6,*) 'Starting main loop ...'
WRITE(6,*)
WRITE(6,*) '     turn      np1        np2     ex1(m)        ey1(m)    &
            &        sigt1(s)        ex2(m)         ey2(m)         sigt2(s)   '
WRITE(6,*) '______________________________________________________&
            &_________________________________________________________________'

! -------------------------------------------------------------------------------------------------
! ****************** TESTING AREA *****************************************************************
! -------------------------------------------------------------------------------------------------
! -------------------------
! *** MODULERAN3 **********
! -------------------------
!
! TEST INTEGER TO ALLOW TESTING OF RAN3 FUNCTION
! INTEGER(I4B) :: DUMMYRAN
!
! DUMMYRAN=34234
! WRITE(*,*) RAN3(DUMMYRAN)
!
!                                       RESULT : OK
! --------------------------------------------------------------------------------------------------
! -----------------------------
! *** MODULEGETINPUT **********
! -----------------------------
!    call getinput(pnumber1,pnumber2,emix1,emiy1,emix2,emiy2,tauhat1,tauhat2,rmsDelta1,rmsBunchLen1,rmsDelta2,rmsBunchLen2)
!
!
!                                       RESULT : OK
! --------------------------------------------------------------------------------------------------
!
! ****************** TESTING AREA **** END *********************************************************
! --------------------------------------------------------------------------------------------------


 !  ! %%%%%%%%%%%%%%%%%%% generating the particle distributions to track  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !    call addparticles(np1,x1,px1,t1,pt1,y1,py1,np10,emix1,emiy1,tauhat1,rmsDelta1,rmsBunchLen1)
 !    call addparticles(np2,x2,px2,t2,pt2,y2,py2,np20,emix2,emiy2,tauhat2,rmsDelta2,rmsBunchLen2)

 !    ! zero out sums
 !    call zerokeep(ykeep1,avgline1,avglinex1,avgliney1,np1)
 !    call zerokeep(ykeep2,avgline2,avglinex2,avgliney2,np2)





!!    calculate transverse emittances
 !    call getemittance(x1,y1,px1,py1,np1,x2,y2,px2,py2,np2,ex1,ex2,ey1,ey2)
!    write(6,*) '_____________________________________________________________________________& 
 !        &__________________________________________'


 !    ! initial emittances
 !    ! write(6,*)
 !    ! write(6,*) 'initial emittances: '
 !    ! write(6,*) 'ex1=',ex1
 !    ! write(6,*) 'ey1=',ey1
 !    ! write(6,*) 'ex2=',ex2
 !    ! write(6,*) 'ey2=',ey2
 !    ! write(6,*)


!!write first turn
 !    kturns=0
 !    iwrite=1
 !    if (writeMountSwitch.eq.1) then
 !       call mountainr(np1,y1,py1,t1,pt1,x1,px1,1,np10,pnumber1,avgline1)
 !       call mountainr(np2,y2,py2,t2,pt2,x2,px2,2,np20,pnumber2,avgline2)
 !       call mountTransv(avglinex1,avgliney1,emix1,emiy1,1)
 !       call mountTransv(avglinex2,avgliney2,emix2,emiy2,2)
 !    endif
 !    call writeEmi(ex1,ey1,t1,pt1,ex1,ey1,t2,pt2,np1,np2)
 !    call writeNb(np1,np2,np10,np20,pnumber1,pnumber2,nLostLum1,nLostLum2,nLostLumSum1,nLostLumSum2,&
 !         nLostDebunch1,nLostDebunch2,nLostDebunchSum1,nLostDebunchSum2,&
 !         nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1,nLostBeta2,nLostBetaSum2,nLostMom2,nLostMomSum2)

 !    if (writeAllCoordSwitch.eq.1) then   !     write all starting conditions to file
 !       call writecoord(np1,x1,px1,y1,py1,t1,pt1,0,1)
 !       call writecoord(np2,x2,px2,y2,py2,t2,pt2,0,2)
 !    endif

 !     Main loop over all turns
!!*******************************************************************************
 !    do kturns = 0,nturns
!!    logic for writes. Write output every nwrite turns
!!-----------------------------------------------------

 !       kcheck = kturns/nwrite
 !       kcheck = kcheck*nwrite
 !       if(kcheck.eq.kturns)then
 !          iwrite = 1
 !       else
 !          iwrite = 0
 !       endif
!!    always write the last turn
 !       if(kturns.eq.nturns)then
 !          write(6,*)' writing last turn'
 !          nwrite=-1
 !          iwrite=1
!!           writeAllCoordSwitch=1
 !       endif

!!    calculate transverse emittances
 !       call getemittance(x1,y1,px1,py1,np1,x2,y2,px2,py2,np2,ex1,ex2,ey1,ey2)

!!    Step through different physics processes:
!!-----------------------------------------------

!!    synchrotron motion
 !       fmix1= 1
 !       if(RFswitch.eq.1) then
!!           write(6,*) 'RF active'
 !          call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fmix1,nLostDebunch1,iwrite)
 !          call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fmix1,nLostDebunch2,iwrite)
 !       endif

!!    betatron motion
 !       if(betatronSwitch.eq.1) then
!!           write(6,*) 'betatron active'
 !          call sptransChrom(np1,x1,px1,y1,py1,pt1)
 !          call sptransChrom(np2,x2,px2,y2,py2,pt2)
 !       endif

 !       call keepTransvProf(avglinex1,avgliney1,x1,y1,np1,emix1,emiy1)
 !       call keepTransvProf(avglinex2,avgliney2,x2,y2,np2,emix2,emiy2)

!!    radiation damping and quantum excitation
 !       if(raddampSwitch.eq.1) then
!!           write(6,*) 'raddamp active'
 !          call raddamp(np1,x1,px1,y1,py1,pt1)
 !          call raddamp(np2,x2,px2,y2,py2,pt2)
 !       endif

!!    ibs. normalize equivalent time to beam 1
 !       if(IBSswitch.eq.1) then
!!           write(6,*) 'IBS active'
 !         call ibslong(np1,y1,py1,t1,pt1,x1,px1,ex1,ey1,np10,pnumber1,pnumber1,avgline1,iwrite)
 !         call ibslong(np2,y2,py2,t2,pt2,x2,px2,ex2,ey2,np20,pnumber2,pnumber1,avgline2,iwrite)
 !      endif

!!    collimation/aperture
 !      if(collimationSwitch.eq.1) then
 !         call collimation(np1,y1,py1,t1,pt1,x1,px1,np10,pnumber1,nLostMom1,nLostBeta1,1)
 !         call collimation(np2,y2,py2,t2,pt2,x2,px2,np20,pnumber2,nLostMom2,nLostBeta2,2)
 !      endif

 !      if(blowupSwitch.eq.1) then
 !         write(101,*) kturns,pxKickFac,pyKickFac
 !         call blowup(np1,py1,px1)
 !         call blowup(np2,py2,px2)
 !      endif

!!      collision between beam 1 and 2
 !      if(collisionSwitch.eq.1) then
!!         call writecoord(np1,x1,px1,y1,py1,t1,pt1,0,1)
!!         call writecoord(np2,x2,px2,y2,py2,t2,pt2,0,2)
 !         if(collRoutine.eq.'6a') then
 !            call collision6a(np1,x1,px1,t1,pt1,y1,py1,np10,pnumber1,ex1,ey1,&
 !                             np2,x2,px2,t2,pt2,y2,py2,np20,pnumber2,ex2,ey2,&
 !                             nLostLum1,nLostLum2,iwrite,lumi,redfac,hglassfacOld)
!!          add rf update to get rid of high amplitude particles in the hope this will help with bunch length growth
!!              fdum =0
!!             call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fdum,nLostDebunch1,iwrite)
!!             call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fdum,nLostDebunch2,iwrite)
 !         elseif(collRoutine.eq.'1d') then
 !            call collision1d(np1,x1,px1,t1,pt1,y1,py1,np10,pnumber1,&
 !   &             np2,x2,px2,t2,pt2,y2,py2,np20,pnumber2,nLostLum1,nLostLum2,iwrite,lumi,redfac,emixy1)
 !         else
 !            write(6,*) 'unkonwn collision routine'
 !            stop
 !         endif
!!         call writecoord(np1,x1,px1,y1,py1,t1,pt1,1,1)
!!         call writecoord(np2,x2,px2,y2,py2,t2,pt2,1,2)
 !      endif

 !      call keepturn(np1,y1,py1,ykeep1)
 !      call keepturn(np2,y2,py2,ykeep2)

!!    first turn already written
 !      if (kturns.eq.0) then
 !         iwrite=0
 !         call writeLumi(lumi,redfac)
 !      endif
!!    write output if desired
 !      if (iwrite.eq.1) then
!!           call writemomentsshort(np1,y1,py1,t1,pt1,x1,px1,1)
!!           call writemomentsshort(np2,y2,py2,t2,pt2,x2,px2,2)
 !         if (writeAllCoordSwitch.eq.1) then   !     write all coordinates to file
 !            call writecoord(np1,x1,px1,y1,py1,t1,pt1,kturns,1)
 !            call writecoord(np2,x2,px2,y2,py2,t2,pt2,kturns,2)
 !         endif
 !         if (writeMountSwitch.eq.1) then
 !            call mountainr(np1,y1,py1,t1,pt1,x1,px1,1,np10,pnumber1,avgline1)
 !            call mountainr(np2,y2,py2,t2,pt2,x2,px2,2,np20,pnumber2,avgline2)
 !            call mountTransv(avglinex1,avgliney1,emix1,emiy1,1)
 !            call mountTransv(avglinex2,avgliney2,emix2,emiy2,2)
 !         endif
 !         call writeEmi(ex1,ey1,t1,pt1,ex2,ey2,t2,pt2,np1,np2)
 !         call writeNb(np1,np2,np10,np20,pnumber1,pnumber2,nLostLum1,nLostLum2,nLostLumSum1,nLostLumSum2,&
 !              nLostDebunch1,nLostDebunch2,nLostDebunchSum1,nLostDebunchSum2,&
 !         nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1,nLostBeta2,nLostBetaSum2,nLostMom2,nLostMomSum2)
 !         if(collisionSwitch.eq.1) call writeLumi(lumi,redfac)
 !      endif

 !   enddo
 !  


    write(6,*)
    write(6,*) 'COLLIDER TIME EVOLUTION: run finished normally'
    write(6,*) '**********************************************'
    
STOP

END PROGRAM collider_time_evolution


!*********************************************************************************
!*********************************************************************************
subroutine keepTransvProf(avglinex,avgliney,x,y,np,emix,emiy)
! bin the transverse coord., store the  profile in arrays
 !     include 'collider_time_evolution_common.f90'
  use commondata
  real, dimension(1:nb) :: x,y,avglinex,avgliney
  real :: xmn,dx,emix,ymn,dy,emiy
  integer :: k,np,nk

      xmn=-5*sqrt(betax*emix)
      ymn=-5*sqrt(betay*emiy)
      dx =-2*xmn/nBins
      dy=-2*ymn/nBins

      do k=1,np
! use nearest integer to find bin in x
         nk = nint((x(k)-xmn)/dx)
         if(nk.le.1)nk=1
         if(nk.gt.nBins)nk=nBins
         avglinex(nk)=avglinex(nk)+1

! same in y
         nk = nint((y(k)-ymn)/dy) ! using DIFFERENT scale in x and y
         if(nk.le.1)nk=1
         if(nk.gt.nBins)nk=nBins
         avgliney(nk)=avgliney(nk)+1
      enddo

  return
end


!*********************************************************************************
!*********************************************************************************

subroutine keepturn(np,y,py,ykeep)
!     include 'collider_time_evolution_common.f90'
  use commondata
  integer :: k,np
  real, dimension(1:nb) :: y,py
  real, dimension(7,nb) :: ykeep
  real :: yk,pk

! accumulates the running sums
      do k=1,np
         yk = y(k)
         pk = py(k)
         ykeep(1,k) = ykeep(1,k) + yk*yk + pk*pk
       enddo
!      write(60,*)ykeep(1,1),ykeep(1,2)
  return
end

!*********************************************************************************
!*********************************************************************************
subroutine sptransChrom(np,x,px,y,py,pt)
! full turn update for x, horizontal
  !    include 'collider_time_evolution_common.f90'
  use commondata
  real, dimension(1:nb) :: x,px,y,py,pt
  real :: psix,psiy,a11,a12,tpdqmin,xk,pkx,xk1,pxk1,yk,pky,yk1,pyk1,coeffchromx,coeffchromy,psi1,ptk
  integer :: np,k,m

      psix = 2*pi*tunex
      psiy = 2*pi*tuney
      coeffchromx = chromx*2*pi/(gamma0-1/gamma0)
      coeffchromy = chromy*2*pi/(gamma0-1/gamma0)

      tpdqmin = 2*pi*dqmin

      do k = 1,np

! rotate in y-py plane
         ptk=pt(k)
         psi1 = psiy + ptk*coeffchromy
         a11 = cos(psi1)
         a12 = sin(psi1)
         yk = y(k)
         pky = py(k)
         yk1  = yk*a11 + pky*a12
         pyk1 = pky*a11 - yk*a12
         y(k)=yk1

! rotate in x-px plane
         psi1 = psix + ptk*coeffchromx
         a11 = cos(psi1)
         a12 = sin(psi1)
         xk = x(k)
         pkx = px(k)
         xk1  = xk*a11 + pkx*a12
         pxk1 = pkx*a11 - xk*a12
         x(k)=xk1

! now have dqmin part - coupling between x and y
         px(k)=pxk1+tpdqmin*yk1
         py(k)=pyk1+tpdqmin*xk1

! thin sextupole kick 
         px(k)= px(k) + 0.5 * k2L * (xk1**2 - yk1**2) - k2Lskew * (xk1 * yk1)
         py(k)= py(k) + k2L * (xk1 * yk1) + 0.5 * k2Lskew * (xk1**2 - yk1**2)

      enddo
!      write(44,*)psi1,ptk*coeffchrom,ptk
  return
end
