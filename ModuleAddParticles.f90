! ---------------------------------------------------------------------------------------------------------------
! AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 1.0 : 
!   AUTHOR    : TOM MERTENS
!   DATE      : 17/06/2015
!   COPYRIGHT : CERN
!
!   DESCRIPTION : 
!       FUNCTION TO GENERATE PARTICLE DISTRIBUTIONS FOR THE SIMULATION
! ---------------------------------------------------------------------------------------------------------------
MODULE ModAddparticles
  USE COMMONDATA
  USE RAN
 
  IMPLICIT NONE

INTERFACE addparticles
  MODULE PROCEDURE addparticles
END INTERFACE

INTERFACE sampleLongMatched
  MODULE PROCEDURE sampleLongMatched
END INTERFACE

CONTAINS 

  SUBROUTINE sampleLongMatched(np,t,pt,ham1sig,sigs,tauhat,ptmax,v00,omega0,hammax)

  ! ------------------------------------------------------------------------------------------------------------------
  ! ******************* DESCRIPTION VARIABLES ************************************************************************
  ! ------------------------------------------------------------------------------------------------------------------
  ! 
  ! *** INTEGER KIND I4B *********************************************************************************************
  ! np          : NUMBER OF MACRO PARTICLES
  !
  ! *** REAL KIND SP *************************************************************************************************
  ! tk        : UNIFORMLY DISTRIBUTED 2 TAUHAT
  ! ptk       : UNIFORMLY DISTRIBUTED 2 PTMAX
  ! ham1sig   : INPUT 1 SIGMA HAMILTONIAN
  ! sigs      : INPUT LONGITUDINAL SIGMA
  ! tauhat    : HALF LENGTH OF THE BUCKET
  ! p1        : UNIFORMLY DISTRIBUTED PHASE RF SYSTEM 1
  ! p2        : UNIFORMLY DISTRIBUTED PHASE RF SYSTEM 2
  ! ptmax     : INPUT MAXIMUM PT
  ! prob      : 
  ! v00       :
  ! test      :
  ! omega0    : INPUT ANGULAR FREQUENCY
  ! ran3      :
  ! ham       :
  ! hammax    :
  !
  ! *** REAL KIND SP  DIMENSION np ***********************************************************************************
  ! t         : INPUT TIME COORDINATE PARTICLES
  ! pt        : INPUT MOMENTUM COORDINATE PARTICLES
  !
  ! ******************* DESCRIPTION VARIABLES **** END ***************************************************************
  ! ------------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(I4B)  :: np

    REAL(SP)      :: tk,ptk,ham1sig,sigs,tauhat,p1,p2,ptmax,prob,v00,test,omega0,ham,hammax

    REAL(SP), DIMENSION(nb) :: t,pt

    
    sigs  = 0.

    DO np=1,nMacro
55      CONTINUE
        tk   = tauhat*(2*ran3(iseed)-1)  
        ptk  = ptmax*(2*ran3(iseed)-1)
        
        p1   = nharm*tk*omega0   ! RANDOM DISTRITUTION OF PHASE FIRST RF SYSTEM
        p2   = nharm2*tk*omega0  ! RANDOM DISTRIBUTION OF PHASE SECOND RF SYSTEM

        ! -------------------------
        ! EXACT HAMILTONIAN:
        ! -------------------------

        ham  = 0.5*ptk*ptk*tcoeff+(cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)

        ! -------------------------------------------------------------------------------------------------
        ! SMALL ANGLE APPROXIMATION:
        ! **************************
        ! ham=0.5*ptk*ptk*tcoeff - p1**2/2*vrf/(nharm*v00*omega0) - p2**2/2*vrf2/(nharm2*v00*omega0)
        ! -------------------------------------------------------------------------------------------------


        IF (ham > hammax) GO TO 55 ! RESTART SAMPLING IF OUTSIDE BUCKET

        prob = EXP(-ham/ham1sig) ! HAM1SIG IS USED AS INPUT IN THE SUBROUTINE
        test = RAN3(iseed)

        IF (prob<test) GO TO 55

        sigs = sigs+tk**2

        ! ------------------------------
        ! WRITING FOR DEBUGGING PURPOSES
        ! WRITE(88,*) tk,ptk,ham
        ! ------------------------------

        t(np)   = tk
        pt(np)  = ptk
      ENDDO

      sigs=clight*sqrt(sigs/nMacro)
  
  RETURN
  END SUBROUTINE sampleLongMatched


  SUBROUTINE fitEmittance(x,px,np,bx,Jmax,emit)

  ! ------------------------------------------------------------------------------------------------------------------
  ! ******************* DESCRIPTION VARIABLES ************************************************************************
  ! ------------------------------------------------------------------------------------------------------------------
  ! 
  ! *** INTEGER KIND I4B *********************************************************************************************
  ! i           : DUMMY INDEX
  ! bin         : BIN NUMBER FOR STATISTICS
  ! np          : NUMBER OF MACRO PARTICLES
  ! maxbin      : MAXIMUM NUMBER OF BINS
  !
  ! *** REAL KIND SP *************************************************************************************************
  !
  ! ------------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER(I4B)  :: i,np,bin,maxbin

    REAL(SP)      :: x(nb),px(nb),binVals(1000),Jx,Jmax,bx,sigma1,sigma0,getNewtonSigma,bin1,emit,binwidth
  
  ! calculate the betatron action for each particle. Then bin the action and fit the resulting bin values to an exponential. Not as robust as calculating backwards from standard deviation as in subroutine

  maxbin=1000
  binwidth=0.025 ! fraction of a sigma used as bin width

  do i=1,maxbin
     binVals(i)=0.0
  enddo
  do i=1,np
     Jx=(x(i)**2+px(i)**2)/(2*bx*refEmxy)      ! calculate the normalized action.
!     write(35,*) Jx
     bin=int(Jx/binwidth)+1
     if (bin<=maxbin) binVals(bin)=binVals(bin)+1.
  enddo
  bin1=binVals(1)
  do i=1,maxbin ! for better numerical accuracy, normalize everything to first bin
     binVals(i)=binVals(i)/bin1
!     write(55,*) binVals(i)
  enddo

! solve f'=0 with Newton's method -- see function getNewtonSigma for details
  sigma0 = 1.
  do i=1,100
     sigma1=getNewtonSigma(sigma0,binVals,Jmax,binwidth)
!     write(6,*) sigma0,sigma1
     if (abs(sigma1-sigma0).le.0.01) exit
     sigma0=sigma1
  enddo
  emit = sigma1 * refEmxy

  return
end subroutine fitEmittance

  SUBROUTINE addparticles(np,x,px,t,pt,y,py,np0,emix,emiy,tauhat,rmsDelta,rmsBunchLen)

  ! ------------------------------------------------------------------------------------------------------------------
  ! ******************* DESCRIPTION VARIABLES ************************************************************************
  ! ------------------------------------------------------------------------------------------------------------------
  ! 
  ! *** INTEGER KIND I4B *********************************************************************************************
  ! i           :
  ! np          :
  ! np0         :
  ! k           :
  ! j           :
  !
  ! *** REAL KIND AP *************************************************************************************************
  ! emix        :
  ! emiy        :
  ! tauhat      : HALF WIDTH OF THE BUNCH LENGTH EXPRESSED IN TIME (ns)
  ! dum1        :
  ! dum2        :
  ! dum3        :
  ! dum4        :
  ! rmsDelta    :
  ! rmsBunchLen :
  ! testemit    :
  ! sigv        :
  ! v00         :
  ! omega0      :
  ! ham         :
  ! hammax      :
  ! hamzero     :
  ! amp         :
  ! facc        :
  ! p1          :
  ! p2          :
  ! pcoeff      :
  ! phik        :
  ! pcoeff3     :
  ! phizero     :
  ! phik0       :
  ! prob        :
  ! ptk         :
  ! ptmax       : MAXIMUM ACCEPTED PT BY THE RF BUCKET
  ! radharm2    :
  ! r1          :
  ! r2          :
  ! vk          :
  ! vzero       :
  ! tk          :
  ! ran3        :
  ! test        :
  ! ampx        :
  ! ampy        :
  ! ampt        :
  ! ampPt       :
  ! epsl        :
  ! sigs0       :
  ! sigs1       :
  ! ham1sig0    :
  ! ham1sig1    :
  ! fPrime      :
  !
  ! *** REAL DIMENSION(nb) KIND AP ***********************************************************************************
  ! NB = 6000000 SET IN MODULECOMMON 
  ! x           :
  ! px          :
  ! t           :
  ! pt          :
  ! y           :
  ! py          :
  !
  ! ------------------------------------------------------------------------------------------------------------------
  !
  ! *** REAL KIND AP *************************************************************************************************
  ! readfile1   : LOGICAL VARIABLE TO KEEP TRACK OF BEAM 1 AND BEAM 2
  !
  ! ******************* DESCRIPTION VARIABLES **** END ***************************************************************
  ! ------------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE
    
    INTEGER(I4B)  ::  i
    INTEGER(I4B)  ::  np,np0,k,j

    REAL(SP)      ::  emix,emiy,tauhat,dum1,dum2,dum3,dum4,rmsDelta,rmsBunchLen,testemit
    REAL(SP)      ::  sigv,v00,omega0,ham,hammax,hamzero,amp,facc,p1,p2,pcoeff2,&
                      phik,pcoeff3,phizero,phik0,prob,ptk,ptmax,radharm2,r1,r2,vk,&
                      vzero,tk,ran3,test,ampx,ampy,ampt,ampPt,&
                      epsl,sigs0,sigs1,ham1sig0,ham1sig1,fPrime


    REAL(SP), DIMENSION(nb) :: x,px,t,pt,y,py
     
    LOGICAL       :: readfile1,readfile1Transv 
    SAVE readfile1, readfile1Transv

! -------------------------------------------------------------------------------------------------------------
! ******************* INITIALIZATION **************************************************************************
! -------------------------------------------------------------------------------------------------------------
    DATA readfile1 /.true./
    DATA readfile1Transv /.true./

 
! WHAT IS THE USE OF THIS ? EXCEPT FOR THE USE OF THE PHI ZERO ?      

! FIRST GET THE ZERO OF THE VOLTAGE
    phizero   = REAL(nharm)/REAL(nharm2)*3*pi/4
    vzero     = ABS(vrf)

    DO k=-1000,1000
       phik   = REAL(nharm)/REAL(nharm2)*(pi+k*pi/2000)
       vk     = vrf * SIN(phik) + vrf2 * SIN(nharm2*phik/nharm)

      IF(ABS(vk) < vzero)THEN
          vzero   = ABS(vk)
          phizero = phik
       ENDIF
    
    ENDDO
      
    phik0 = phizero
    DO k=-1000,1000
      phik = phik0 + k/100000.
      vk   = vrf * SIN(phik) + vrf2 * SIN(nharm2*phik/nharm)

      IF (abs(vk) < vzero) THEN
        vzero   = abs(vk)
        phizero = phik
      ENDIF
    ENDDO

    np          = 0
    pcoeff2     = pcoeff/omegarf
    pcoeff3     = pcoeff2/omegarf
    radharm2    = vrf2/vrf
    v00         = aatom * 938.e6/qatom ! VOLTAGE FOR UNIT CHARGE IN GAMMA
    sigv        = 0
    omega0      = 2 * pi/trev
    p1          = nharm*tauhat*omega0 ! tauhat IS HALF THE TIME WIDTH OF THE RF BUCKET
    p2          = nharm2*tauhat*omega0
    hammax      = (COS(p1)-1)*vrf/(nharm*v00*omega0) + (COS(p2)-1)*vrf2/(nharm2*v00*omega0) ! MAX HAMILTONIAN
    p2          = phizero*nharm2/nharm
    HAMZERO     = (COS(phizero)-1)*vrf/(nharm*V00*omega0) + (COS(P2)-1)*vrf2/(nharm2*V00*omega0) ! ZERO HAMILTONIAN
    ptmax       = SQRT(2*hammax/tcoeff) ! MAXIMALLY ALLOWED PT TO FALL INTO THE RF BUCKET DEFINED BY THE HAMILTONIAN

! GET LONGITUDINAL COORDINATES
! IF SWITCH IS SET TO
! 0: PARABOLIC DISTRIBUTION WITH SMOKE RING
! 1: READ FROM FILE
! 2: GAUSSIAN

    WRITE(6,*) 'generating longitudinal coordinates'
    SELECT CASE (longCoordMethod)
      CASE (0)    ! SMOKE RING DIST.
          DO np = 1,nMacro

 43          CONTINUE
             tk   = tauhat * (2*RAN3(iseed)-1)    ! RANDOMIZE A VALUE AROUND TAUHAT
             ptk  = ptmax  * (2*RAN3(iseed)-1)    ! RANDOMIZE A VALUE AROUND PTMAX

             ! -----------------------------------------
             ! ******** REMINDER (ModuleGetinput) ******
             ! -----------------------------------------
             !
             ! pcoeff = qatom*v1*omegarf/(938.e6*aatom)
             ! tcoeff = trev*eta/(beta**2*gamma0)
             ! rcoeff = pcoeff/tcoeff
             !
             ! -----------------------------------------
             !
             ! **********************************       DISTRIBUTION OF PHASE      ****************************************************
             !
             ! |--------------------------------------------------- nharm = 4 ------------------------------------|
             !
             !        +                         +                       +                       +                       +
             !     +     +                   +     +                 +     +                 +     +                 +     +
             !   +         +               +         +             +         +             +         +             +         +
             !  *-------------------------*-----------------------*-----------------------*-----------------------*------------------
             !               +          +              +         +             +         +             +         +
             !                 +      +                  +     +                 +     +                 +     +
             !                     +                        +                       +                       +    
             !  |---------- tk*omega0-----|
             !
             !  |------------------------------------ nharm*tk*omega0 --------------------------------------------|
             !
             ! **************************************************************************************************************************

             p1   = nharm*tk*omega0       ! tauhat IS HALF THE TIME WIDTH OF THE RF BUCKET, HERE A RANDOM P IS GENERATED AROUND THIS VALUE
             p2   = nharm2*tk*omega0

            ! ------------------------------------------------------------
            ! LONGITUINAL HAMILTONIAN WITH THESE RANDOMLY GENERATED VALUES
            ! ------------------------------------------------------------
             ham  = 0.5*ptk*ptk*tcoeff + (COS(p1)-1)*vrf/(nharm*v00*omega0) + (COS(p2)-1)*vrf2/(nharm2*v00*omega0)  
             
             IF (ham >= hammax) GO TO 43  ! IF HAMILTONIAN LARGER THAN THE ONE ACCEPTED BY THE BUCKET REJECT AND RE-GENERATE RANDOM PARTICLE

             prob = 1-(ham/hammax)
             prob = prob**power           ! POWER : BUNCH SHAPE PARAMETER, ONLY USED WITH SMOKE RING (GIVEN IN THE INPUT FILE)
             
             IF ((ABS(p1)>=phizero).and.(ham<=hamzero)) THEN
                prob = prob*(ham/hamzero)**alint  ! ALINT IS SMOKE RING PARAMETER (GIVEN IN THE INPUT FILE)
             ENDIF

             test = RAN3(iseed)

             IF (pro <= test) GO TO 43

             pt(np) = ptk
             t(np)  = tk

          ENDDO

       CASE (1)                                   ! READ FROM FILE
          IF (readfile1) THEN                     ! LOGICAL VARIABLE TO KEEP TRACK OF BEAM 1 AND BEAM 2
             OPEN(unit=10,file=trim(filebeam1))
             readfile1=.false.
             WRITE(6,*) 'reading longitudinal starting coordinates for beam1 from ', filebeam1
          ELSE
             OPEN(unit=10,file=trim(filebeam2))
             WRITE(6,*) 'reading longitudinal starting coordinates for beam2 from ', filebeam2
          ENDIF

          DO np = 1,nMacro
             READ(10,*) dum1,dum2,dum3,dum4,t(np),pt(np)
          ENDDO

          CLOSE(10)

       CASE (2) ! GENERATE BI-GAUSSIAN (POSSIBLY UNMATCHED, CAN ONLY BE MATCHED IN SMALL ANGLE APPROXIMATION)
          ampPt =   gamma0*rmsDelta
          ampt  =   rmsBunchLen/clight

56        CONTINUE
            ! --------------------------------------------------------------
            ! WRITING OUT THE ABOVE VARIABLE VALUES FOR DEBUGGING PURPOSES 
            ! WRITE(6,*) ampt,ampPt
            ! --------------------------------------------------------------
            ! r1 AND r2 ARE UNIFORM ON (-1,1)
            !---------------------------------------------------------------

            r1   = 2*ran3(iseed)-1
            r2   = 2*ran3(iseed)-1
            amp  = r1*r1+r2*r2    ! IS THE SQUARED VALUE OF THE DISTANCE FROM THE CENTER

            IF (amp >= 1) GO TO 56

            facc = sqrt(-2.*log(amp)/amp)
            tk   = ampt*r1*facc
            ptk  = ampPt*r2*facc
            p1   = nharm*tk*omega0
            p2   = nharm2*tk*omega0
            ham  = 0.5*ptk*ptk*tcoeff+(COS(p1)-1)*vrf/(nharm*v00*omega0) + (COS(p2)-1)*vrf2/(nharm2*v00*omega0)

            IF (ham >= hammax) GO TO 56      ! RESTART SAMPLING IF OUTSIDE BUCKET

            IF (ABS(tk) >= tauhat) GO TO 56  ! RESTART SAMPLING IF T IS TOO LARGE
            
            t(np)  = tk
            pt(np) = ptk

       CASE (3) ! GENERATE MATCHED 'PSEUDO-GAUSSIAN' PHASE SPACE
          ! ASSUMING DISTRIBUTION FUNCTION FOR HAMILTONIAN RHO(H)=1/H*EXP(-H/H), SO RHO(PT,T)~EXP(-C1*PT^2-C2*(1-COS(C3*T)). AT 2ND ORDER IN T THIS IS A BI-GAUSSIAN
          ! AVERAGE H CAN NOT BE CALCULATED ANALYTICALLY FOR THE EXACT HAMILTONIAN, ONLY FOR SMALL ANGLE APPROXIMATION.
          ! USE FIRST THIS VALUE TO CALCULATE BUNCH LENGTH, THEN USE NEWTON'S METHOD TO CALCULATE NEW H AND ITERATE UNTIL CONVERGENCE

          ! GET AVERAGE H IN SMALL OSC. APPROX. SEE DERIVATION IN CHECK_PHASE_SPACE_MATCHING.NB
          
          ham1sig0  = -nharm2*omega0*vrf2/(clight**2*v00)*rmsBunchLen**2
          
          CALL sampleLongMatched(np,t,pt,ham1sig0,sigs0,tauhat,ptmax,v00,omega0,hammax)
          
          ham1sig1  = ham1sig0*(rmsBunchLen/sigs0)**2
          
          CALL sampleLongMatched(np,t,pt,ham1sig1,sigs1,tauhat,ptmax,v00,omega0,hammax)

          ! USE NEWTON'S METHOD AND ITERATE UNTIL WE GET THE RIGHT BUNCH LENGTH
          ! CONSIDER SIGS=F(H) AS FUNCTION OF HAM1SIG. WE WANT TO FIND ZERO OF F(H)-SRMS
          ! ITERATE H1=H0 - (F(H0)-SRMS)/F'(H0)
          ! APPROXIMATE F'=(F(H1)-F(H0))/(H1-H0)

          DO WHILE (ABS(sigs1/rmsBunchLen-1.) >= bunchLenPrecis)
            WRITE(6,*) '       iterating: bunch length = ',sigs1,' m.'
            fPrime     = (sigs1-sigs0)/(ham1sig1-ham1sig0)
            ham1sig0   = ham1sig1
            sigs0      = sigs1
            ham1sig1   = ham1sig0-(sigs0-rmsBunchLen)/fPrime
            CALL sampleLongMatched(np,t,pt,ham1sig1,sigs1,tauhat,ptmax,v00,omega0,hammax)
            ! ---------------------------------------------------------------------------
            ! WRITING OUT SOME VARIABLES FOR DEBUGGING PURPOSES
            ! WRITE(*,*) ham1sig0, ham1sig1, sigs0, sigs1
            ! ----------------------------------------------------------------------------
          ENDDO

          WRITE(6,*) '       iterating: bunch length = ',sigs1,' m.'
      
      CASE (DEFAULT)
          WRITE(6,*) 'unknown method for long. coord.'
      
      END SELECT
      
      sigv = sigv + pt(np)**2

    WRITE(6,*) transvCoordMethod
    !------------------------------------------------------------------------------------------  
    ! GENERATING TRANSVERSE COORDINATE
    !------------------------------------------------------------------------------------------       
    SELECT CASE (transvCoordMethod)
      CASE (2)  ! SAMPLE DOUBLE GAUSSIAN
          ! ----------------------------------
          ! DETERMINE TRANSVERSE RMS AMPLITUDE
          ! ---------------------------------- 
          ampx = sqrt(betax*emix)
          ampy = sqrt(betay*emiy)
          
          WRITE(6,*)
          WRITE(6,*) 'generating transverse coordinates - double Gaussian'

          DO np = 1,nMacro
44          CONTINUE
            ! -------------------------------
            ! R1 AND R2 ARE UNIFORM ON (-1,1)
            ! -------------------------------
            r1  = 2*ran3(iseed)-1
            r2  = 2*ran3(iseed)-1
            amp = r1*r1+r2*r2

            IF ((amp >= 1).or.(amp <= 3.e-6)) GO TO 44

            facc = SQRT(-2.*LOG(amp)/amp)

            ! ---------------------------------------------
            ! INJECT WITH NO ANGLE ERROR
            !        X(NP) = XINJECT + AMPX*R1*FACC
            ! TRANSVERSE KICK WILL BE GIVEN ON TURN NTURNON
            ! ---------------------------------------------
            
            x(np) = ampx*r1*facc

            ! ---------------------------------------------
            ! PX HAS SAME AMPLITUDE AS X IE PX = BETA_L*X'
            ! ---------------------------------------------

            px(np) = ampx*r2*facc

            ! ----------------------------------------------------------------------------
            ! REJECT IF X > INITIAL CUT OFF (GIVEN IN SIGMAS WITH THE REFERENCE EMITTANCE)
            ! ----------------------------------------------------------------------------

            IF (SQRT(x(np)**2+px(np)**2) >= xcut) GO TO 44 ! XCUT = CUTOFFAMPL*SQRT(BETAX*REFEMXY) WHERE CUTOFFAMPL IS DEFINED IN THE INPUT FILE - COLLIMATION
            
45          CONTINUE
            ! -------------------------------
            ! R1 AND R2 ARE UNIFORM ON (-1,1)
            ! -------------------------------
            r1  = 2*ran3(iseed)-1
            r2  = 2*ran3(iseed)-1
            amp = r1*r1+r2*r2

            IF ((amp >= 1).or.(amp <= 3.e-6)) GO TO 45 ! OUTSIDE REITERATE

            facc    = SQRT(-2.*log(amp)/amp)
            y(np)   = ampy*r1*facc
            py(np)  = ampy*r2*facc

            ! ----------------------------------------------------------------------------
            ! REJECT IF X > INITIAL CUT OFF (GIVEN IN SIGMAS WITH THE REFERENCE EMITTANCE)
            ! ----------------------------------------------------------------------------

            IF (sqrt(y(np)**2+py(np)**2) >= ycut) GO TO 45 ! YCUT  = CUTOFFAMPL*SQRT(BETAY*REFEMXY) WHERE CUTOFFAMPL IS DEFINED IN THE INPUT FILE - COLLIMATION
            
            ! -------------------------------------------------------------------------------------------------------------------
            ! ANOTHER OPTION FOR COLLIMATION
            !
            ! REJECT IF SKEW AMPLITUDE > CUTOFF
            !          
            ! IF ((y(np)**2+py(np)**2)/(2*betay*refEmxy)+(x(np)**2+px(np)**2)/(2*betax*refEmxy) >= cutoffAmpl) GO TO 44
          ENDDO

      CASE (1)                     ! READ FROM FILE
          IF (readfile1Transv) THEN ! logical variable to keep track of beam 1 and beam 2
             open(unit=10,file=trim(filebeam1))
             readfile1Transv=.false.
             WRITE(6,*) 'reading transverse starting coordinates for beam1 from ', filebeam1
          ELSE
             open(unit=10,file=trim(filebeam2))
             WRITE(6,*) 'reading transverse starting coordinates for beam2 from ', filebeam2
          ENDIF

          DO np = 1,nMacro
             READ(10,*) y(np),py(np),x(np),px(np),dum1,dum2
          ENDDO
          CLOSE(10)
      CASE DEFAULT
          WRITE(6,*) 'Unknown method for transverse coordinates! Should be 1 or 2.'
          STOP
    END SELECT

    sigv = sqrt(sigv/np)

! ------------------------------------------------------------------------------------------
! DEBUGGING
! *******************************************************************************************
! WRITE(6,*)' rms of pt = ',sigv
! ------------------------------------------------------------------------------------------

! ------------------------------------------------------------------------------------------
! INITIAL NUMBER OF MACRO-PARTICLES FOR NORMALIZATION OF CURRENT ETC
! ------------------------------------------------------------------------------------------
    np    = np-1
    np0   = np

    CALL fitemittance(x,px,np,betax,cutoffAmpl,testemit)
    WRITE(6,*) 'emitx from fit = ',testemit
    CALL fitemittance(y,py,np,betay,cutoffAmpl,testemit)
    WRITE(6,*) 'emity from fit = ',testemit

! compare with emittance from standard deviation
       call emitStDev(x,px,np,betax,testemit)
       write(6,*) 'emitx from sigma = ',testemit
       call emitStDev(y,py,np,betay,testemit)
       write(6,*) 'emity from sigma = ',testemit

! calculate emittance from standard deviation backwards
       call calcEmitWCut(x,px,np,xcut,betax,testemit)
       write(6,*) 'emitx from calc. = ', testemit
       call calcEmitWCut(y,py,np,ycut,betay,testemit)
       write(6,*) 'emity from calc. = ', testemit
       write(6,*)

       return
     end subroutine addparticles

end module addparticles