! ---------------------------------------------------------------------------------------------------------------
! AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 1.0 : 
!   AUTHOR    : TOM MERTENS
!   DATE      : 16/06/2015
!   COPYRIGHT : CERN
!
!   DESCRIPTION : 
!       FUNCTION TO READ IN THE INPUT FILES CONTAINING THE SIMULATION INPUT ('collider_time_evolution.in')
!
!   REFERENCE   : Chao, Tigner: Handbook of Accelerator physics and engineering, (1998)
! ---------------------------------------------------------------------------------------------------------------
MODULE ModGetInput
	
	USE COMMONDATA
!  USE f90_unix_proc ! MODULE TO USE WITH NAGFOR COMPILER TO ACCESS THE CALL SYSTEM COMMAND

	IMPLICIT NONE

INTERFACE getinput
  MODULE PROCEDURE getinput
END INTERFACE

CONTAINS
	SUBROUTINE getinput(pnumber1,pnumber2,emix1,emiy1,emix2,emiy2,tauhat1,tauhat2,rmsDelta1,rmsBunchLen1,rmsDelta2,rmsBunchLen2)

  ! ------------------------------------------------------------------------------------------------------------------
  ! ******************* DESCRIPTION VARIABLES ************************************************************************
  ! ------------------------------------------------------------------------------------------------------------------
  !   
  ! *** INTEGER KIND I1B *********************************************************************************************
  ! IOstatus          : FLAG TO CHECK WHEN END OF FILE IS REACHED
  !
  ! *** INTEGER KIND I4B *********************************************************************************************
  ! k                 :
  ! nfft              :
  ! m                 :
  ! i                 :
  ! j                 :
  !
  ! *** REAL KIND SP *** OUTPUT **************************************************************************************
  ! pnumber1          : STARTING BUNCH POPULATION BUNCH BEAM 1
  ! pnumber2          : STARTING BUNCH POPULATION BUNCH BEAM 2
  ! emix1             : GEOMETRIC 1 SIGMA STARTING HORIZONTAL EMITTANCE BUNCH BEAM 1 (METERS)
  ! emiy1             : GEOMETRIC 1 SIGMA STARTING VERTICAL EMITTANCE BUNCH BEAM 1 (METERS)
  ! emix2             : GEOMETRIC 1 SIGMA STARTING HORIZONTAL EMITTANCE BUNCH BEAM 2 (METERS)
  ! emiy2             : GEOMETRIC 1 SIGMA STARTING VERTICAL EMITTANCE BUNCH BEAM 2 (METERS)
  ! tauhat1           : HALF BUCKET LENGTH (SECONDS) FOR BEAM 1 OF THE RF SYSTEM WITH SHORTEST WAVELENGTH, NOTE: FOR SINGLE RF, tauhat=thib/2
  ! tauhat2           : HALF BUCKET LENGTH (SECONDS) FOR BEAM 2 OF THE RF SYSTEM WITH SHORTEST WAVELENGTH, NOTE: FOR SINGLE RF, tauhat=thib/2
  ! rmsDelta1         : RMS OF DELTAP/P0 FOR BEAM 1 (ONLY USED WITH LONGCOORDMETHOD=2, SEE COMMON MODULE)
  ! rmsBunchLen1      : RMS OF BUNCH LENGTH FOR BEAM 1 (IN METRES, USED ONLY WITH LONGCOORDMETHOD=2 AND LONGCOORDMETHOD=3)
  ! rmsDelta2         : RMS OF DELTAP/P0 FOR BEAM 2 (ONLY USED WITH LONGCOORDMETHOD=2, SEE COMMON MODULE)
  ! rmsBunchLen2      : RMS OF BUNCH LENGTH FOR BEAM 2 (IN METRES, USED ONLY WITH LONGCOORDMETHOD=2 AND LONGCOORDMETHOD=3)
  !
  ! *** REAL KIND SP *************************************************************************************************
  ! eta               : SLIPFACTOR OF THE MACHINE
  ! tunesexact        :
  ! rcoeff            :
  ! tunes             :
  ! v1                :
  ! s                 :
  ! dgamma_hat        :
  ! rIon              :
  ! CalphaE3C         :
  ! I2                : RADIATION DAMPING INTEGRAL
  ! I4x               : RADIATION DAMPING INTEGRAL
  ! I4y               : RADIATION DAMPING INTEGRAL
  ! I3                : RADIATION DAMPING INTEGRAL
  ! sigEoE0           :
  ! Cq                :
  ! alfax             : RADIATION DAMPING PARAMETER - METHOD APPROX
  ! alfay             : RADIATION DAMPING PARAMETER - METHOD APPROX
  ! gamax             : RADIATION DAMPING PARAMETER - METHOD LATTIC
  ! gamay             : RADIATION DAMPING PARAMETER - METHOD LATTIC
  ! Dx                : RADIATION DAMPING PARAMETER - METHOD APPROX
  ! Dy                : RADIATION DAMPING PARAMETER - METHOD APPROX
  ! DxP               : RADIATION DAMPING PARAMETER - METHOD APPROX
  ! DyP               : RADIATION DAMPING PARAMETER - METHOD APPROX
  ! Hx                : RADIATION DAMPING PARAMETER - METHOD APPROX
  ! Hy                : RADIATION DAMPING PARAMETER - METHOD APPROX
  ! I5x               : RADIATION DAMPING INTEGRAL
  ! I5y               : RADIATION DAMPING INTEGRAL
  ! Jx                :
  ! Jy                :
  ! ki                :
  ! rhoi              :
  !
  ! *** CHARACTER ARRAY LEN = 4 **************************************************************************************
  ! col1              :
  ! col2              :
  ! col3              :
  ! col4              :
  ! col5              :
  ! col6              :
  ! col7              :
  ! col8              :
  ! col9              :
  ! col10             :
  ! col11             :
  ! col12             :
  ! col13             :
  ! col14             :
  ! col15             :
  !
  ! *** CHARACTER ARRAY LEN = 150 ************************************************************************************
  ! twissFile         : TFS TYPE TWISS OUTPUT FILE GENERATED BY MADX USING THE DESIRED OPTICS
  ! ibsIntFile        : FILE CONTAINING PRECALCULATED IBS INTEGRALS, ALLOWS TO INTERPOLATE
  ! gTabFile          : TABLE CONTAINING G FUNCTION USED BY BANE TO CALCULATE IBS LIFE TIMES
  ! 
  ! ------------------------------------------------------------------------------------------------------------------
  ! ******************* DESCRIPTION VARIABLES **** END ***************************************************************
  ! ------------------------------------------------------------------------------------------------------------------
    IMPLICIT NONE
		REAL(SP) ,    INTENT(OUT)     :: pnumber1,pnumber2,emix1,emiy1,emix2,emiy2,tauhat1,tauhat2, &
							                       rmsDelta1,rmsBunchLen1,rmsDelta2,rmsBunchLen2

		REAL(SP)                      :: eta,tunesexact,rcoeff,tunes, &
           		                       v1,s,dgamma_hat,rIon,CalphaE3C,I2,I4x,I4y,I3, &
           		                       sigEoE0,Cq,alfax,alfay,gamax,gamay,Dx,Dy,DxP,DyP, &
           		                       Hx,Hy,I5x,I5y,Jx,Jy,ki,rhoi

    INTEGER(I1B)                  :: IOstatus              
	  INTEGER(I4B)                  :: k,nfft,m,i,j
    
    CHARACTER(len=4)              :: col1,col2,col3,col4,col5,col6,col7,col8,col9, &
          				                   col10,col11,col12,col13,col14,col15
   		
		CHARACTER(len=150)            :: twissFile,ibsIntFile,gTabFile

! -------------------------------------------------------------------------------------------------------------
! ******************* INITIALIZATION **************************************************************************
! -------------------------------------------------------------------------------------------------------------

! -----------------------------------------------
! INITIALIZATION OF REAL OUTPUT VARIABLES
! -----------------------------------------------   
  pnumber1      = 0.0
  pnumber2      = 0.0
  emix1         = 0.0
  emiy1         = 0.0
  emix2         = 0.0
  emiy2         = 0.0
  tauhat1       = 0.0
  tauhat2       = 0.0
  rmsDelta1     = 0.0
  rmsBunchLen1  = 0.0
  rmsDelta2     = 0.0
  rmsBunchLen2  = 0.0

! -----------------------------------------------
! INITIALIZATION OF REAL PRIVATE VARIABLES
! ----------------------------------------------- 
  eta             = 0.0
  tunesexact      = 0.0
  rcoeff          = 0.0
  tunes           = 0.0
  v1              = 0.0
  s               = 0.0
  dgamma_hat      = 0.0
  rIon            = 0.0
  CalphaE3C       = 0.0
  I2              = 0.0
  I4x             = 0.0
  I4y             = 0.0
  I3              = 0.0
  sigEoE0         = 0.0
  Cq              = 0.0
  alfax           = 0.0
  alfay           = 0.0
  gamax           = 0.0
  gamay           = 0.0
  Dx              = 0.0
  Dy              = 0.0
  DxP             = 0.0
  DyP             = 0.0
  Hx              = 0.0
  Hy              = 0.0
  I5x             = 0.0
  I5y             = 0.0
  Jx              = 0.0
  Jy              = 0.0
  ki              = 0.0
  rhoi            = 0.0
  
! -------------------------------------------------
! INITIALIZATION OF INTEGER (I1B) PRIVATE VARIABLES
! -------------------------------------------------
  IOstatus=0             ! FLAG TO CHECK WHEN END OF FILE IS REACHED

! -------------------------------------------------
! INITIALIZATION OF INTEGER (I4B) PRIVATE VARIABLES
! -------------------------------------------------
  k       = 0
  nfft    = 0
  m       = 0
  i       = 0
  j       = 0

! -----------------------------------------------------------
! INITIALIZATION OF CHARACTER ARRAY (LEN=4) PRIVATE VARIABLES
! -----------------------------------------------------------
  col1  = 'XXXX'
  col2  = 'XXXX'
  col3  = 'XXXX'
  col4  = 'XXXX'
  col5  = 'XXXX'
  col6  = 'XXXX'
  col7  = 'XXXX'
  col8  = 'XXXX'
  col9  = 'XXXX'
  col10 = 'XXXX'
  col11 = 'XXXX'
  col12 = 'XXXX'
  col13 = 'XXXX'
  col14 = 'XXXX'
  col15 = 'XXXX'
 
! -------------------------------------------------------------
! INITIALIZATION OF CHARACTER ARRAY (LEN=150) PRIVATE VARIABLES
! -------------------------------------------------------------
!  twissFile   = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX &
!                 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX &
!                 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
!
! ibsIntFile  = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX &
!                 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX &
!                 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
!
!  gTabFile    = 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX &
!                 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX &
!                 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'

! -------------------------------------------------------------------------------------------------------------
! ******************* INITIALIZATION **** END *****************************************************************
! -------------------------------------------------------------------------------------------------------------

! -------------------------------------------
! OPENING INPUT FILE WITH SIMULATION SETTINGS
! -------------------------------------------
  
  OPEN(unit=10,file='collider_time_evolution.in',status='unknown')

  READ(10,*) col1
  WRITE(6,*) col1

! ------------------------------------------------------------
! READ THE SWITCHES TO TURN ON THE VARIOUS AVAILABLE PROCESSES
! ------------------------------------------------------------

  READ(10,*) RFswitch       		!(SET TO 1 TO ACTIVATE SYNCHROTRON MOTION)
  READ(10,*) betatronSwitch 		!(SET TO 1 TO ACTIVATE BETATRON MOTION)
  READ(10,*) raddampSwitch  		!(SET TO 1 TO ACTIVATE RADIATION DAMPING AND QUANTUM EXCITATION)
  READ(10,*) IBSswitch      		!(SET TO 1 TO ACTIVATE IBS)
  READ(10,*) collimationSwitch 	!(SET TO 1 TO ACTIVATE LOSSES ON APERTURE CUTS)
  READ(10,*) blowupSwitch 		  !(SET TO 1 TO ACTIVATE ARTIFICIAL BLOWUP - ADT LOSS MAPS ETC)
  READ(10,*) collisionSwitch 		!(SET TO 1 TO ACTIVATE LUMINOSITY)
  READ(10,*) levellingSwitch 		!(SET TO 1 TO ACTIVATE BETA LEVELLING OF LUMINOSITY IN ONE OR MORE IPS)
  READ(10,*) col1               ! TO READ SEPARATION LINE
   
! --------------------------------------------------------------------------
! **************************************************************************
! SOME OF THE VARIABLES BELOW ARE DEFINED AND EXPLAINED IN THE COMMON MODULE
! **************************************************************************
! --------------------------------------------------------------------------

! ------------------------------------------------------------
! READ THE SWITCHES TO SELECT OUTPUT FILES TO BE GENERATED
! ------------------------------------------------------------
  READ(10,*) writeLuminositySwitch    ! SET TO 1 TO WRITE OUT THE LUMINOSITY
  READ(10,*) writeIntensitySwitch     ! SET TO 1 TO WRITE OUT INTENSITIES
  READ(10,*) writeEmittanceSwitch     ! SET TO 1 TO WRITE OUT EMITTANCES
  READ(10,*) writeIBSSwitch           ! SET TO 1 TO WRITE OUT IBS
  READ(10,*) writeCollImpactsxbetatron ! SET TO 1 TO WRITE OUT COLLIMATION IMPACT ON X BETATRON
  READ(10,*) writeCollImpactsybetatron ! SET TO 1 TO WRITE OUT COLLIMATION IMPACT ON Y BETATRON
  READ(10,*) writeCollImpactsxmomentum ! SET TO 1 TO WRITE OUT COLLIMATION IMPATC ON X MOMENTUM
  READ(10,*) WRITEAllCoordSwitch      ! SET TO 1 TO ACTIVATE A DUMP OF ALL PARTICLE COORINATES ON EVERY TURN WHERE OUTPUT IS WRITTEN 
  READ(10,*) WRITEMountSwitch         ! SET TO 1 TO ACTIVATE MOUNTAIN RANGE OUTPUT OF DISTRIBUTION SHAPE IN TRANSV AND LONG PLANES
  READ(10,*) WRITEEqTimeSwitch        ! SET TO 1 TO WRITE REAL TIME TO SEPARATE FILE
  READ(10,*) WRITEScreenSummarySwitch ! SET TO 1 TO WRITE A SUMMARY OF THE SIMULATION SETTINGS TO THE SCREEN AT START OF SIMULATION
  READ(10,*) col1                     ! TO READ SEPARATION LINE

! ------------------------------------------------------------
! READ SIMULATION INPUT SETTINGS
! ------------------------------------------------------------
 	READ(10,*) nturns,nMacro,timeRatio,nWRITE,iseed
  WRITE(6,*) '    nturns=',nturns,', # macro part.=',nMacro

! ------------------------------------------------------------
! READ IN MACHINE PARAMETERS
! ------------------------------------------------------------
 	READ(10,*) gammat,circ,gamma0
  READ(10,*) vrf,nharm,vrf2,nharm2
  READ(10,*) tunex,tuney,chromx,chromy,dqmin,k2L,k2Lskew
! ------------------------------------------------------------
! SMOOTH APPROXIMATION OF LATTICE FUNCTIONS FOR LATER USE
! ------------------------------------------------------------
  IF (tunex >= 1.0E-9) THEN ! CHECK IF TUNEX IS POSITIVE
    betax=circ/(2*pi*tunex)
  ELSE
    WRITE(6,*) 'tunex = 0, division by zero'
    STOP
  END IF

  IF (tuney >= 1.0E-9) THEN ! CHECK IF TUNEY IS POSITIVE
    betay=circ/(2*pi*tuney)
  ELSE
    WRITE(6,*) 'tuney = 0, division by zero'
    STOP
  END IF

! ------------------------------------------------------------
! WRITE SMOOTH APPROXIMATION TO SCREEN
! ------------------------------------------------------------
  WRITE(6,*)
  WRITE(6,*) 'reference beta functions x,y (m) = ',betax,betay
  WRITE(6,*)

! ------------------------------------------------------------
! CONTINUE TO READ SIMULATION INPUT SETTINGS
! ------------------------------------------------------------
 	READ(10,*) aatom,qatom,thib     
  
! --------------------------------------------------------------------------------------------
! OPTICS TWISSFILE - USED FOR IBS WITH IBSMETHOD='PIWLATTICE', 'MODPIWLAT', OR 'NAGAITSEV' AND 
! FOR RADIATION DAMPING WITH RADMETHOD='LATTIC'FILES WITH LONG. 
! STARTING CONDITIONS, USED ONLY IF LONGCOORDMETHOD=1
! ---------------------------------------------------------------------------------------------
 89   FORMAT(A200)            ! SET FORMAT FOR READING LINES
      READ(10,89) twissFile   ! READING IN OF THE TWISS FILE GENERATED BY MAD
      READ(10,*) col1         ! TO READ SEPARATION LINE

! ------------------------------------------------------------
! CONTINUE TO READ SIMULATION INPUT SETTINGS
! ------------------------------------------------------------
 	READ(10,*) emix1,emiy1,pnumber1 ! INITIAL TRANSVERSE EMITTANCES AND PARTICLES IN BUNCH FOR BEAM 1
  READ(10,*) emix2,emiy2,pnumber2 ! INITIAL TRANSVERSE EMITTANCES AND PARTICLES IN BUNCH FOR BEAM 2
 	READ(10,*) longCoordMethod, transvCoordMethod ! READING METHODS TO GENERATE LONGITUDINAL AND TRANSVERSE PARTICLE DISTRIBUTIONS
  READ(10,*) rmsBunchLen1,rmsBunchLen2
  READ(10,*) rmsDelta1,rmsDelta2
 	READ(10,*) tauhat1,tauhat2,bunchLenPrecis,power,alint
 	READ(10,89) filebeam1 ! WHEN USING EXTERNALLY GENERATED INITIAL CONDITIONS AVAILABLE IN EXTERNAL FILE DEFINED IN 'collider_time_evolution.in'
 	READ(10,89) filebeam2 ! WHEN USING EXTERNALLY GENERATED INITIAL CONDITIONS AVAILABLE IN EXTERNAL FILE DEFINED IN 'collider_time_evolution.in'
 	READ(10,*) col1       ! TO READ SEPARATION LINE
  READ(10,*) radMethod  ! RADIATION DAMPING
 	READ(10,*) tradlong,tradperp,siglong,sigperp ! FOR MANUAL RADIATION DAMPING
 	READ(10,*) rho0       ! READING THE DIPOLE BENDING RADIUS
  READ(10,*) col1       ! TO READ SEPARATION LINE

  READ(10,*) IBSMETHOD    !	READING THE METHOD TO CALCULATE EMITTANCE CHANGES DUE TO IBS
  READ(10,*) coupleIBS,coulombLog ! COULOMBLOG IS A VALUE USED IN NAGAITSEV IBS ROUTINE
  READ(10,*) fracibstot, nBins 
  READ(10,89) ibsIntFile  ! FILE WITH INTERPOLATED VALUES FOR IBS, USED WITH IBSMETHOD='INTERPOLAT'
  READ(10,89) gTabFile    ! FILE WITH TABULATED VALUES OF BANE'S G-FUNCTION. USED ONLY WITH IBSMETHOD='BANEAPPROX'  
  READ(10,*) col1         ! TO READ SEPARATION LINE

! ------------------------------------------------------------
! CONTINUE TO READ SIMULATION INPUT SETTINGS
! VARIABLES RELATED TO THE COLLIMATION PROCESS
! ------------------------------------------------------------
  READ(10,*) refEmxy,cutoffAmpl,collimAvgSwitch,emitMethod
  WRITE(6,*) 'Emittance generating method = ',emitMethod
  READ(10,*) nSigCutBeta,nSigCutMom,betaxMom,dispxMom
  READ(10,*) col1 ! TO READ SEPARATION LINE


  xcut  = cutoffAmpl*sqrt(betax*refEmxy) ! CUTOFF IN METRES OF THE TRANSVERSE DISTRIBUTIONS X=SQRT(BETA REFEMXY) FOR MATCHED DISTRIBUTION
  ycut  = cutoffAmpl*sqrt(betay*refEmxy) ! CUTOFF IN METRES OF THE TRANSVERSE DISTRIBUTIONS X=SQRT(BETA REFEMXY) FOR MATCHED DISTRIBUTION

! ------------------------------------------------------------
! CONTINUE TO READ SIMULATION INPUT SETTINGS
! VARIABLES RELATED TO THE BLOWUP PROCESS
! ------------------------------------------------------------
  READ(10,*) pxKickFac,pyKickFac,blowupMethod
  READ(10,*) col1 ! TO READ SEPARATION LINE

! ------------------------------------------------------------
! CONTINUE TO READ SIMULATION INPUT SETTINGS
! VARIABLES RELATED TO THE COLLISION PROCESS
! ------------------------------------------------------------

  READ(10,*) collRoutine
  READ(10,*) nIPs,sigI
  READ(10,*) nbunches,longIntBins,angleSwitch

! --------------------------------------------------------------------------
! ADD LUMIMAX AND LEVIPSWITCH TO INPUT PARAMETERS FOR BETA* LEVELLING PER IP
! --------------------------------------------------------------------------
      DO m=1,nIPs  ! LOOP OVER ALL IPS, READ BETA*, CROSSING ANGLE/2 AND MULTIPLICITY
         READ(10,*) betaS(m),theta(m),IPmult(m),lumimax(m),levIPSwitch(m)
         WRITE(6,*) lumimax(m)
         betaSMin(m) = betaS(m)
      ENDDO

!---------------------------------
! CALCULATING DYNAMICAL PARAMETERS
!---------------------------------
  beta        =   sqrt(1-1/gamma0**2)         ! RELATIVISTIC BETA
  vrev        =   clight*beta                 ! RELATIVISTIC SPEED
  trev        =   circ/vrev                   ! REVOLUTION TIME FOR THE PARTICLE
  ! NEED TO CHECK BECAUSE BELOW IS FOR FIRST HARMONIC RF SYSTEM AND NOT FOR SECOND ONE WHICH IS THE LHC'S
  trf         =   trev/nharm                  ! tref \approx 2.5E-9 REVOLUTION TIME DIVIDED BY HARMONIC NUMBER LHC = 35640 FOR THE 
  frf         =   1/trf                       ! frf = 400 MHz RF FREQUENCY      
  omegarf     =   2*pi/trf                    ! RF ANGELUR FREQUENCY
  eta         =   1/gammat**2 - 1/gamma0**2   ! SLIPFACTOR SEE Lee (2.162)

!-------------------------------------
! CALCULATE EQUIVALENT SIMULATION TIME
!-------------------------------------
  eqtime      =   nturns * trev * timeRatio/3600.     ! RETURNS THE REAL TIME ELAPSED IN HOURS

  WRITE(6,*)'    equivalent time (hours) = ', eqtime
  WRITE(6,*)

  IF (WRITEEqTimeSwitch == 1) THEN
    OPEN(unit=34,file='eqtime.out')
    WRITE(34,*)' equivalent time (hours) = ', eqtime
  END IF

!--------------------------------------------
! WRITE SUMMARY ON SCREEN AT SIMULATION START
!--------------------------------------------
Wsumm:    IF (WRITEScreenSummarySwitch==1) THEN
            IF (RFswitch == 1) THEN
              WRITE(6,*) '    RF motion is:         ON'
            ELSE
              WRITE(6,*) '    RF motion is:         OFF'
            END IF
          
            WRITE(6,*) ! WRITING AN EMPY LINE TO THE SCREEN
          
            IF (betatronSwitch == 1) THEN
               WRITE(6,*) '    betatron motion is:   ON'
               WRITE(6,*) '        -->reference beta functions:'
               WRITE(6,*) '          circ/(2*pi*tunex) = ',circ/(2*pi*tunex),',  circ/(2*pi*tuney) = ',circ/(2*pi*tuney)
            ELSE
               WRITE(6,*) '    betatron motion is:   OFF'
            END IF
          
            WRITE(6,*) ! WRITING AN EMPY LINE TO THE SCREEN
          
            IF (collimationSwitch == 1) THEN
               WRITE(6,*) '    Collimation is:        ON'
               WRITE(6,*) '       --> number of aperture limitations:',nAperCuts
            ELSE
               WRITE(6,*) '    Collimation is:       OFF'
            END IF
          
            WRITE(6,*) ! WRITING AN EMPY LINE TO THE SCREEN
          
            IF (blowupSwitch == 1) THEN
               WRITE(6,*) '    Blowup is:             ON'
            ELSE
               WRITE(6,*) '    Blowup is:            OFF'
            END IF
          
            WRITE(6,*) ! WRITING AN EMPY LINE TO THE SCREEN
          
            IF (collisionSwitch == 1) THEN
               WRITE(6,*) '    Collisions are:       ON'
               WRITE(6,*) '       --> using collision routine:',collRoutine
          
               DO m=1,nIPs
                  WRITE(6,*)'       -->',IPmult(m),' IPs with beta*=',betaS(m), 'm, phi/2=',theta(m)
               END DO
            ELSE
               WRITE(6,*) '    Collisions are:       OFF'
            END IF
          
            WRITE(6,*) ! WRITING AN EMPY LINE TO THE SCREEN

! --------------
! BETA LEVELLING
! -------------- 
  
            IF (levellingSwitch.eq.1) THEN  
              WRITE(6,*) 'Levelling is :        ON'
              DO m=1,nIPs
                WRITE(6,*) '       ----> IP ',m, ' is levelled: ', levIPSwitch(m), ' levelling value for luminosity: ',lumimax(m)
              END DO
            ELSE
              WRITE(6,*) 'Levelling is :        OFF'
            END IF

            WRITE(6,*) 

            WRITE(6,*)'       SUMMARY OF SIMULATION PARAMETERS'
            WRITE(6,*)'eta = ',eta

          END IF Wsumm

! ----------------------------------------------------
! CLOSING THE INPUT FILE 'collider_time_evolution.in'
! ----------------------------------------------------
  CLOSE(10)
  WRITE(6,*)  ! WRITING AN EMPY LINE TO THE SCREEN

! ----------------------------------------------------
! OBSOLETE CODE ???
! ----------------------------------------------------
! makesure nresamp is a power of 2 ??? where is this nresamp defined the first time ???
!     nfft=2
!      do k=1,1000
!         nfft = 2*nfft
!         if(nfft.ge.nresamp)go to 22
!      enddo
! 22   continue
!      WRITE(6,*)'nresamp,nfft',nresamp,nfft
!      nresamp=nfft


!      
! longitudinal eom, t(k) = arrival time with synchronous point at pi rf radians
!                  pt(k) = gamma-gamma0
! d pt(k)/dn = e*qatom*vrf*omegarf*(t(k)-trf/2)/(aatom*pmass*clight**2)
! d  t(k)/dn = trev*eta*pt(k)/(beta0**2*gamma0)
! coefficients  d pt(k)/dn = pcoeff*(t(k)-trf/2) , d  t(k)/dn = tcoeff*pt(k)
      IF (vrf.eq.0) THEN
         v1=vrf2
      ELSE
         v1=vrf
      END IF
      pcoeff = qatom*v1*omegarf/(938.e6*aatom)
      tcoeff = trev*eta/(beta**2*gamma0)
      rcoeff = pcoeff/tcoeff
      if(rcoeff.gt.0)THEN
         WRITE(6,*)'wrong sign of voltage for sign of eta'
         STOP
      END IF
!      WRITE(6,*)' tuney,mode',tuney,mode
! get the frequency of the slowest betatron sideband
!      WRITE(6,*)'lowest sideband (MHz) = ',1.e-6*(tuney-mode)/trev
! given gap volts get (gamma-gamma0) corresponding to amplitude in arrival time variation

!---------
!      dgamma_hat = tauhat1*sqrt(-rcoeff)
!-----------
!tauhat1=half bucket length (seconds) for beam 1 of the RF system with shortest wavelength
!      WRITE(6,*)' amp of de/e = ',dgamma_hat/gamma0
!      WRITE(6,*)'max change in tau per turn = ',dgamma_hat*tcoeff
! synchrotron tune
! ---------
!      tunes = sqrt(-pcoeff*tcoeff)/(2*pi)
!-----------
!      WRITE(6,*) tunes,pcoeff,vrf,omegarf,qatom
!      stop
!------------
!      tunesexact = acos(1+pcoeff*tcoeff/2)/(2*pi)
!-------------
!      WRITE(6,*)' synchrotron frequency = ',tunes/trev
!      WRITE(6,*)' revolution period = ',trev

!      WRITE(6,*)' exact synchrotron tune = ',tunesexact
! numerical parms
!      WRITE(6,*)'tbin/sigtsmooth=', 2.506628*(thib-tlob)/
!     : (taupart*nresamp)
      CALL SYSTEM('rm intensity.out')
      CALL SYSTEM('rm luminosity.out')
      CALL SYSTEM('rm ibs.out')
      CALL SYSTEM('rm mountain*.out')
      CALL SYSTEM('rm emittance.out')
      CALL SYSTEM('rm eqtime.out')
      CALL SYSTEM('rm collimator_impacts_x_betatron.out')
      CALL SYSTEM('rm collimator_impacts_x_momentum.out')
      CALL SYSTEM('rm collimator_impacts_y_betatron.out')



!      open(unit=32,file='tran.xvz',status='unknown')
!      open(unit=33,file='csmon.out',status='unknown')
!      open(unit=52,file='nLostBurnoff.out',status='unknown')

!      open(unit=34,file='tran.cs',status='unknown')
!      open(unit=44,file='tran.sig',status='unknown')

!      call getimped !! used for cooling - uncomment later
! get fiducial normalizations
!      WRITE(6,*)' dpwake for 1 meter offset with step wake = 1'
!      WRITE(6,*)(circ/(2*pi*tuney))*1.6e-19*qatom*pnumber1/
!     :          (beta*gamma0*938.e6*aatom/qatom)

! normalization for longitudinal profiles
 !     nWRITE0 = nWRITE
! equivalent time in machine of simulation, rescaling with number of sim turns per real turns

! --------------------------------------------------
! RADIATION DAMPING
! --------------------------------------------------
! READ MAD TWISS FILE WITH LATTICE (FOR SOME METHODS OF IBS AND RADIATION DAMPING ONLY). ASSUME FIXED FORMAT OF TFS TABLE FOR SIMPLICITY
lat:        IF (((raddampSwitch == 1) .and. (radMethod /= "manual")) .or. &
                ((ibsSwitch == 1) .and. ((ibsMethod == 'piwLattice') .or. &
                (ibsMethod == 'modPiwLatt') .or. (ibsMethod == 'baneApprox') .or. (ibsMethod == 'nagaitsev')))) THEN
              
                open(unit=99,file=trim(twissFile)) ! OPENING TWISS FILE FOR READING OPTICAL LATTICE
                col1="||"
              
              ! READ UNTIL HEADER OF TWISS FILE

              DO i=1,45
                 READ(99,*) col1
              END DO

              WRITE(6,*) '    Reading twiss file ',trim(twissFile)
              WRITE(6,*)
              
              READ(99,*) col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15 ! READ COLUMN HEADERS

              ! CHECK THAT WE HAVE THE RIGHT COLUMNS IN TWISS FILE

              IF ((col2 == "NAME").and.(col3 == "S   ").and.(col4 == "L   ").and.(col5 == "BETX").and.(col6 == "BETY")&
                   .and.(col7 == "ALFX").and.(col8 == "ALFY").and.(col9 == "DX  ").and.(col10 == "DPX ")&
                   .and.(col11 == "DY  ").and.(col12 == "DPY ").and.(col13 == "ANGL").and.(col14 == "K1L ")&
                   .and.(col15 == "K1S ")) THEN
              ELSE
                 WRITE(6,*) "Error in format of twiss file ", trim(twissFile)
                 WRITE(6,*) "Columns must be NAME|S   |L   |BETX|BETY|ALFX|ALFY|DX  |DPX |DY  |DPY |ANGLE|K1L|K1S"
                 STOP
              END IF

              READ(99,*) col1        ! READ ONE MORE LINE TO GET TO THE ACTUAL TABLE

              ! READ IN LATTICE

              nElem     = 1                ! COUNT THE NUMBER OF ELEMENTS IN LATTICE

              DO WHILE (IOstatus == 0)     ! READ ELEMENTS UNTIL END OF FILE
                 READ(99,*,iostat=IOstatus) &
                      col1,S,leng(nElem),betx(nElem),bety(nElem),alfx(nElem),alfy(nElem),dispx(nElem),dispxP(nElem),dispy(nElem),&
                      dispyP(nElem),bendAng(nElem),K1L(nElem),K1S(nElem)
                  IF (leng(nElem) /=0 ) nElem = nElem+1
              END DO

              nElem = nElem-1
              CLOSE(99)   ! CLOSING THE TWISS FILE 
           END IF lat

ibss: IF (ibsSwitch == 1) THEN
         WRITE(6,*) '    IBS is:               ON'
         WRITE(6,*) '       -->using method: ',ibsMethod

ibscase:  SELECT CASE (ibsMethod)

          CASE ('baneApprox')
            ! --------------------------------------------- 
            ! READ IN TABULATED VALUES OF BANE'S G-FUNCTION
            ! ---------------------------------------------
            OPEN(unit=99,file=gTabFile) ! OPENING THE G-FUNCTION PRE-EVALUATED TABLE FILE
            
            DO i=1,1501
               READ(99,*) gTab(i,1),gTab(i,2)
            ENDDO
            
            CLOSE(99) ! CLOSE THE G-FUNCTION TABLE FILE

          CASE ('interpolat')
            ! --------------------------------------------- 
            ! READ IN TABULATED VALUES OF IBS LIFETIMES
            ! ---------------------------------------------
            OPEN(unit=10,file=trim(ibsIntFile))
            WRITE(6,*) '       -->reading tabulated IBS values from ',trim(ibsIntFile)
            READ(10,*) xMin,xMax,xBins
            READ(10,*) yMin,yMax,yBins
            READ(10,*) zMin,zMax,zBins
            READ(10,*)  Ap,Ax,Ay
            DO i=1,xBins
               DO j=1,yBins
                  DO k=1,zBins
                     READ(10,*) Ap(i,j,k),Ax(i,j,k),Ay(i,j,k)
                  END DO
               END DO
            END DO
            CLOSE(10) ! CLOSE IBS INTERPOLAT FILE
         END SELECT ibscase

      ELSE

         WRITE(6,*) '    IBS is:               OFF'
      
      END IF ibss

      WRITE(6,*)

      IF (raddampSwitch == 1) THEN
         WRITE(6,*) '    Radiation damping is: ON'
         WRITE(6,*) '       -->using method: ',radMethod
      ELSE
         WRITE(6,*) '    Radiation damping is: OFF'
      END IF


manrads: IF ((raddampSwitch == 1).and.(radMethod /= "manual")) THEN

      ! ------------------------------------------------------------------------------------------
      ! CALCUALTE RADIATION DAMPING TIMES AND EQUILIBRIUM EMITTANCES IF NOT GIVEN MANUALLY
      ! Reference: Chao, Tigner: Handbook of Accelerator physics and engineering, (1998) page 186
      ! ------------------------------------------------------------------------------------------

radmcase: SELECT CASE (radMethod)
          CASE ('approx') 
          ! ---------------------------------------------------------------------------------------------------------------------------------
          ! APPROX: SMOOTH LATTICE, ASSUME I4_J<<I2 FOR ALL J. SEE HANDBOOK P186 AND PUT I4_J=0, I2=2*PI/RHO0 AND C_ALPHA=RION/(3*C^5*MION^3)
          ! ADDITIONAL DIV BY 2 COMPARED TO HANDBOOK TO GET DAMPING TIME FOR EMITTANCE INSTEAD OF AMPLITUDE. SEE  PAPER DERIVATION.
          ! ----------------------------------------------------------------------------------------------------------------------------------
            gamax = (1.+alfax**2)/betax
            gamay = (1.+alfay**2)/betay
            Dx    = circ/(2*pi*gammat**2)
            Dy    = 0.
            DxP   = 0.1             ! SHOULD FIND AN APPROXIMATION FORMULA. HOWEVER NOT VERY IMPORTANT
            DyP   = 0.
            Hx    = (betax*DxP+2*alfax*Dx*DxP+gamax*Dx)
            Hy    = (betay*DyP+2*alfay*Dy*DyP+gamay*Dy)
          ! --------------------------------------------------
          ! DEFINE SMOOTH APPROXIMATION OF RADIATION INTEGRALS
          ! --------------------------------------------------
            I2    = 2*pi/rho0
            I3    = 2*pi/rho0**2
            I4x   = 0.0
            I4y   = 0.0
            I5x   = Hx*2*pi/rho0**2
            I5y   = Hy*2*pi/rho0**2

            CASE ('lattic') 
          ! -----------------------------------------------------------
          ! CALCULATE RADIATION INTEGRALS OVER LATTICE USING TWISS FILE
          ! -----------------------------------------------------------
            DO i=1,nElem
                IF (bendAng(i) /= 0.) THEN
                  rhoi  =   leng(i)/bendAng(i)
                  ki    =   K1L(i)/leng(i)
                  I2    =   I2+1./rhoi**2 *leng(i)
                  I3    =   I3+1./rhoi**3 *leng(i)
          ! ---------------------------------------------------------------------------------------------------------------------------------------------
          ! I4x=I4x+(dispx(i)/rhoi**3*(1.+2.*rhoi**2*ki)+2.*dispx(i)*ki/rhoi)*leng(i) corrected to equations in accelerator handbook second edition p 220
          ! I4y=I4y+(dispy(i)/rhoi**3*(1.-2.*rhoi**2*ki)-2.*dispy(i)*ki/rhoi)*leng(i)
          ! ---------------------------------------------------------------------------------------------------------------------------------------------
                  I4x   =   I4x+(dispy(i)/rhoi**3*(1.+2.*rhoi**2*ki)+2.*dispy(i)*K1S(i)/rhoi)*leng(i)
                  I4y   =   0
                  gamax =   (1.+alfx(i)**2)/betx(i)
                  gamay =   (1.+alfy(i)**2)/bety(i)
                  Hx    =   betx(i)*dispxP(i)**2 + 2*alfx(i)*dispx(i)*dispxP(i)+gamax*dispx(i)**2
                  Hy    =   bety(i)*dispyP(i)**2 + 2*alfy(i)*dispy(i)*dispyP(i)+gamay*dispy(i)**2
                  I5x   =   I5x+Hx*2*pi/rho0**2*leng(i)
                  I5y   =   I5y+Hy*2*pi/rho0**2*leng(i)
                END IF

            END DO

            WRITE(6,*) "          calculating radiation integrals: "
            WRITE(6,*) "          I2=",I2,"I3=",I3,"I4x=",I4x,"I4y=",I4y,"I5x=",I5x,"I5y=",I5y

         END SELECT radmcase

         rIon       =   (qatom**2/aatom)*1.54e-18
         CalphaE3C  =   rIon*gamma0**3/(3*trev)  ! C_alpha*E^3/C in handbook formulas p 186

         tradperp   =   1./(CalphaE3C*I2*(1.-I4x/I2))/2. ! eq 9, but div by 2 to get rise time for emittance and not amplitude
         tradlong   =   1./(CalphaE3C*I2*(2.+(I4x+I4y)/I2))/2.

         Cq         =   55./(32.*sqrt(3.0))*1.0546e-34*clight / (938.e6*aatom*1.60218e-19)
         sigEoE0    =   Cq*gamma0**2*I3/(2*I2+I4x+I4y)  ! eq 18
         siglong    =   gamma0*sigEoE0 ! sigma of pt
         Jx         =   1.-I4x/I2
         Jy         =   1.-I4y/I2
         sigperp    =   Cq*gamma0**2*I5x/(Jx*I2)


         WRITE(6,*) '          calculating radiation damping times and equilibrium beam sizes:'
         WRITE(6,*) '          Tx=',tradperp/3600,'h,  Ts=',tradlong/3600,'h,  sigPt=',siglong,'ex=',sigperp,'m'
      END IF manrads

      WRITE(6,*)
      WRITE(6,*)'... collider_time_evolution.in read'
      WRITE(6,*)

! OPEN OUTPUT FILES, WRITE HEADERS
      
      ! ---------------------------------------------------------------
      ! IF WRITELUMINOSITY SWITCH IS ON OPEN THE LUMINOSITY OUTPUT FILE
      ! ---------------------------------------------------------------

      IF (writeLuminositySwitch == 1) THEN 
        OPEN(unit=80,file='luminosity.out')
        WRITE(80,*) 'sim.turn     t(hours)        L(cm^-2 s^-1)  reduction factor (hourglass,crossing angle) '
        WRITE(80,*) '----------------------------------------------------------------------------------------'
        CLOSE(80)
      END IF

      IF (writeIntensitySwitch == 1) THEN
        OPEN(unit=80,file='intensity.out')
        WRITE(80,*)' sim.turn     t(hours)        N1_macro    N1_real        NlostLum1     Sum   NlostDebunch1  Sum   NLostBet1&
             &   Sum   NlostMom1     Sum     N2_macro&
             &     N2_real     NlostLum2     Sum   NlostDebunch2  Sum&
             & NLostBet2   Sum   NlostMom2     Sum'
        WRITE(80,*) ' ---------------------------------------------------------------------------------------------------------&
             &-----------------------------------------------------------------------------&
             &---------------------------------------------------------'
        CLOSE(80)
      END IF

      IF (writeEmittanceSwitch == 1) THEN
        OPEN(unit=80,file='emittance.out')
        WRITE(80,*) 'sim.turn     t(hours)         ex1(m)           ey1(m)&
          &           el1(eV/s/charge) sig_T1           sig_dP/P_2      ex2(m&
          &)           ey2(m)           el2(eV/s/charge) sig_T1           sig&
          &_dP/P_2'
         WRITE(80,*) '-----------------------------------------------------&
          &-----------------------------------------------------------------'
         CLOSE(80)
      END IF
    
      IF (writeIBSSwitch == 1) THEN   
         OPEN(unit=80,file='ibs.out')
         WRITE(80,*) '    sim.turn     t(hours)      Tp(hours)       Tx(hou&
            &rs)        Ty(hours)'
         WRITE(80,*) '-----------------------------------------------------&
            &-----------------------------------------------------------------'
         CLOSE(80)
      END IF

      IF (writeCollImpactsxbetatron == 1) THEN 
        OPEN(unit=80, file='collimator_impacts_x_betatron.out')
        WRITE(80,*)'sim.turn     beam         t(hours)     impact_par(m)    impact_par(sigma)'
        WRITE(80,*) '-------------------------------------------------------------------------'
        CLOSE(80)
      END IF 

      IF (writeCollImpactsybetatron == 1) THEN 
        OPEN(unit=80, file='collimator_impacts_y_betatron.out')
        WRITE(80,*)'sim.turn     beam         t(hours)     impact_par(m)    impact_par(sigma)'
        WRITE(80,*) '-------------------------------------------------------------------------'
        CLOSE(80)
      END IF

      IF (writeCollImpactsxmomentum == 1) THEN 
        OPEN(unit=80, file='collimator_impacts_x_momentum.out')
        WRITE(80,*)'sim.turn     beam         t(hours)     impact_par(m)    impact_par(sigma)'
        WRITE(80,*) '-------------------------------------------------------------------------'
        CLOSE(80)
      END IF

  END subroutine getinput
END MODULE ModGetInput