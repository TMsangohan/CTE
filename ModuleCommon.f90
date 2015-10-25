! ---------------------------------------------------------------------------------------------------------------
! AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 1.0 : 
!   AUTHOR    : TOM MERTENS
!   DATE      : 11/06/2015
!   COPYRIGHT : CERN
!
!   DESCRIPTION : 
!       DEFINES COMMON (GLOBAL) VARIABLES 
! ---------------------------------------------------------------------------------------------------------------

MODULE commondata

IMPLICIT NONE

! ---------------------------------------------------------------------------------------------------------------
! ***************************************************************************************************************
!
! CODE DESCRIPTION: SELECTED_INT_KIND(R)
! **************************************
!
! RETURN THE KIND VALUE OF THE SMALLEST INTEGER TYPE THAT CAN REPRESENT ALL VALUES 
! RANGING FROM -10^R (EXCLUSIVE) TO 10^R (EXCLUSIVE). 
! IF THERE IS NO INTEGER KIND THAT ACCOMMODATES THIS RANGE, SELECTED_INT_KIND 
!
! RETURN VALUE
! ------------
! -1 = ALLOWS FOR SAFETY CHECKS. 
!
! ***************************************************************************************************************
! ---------------------------------------------------------------------------------------------------------------

      INTEGER,      PARAMETER :: I4B   = SELECTED_INT_KIND(9)
      INTEGER,      PARAMETER :: I1B   = SELECTED_INT_KIND(2)

! ---------------------------------------------------------------------------------------------------------------
! ***************************************************************************************************************
!
! CODE DESCRIPTION: SELECTED_REAL_KIND(P,R,RADIX)
! ***********************************************
!
! RETURNS THE KIND VALUE OF A REAL DATA TYPE WITH 
! DECIMAL PRECISION OF AT LEAST P DIGITS, EXPONENT RANGE OF AT LEAST R, AND WITH A RADIX OF RADIX. 
!
! RETURN VALUE:
! -------------
! SELECTED_REAL_KIND RETURNS THE VALUE OF THE KIND TYPE PARAMETER OF A REAL DATA TYPE WITH 
! DECIMAL PRECISION OF AT LEAST P DIGITS, A DECIMAL EXPONENT RANGE OF AT LEAST R, AND WITH THE REQUESTED RADIX. 
! IF THE RADIX PARAMETER IS ABSENT, REAL KINDS WITH ANY RADIX CAN BE RETURNED. 
! IF MORE THAN ONE REAL DATA TYPE MEET THE CRITERIA, THE KIND OF THE DATA TYPE WITH THE SMALLEST DECIMAL PRECISION IS RETURNED. 
! IF NO REAL DATA TYPE MATCHES THE CRITERIA, THE RESULT IS
!
! -1 IF THE PROCESSOR DOES NOT SUPPORT A REAL DATA TYPE WITH A
!     PRECISION GREATER THAN OR EQUAL TO P, BUT THE R AND RADIX REQUIREMENTS CAN BE FULFILLED 
!
! -2 IF THE PROCESSOR DOES NOT SUPPORT A REAL TYPE WITH AN EXPONENT
!     RANGE GREATER THAN OR EQUAL TO R, BUT P AND RADIX ARE FULFILLABLE 
!
! -3 IF RADIX BUT NOT P AND R REQUIREMENTS ARE FULFILLABLE 
!
! -4 IF RADIX AND EITHER P OR R REQUIREMENTS ARE FULFILLABLE 
!
! ***************************************************************************************************************
! ---------------------------------------------------------------------------------------------------------------

      INTEGER,         PARAMETER :: SP    = SELECTED_REAL_KIND(6,9)
      INTEGER,         PARAMETER :: DP    = SELECTED_REAL_KIND(6,15)


! ------------------------------------------------------------------------------------------------------------------
! ******************* DESCRIPTION VARIABLES*************************************************************************
! ------------------------------------------------------------------------------------------------------------------
!
!  *** INTEGER KIND I1B ********************************************************************************************
! RFswitch          : SET TO 1 TO ACTIVATE SYNCHROTRON MOTION
! betatronSwitch    : SET TO 1 TO ACTIVATE BETATRON MOTION
! raddampSwitch     : SET TO 1 TO ACTIVATE RADIATION DAMPING AND QUANTUM EXCITATION
! IBSswitch         : SET TO 1 TO ACTIVATE IBS
! collimationSwitch : SET TO 1 TO ACTIVATE LOSSES ON APERTURE CUTS
! blowupSwitch      : SET TO 1 TO ACTIVATE ARTIFICIAL BLOWUP - ADT LOSS MAPS ETC
! collisionSwitch   : SET TO 1 TO ACTIVATE LUMINOSITY
! levellingSwitch   : SET TO 1 TO ACTIVATE BETA LEVELLING OF LUMINOSITY IN ONE OR MORE IPS
! writeEqTimeSwitch : SET TO 1 TO WRITE REAL TIME TO SEPARATE FILE
! writeScreenSummarySwitch : SET TO 1 TO WRITE A SUMMARY OF THE SIMULATION SETTINGS TO THE SCREEN AT START OF SIMULATION
! longCoordMethod   : SWITCH FOR LONGITUDINAL STARTING CONDITIONS POSSIBLE VALUES ARE
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       0: PARABOLIC WITH SMOKERING,
!       1: READ FROM FILE
!       2: GENERATE BI-GAUSSIAN, GIVEN ENERGY SPREAD AND BUNCH LENGTH. CAN ONLY BE MATCHED IN SMALL-ANGLE APPROXIMATION
!       3: "PSEUDO-GAUSSIAN", DISTRIBUTION FOR HAMILTONIAN RHO(H)=1/H*EXP(-H/H), SO RHO(PT,T)~EXP(-C1*PT^2-C2*(1-COS(C3*T)). 
!             AT 2ND ORDER IN T THIS IS A BI-GAUSSIAN THIS CREATES AN EXACTLY MATCHED PHASE-SPACE WHICH 
!             DOES NOT FILAMENT (SEE LEE EQ. 3.29 - TO CHECK IN ACCELERATOR HANDBOOK)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! transvCoordMethod : SAME AS FOR LONGCOORDMETHOD, BUT ONLY OPTIONS 1 OR 2
! coupleIBS         : SWITCH FOR TRANSVERSE IBS COUPLING (0 GIVES SEPARATE GROWTH RATES, 1 GIVES SAME ALFA IN X AND Y)
! collimAvgSwitch   : SET TO 1 TO COLLIMATE THE MAX AMPLITUDE OVER ALL PHASES sqrt(x^2_p^2) AND TO 0 TO COLLIMATE ONLY X, THUS ACCOUNTING FOR THE BETATRON PHASE
! nIPs              : NUMBER OF IPS WITH DIFFERENT PARAMETERS (IPS WITH IDENTICAL PARAMETERS ARE GIVEN WITH )
! angleSwitch       : SWITCH TO CHANGE CROSSING ANGLE PLANE ON ODD TURNS (1 ALTERNATING ANGLE, 0 NOT ALTERNATING)
! IPmult            : NUMBER OF IPS WITH THIS BETA* AND CROSSING ANGLE
! levIPSwitch       : DEFINES IF THE IP IS LEVVELED (SET TO ONE) OR NOT
!
!  *** INTEGER KIND I4B ********************************************************************************************
! nturns            : NUMBER OF TURNS TO BE SIMULATED
! nMacro            : NUMBER OF MACRO PARTICLES PER BUNCH
! timeRatio         : NUMBER OF REAL MACHINE TURNS PER SIMULATION TURN
! nwrite            : NUMBER OF TURNS BETWEEN WRITING TO FILE OF GENERATED OUTPUT
! iseed             : RANDOM NUMBER SEED
! nBins             : NUMBER OF BINS USED FOR LONGITUDINAL IBS AND MOUNTAIN RANGE PLOTS
! nbunches          : NUMBER OF BUNCHES IN THE BEAMS
! longIntBins       : NUMBER OF BINS FOR NUMERIC INTEGRATION OF HOURGLASS EFFECT
! nElem             : NUMBER OF ELEMENTS IN THE TWISS FILE
! xbins             : NUMBER OF BINS IN THE DATA RANGE IN THE IBS INTERPOLAT FILE
! ybins             : NUMBER OF BINS IN THE DATA RANGE IN THE IBS INTERPOLAT FILE
! zbins             : NUMBER OF BINS IN THE DATA RANGE IN THE IBS INTERPOLAT FILE
!
!  *** REAL KIND SP ************************************************************************************************
! gammat            : TRANSITION (RELATIVISTIC) GAMMA OF THE MACHINE
! circ              : RING CIRCUMFERENCE
! Gamma0            : (RELATIVISTIC) GAMMA OF BEAM
! vrf               : VOLTAGE OF FIRST RF SYSTEM (VOLT) - SET VOLTAGE TO ZERO OF ONLY ONE SYSTEM AVAILABLE
! nharm             : HARMONIC NUMBER OF FIRST RF SYSTEM
! vrf2              : VOLTAGE OF SECOND RF SYSTEM
! nharm2            : HARMONIC NUMBER OF SECOND RF SYSTEM
! tunex             : BETATRON TUNE HOR.
! tuney             : BETATRON TUNE VER.
! chromx            : CHROMATICITY HOR.
! chromy            : CHROMATICITY VER.
! dqmin             : LINEAR COUPLING TERM BETWEEN TRANSVERSE PLANES
! k2L               : THIN SEXTUPOLE STRENGTH
! k2Lskew           : THIN SKEW SEXTUPOLE STRENGTH
! aatom             : MASS NUMBER OF ION SPECIES
! qatom             : CHARGE NUMBER OF ION SPECIES
! thib              : [-thib/2:thib/2] IS THE INTERVAL USED FOR LOSSES. PARTICLES WITH ABS(T)>THIB/2 ARE CONSIDERED LOST.
!                     FOR SINGLE RF, THIB=RF PERIOD. IF THIB = 0.0, AN INFINITE BOUNDARY IS ASSUMED AND NO LOSSES OCCUR LONGITUDINALLY.
! bunchLenPrecis    : PRECISION IN SAMPLING OF BUNCH LENGTH [UPPER BOUND ON ABS(SAMPLED/WANTED -1)]. EX 0.01 GIVES 1% ACCURACY
! power             : BUNCH SHAPE PARAMETER, ONLY USED WITH SMOKE RING
! alint             : SMOKE RING PARMAMETER, ONLY USED WITH SMOKE RING
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FOR MANUAL RADIATION DAMPING
!
! tradlong          : LONGITUDINAL RADIATION DAMPING TIME (SECONDS). USED IF RADMETHOD=MANUAL
! tradperp          : TRANSVERSE RADIATION DAMPING TIME. USED IF RADMETHOD=MANUAL
! siglong           : EQUILIBRIUM SIGMA OF PT FROM RADIATION DAMPING AND QUANTUM EXCITATION. USED IF RADMETHOD=MANUAL
! sigperp           : EQUILIBRIUM TRANSVERSE BEAM SIZE FROM RADIATION DAMPING AND QUANTUM EXCITATION (METRES). USED IF RADMETHOD=MANUAL
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%
! rho0              : DIPOLE BENDING RADUUS IN METRES. USED IF RADMETHOD=APPROX
! COULOMBLOG        : VALUE USED IN NAGAITSEV IBS ROUTINE
! fracibstot        : FRACTION OF IBS STRENGTH THAT SHOULD BE USED ALLOWING SOME FLEXIBILITY IN THE MODELLING
! refEmxy           : REFERENCE GEOMETRIC TRANSVERSE EMITTANCE USED FOR THE COLLIMATION/APERTURE CUTS IN BETATRON SPACE (METRES) - COLLIMATION
! cutoffAmpl        : AMPLITUDE IN REFERENCE sigma =sqrt(beta*refemxy) AT WHICH INITIAL TRANSVERSE DISTR IS CUT (DIMENSIONLESS) - COLLIMATION
! nSigCutBeta       : BETATRON CUT IN SIGMA (APPLIED BOTH IN X AND Y. NO SKEW IMPLEMENTED---THIS COMPLICATES EMITTANCE CALCULATION) (DIMENSIONLESS) - COLLIMATION
! nSigCutMom        : NUMBER OF BETATRON SIGMA AT WHICH THE MOMENTUM COLLIMATOR IS PLACED - COLLIMATION
! betaxMom          : HORIZONTAL BETA FUNCTION AT MOMENTUM COLLIMATOR (METRES) - COLLIMATION
! dispxMom          : HORIZONTAL DISPERSION - COLLIMATION
! xcut              : CUTOFF IN METRES OF THE TRANSVERSE DISTRIBUTIONS X=SQRT(BETA REFEMXY) FOR MATCHED DISTRIBUTION
! ycut              : CUTOFF IN METRES OF THE TRANSVERSE DISTRIBUTIONS X=SQRT(BETA REFEMXY) FOR MATCHED DISTRIBUTION
! pxKickFac         : KICKS GIVEN BY BLOWUP METHOD TO PX ON EVERY TURN. UNIT M (NORM. COORD.)
! pyKickFac         : KICKS GIVEN BY BLOWUP METHOD TO PY ON EVERY TURN. UNIT M (NORM. COORD.)
! sigI              : TOTAL INTERACTION CROSS SECTION IN COLLISIONS (BARN)
! beta              : RELATIVISTIC BETA 
! vrev              : RELATIVISTIC VELOCITY OF THE REFERENCE PARTICLE
! trev              : REVOLUTION TIME OF THE REFERENCE PARTICLE
! trf               : RF PERIOD TAKING INTO ACCOUNT THE CHOSEN HARMONIC NUMBER I.E. RETURNING THE PERIOD OF THE BUCKETS
! frf               : RF FREQUENCY
! omegarf           : RF ANGULAR FREQUENCY
! eqTime            : REAL TIME ELAPSED EXPRESSED IN HOURS
! xmin              : BOUNDARY OF TABLE IN PRE-CALCULATED IBS GROWTH RATES FILE FOR INTERPOLAT
! xmax              : BOUNDARY OF TABLE IN PRE-CALCULATED IBS GROWTH RATES FILE FOR INTERPOLAT
! ymin              : BOUNDARY OF TABLE IN PRE-CALCULATED IBS GROWTH RATES FILE FOR INTERPOLAT
! ymax              : BOUNDARY OF TABLE IN PRE-CALCULATED IBS GROWTH RATES FILE FOR INTERPOLAT
! zmin              : BOUNDARY OF TABLE IN PRE-CALCULATED IBS GROWTH RATES FILE FOR INTERPOLAT
! zmax              : BOUNDARY OF TABLE IN PRE-CALCULATED IBS GROWTH RATES FILE FOR INTERPOLAT
! pnumInterp        : NUMBER OF INTERPOLATION POINTS IN THE IBS INTERPOLAT FILE
!
! *** CHARACTER DIMENSION 2 ***************************************************************************************
! collRoutine       : CHOOSE COLLISION MODEL. POSSIBLE VALUES: 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!             6a
!             1d 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
! *** CHARACTER DIMENSION 5 ***************************************************************************************
! emitMethod        : METHOD USED TO CALCULATE TRANSVERSE EMITTANCE. CAN BE
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!           stdev : JUST CALCULATE IT FROM THE STANDARD DEVIATION OF THE COORDINATES
!           stcal : CALCULATE STANDARD DEVATION, THEN ASSUMUE GAUSSIAN WITH CUT TAILS. CALCULATE BACKWARDS FROM THE CUT THE STANDARD DEV. OF THE UNCUT GAUSSIAN
!           exfit : CALCULATE BETATRON ACTIONS FOR ALL PART., FIT EXPONENTIAL TO GET EMITTANCE
!           IF TAILS ARE NOT CUT, STDEV AND STCAL SHOULD BE EQUIVALENT.
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! blowupMethod      : METHOD USED FOR ARTIFICIAL BLOWUP CHOICES ARE
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!           unifo
!           gauss
!           unSum
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! *** CHARACTER DIMENSION 6 ***************************************************************************************
!   radMethod       : METHOD TO USE FOR RADIATION DAMPING TIMES. POSSIBLE VALUES ARE
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       manual : DAMPING TIMES AND EQ. EMITTANCES ARE WRITTEN IN INPUT FILE
!       lattic : DAMPING TIMES AND EQ. EMITTANCES CALCULATED FROM RADIATION INTEGRAL OVER LATTICE. TWISS FILE REQUIRED
!       approx : DAMPING TIMES AND EQ. EMITTANCES APPROXIMATED. SMOOTH LATTICE ASSUMED, I4/I2=0. SAME FORMULA AS IN ODE MODEL. RHO0 REQUIRED
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! *** CHARACTER DIMENSION 10 ***************************************************************************************
!   ibsMethod   = SWITCH FOR IBS, POSSIBLE VALUES ARE
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!               PIWISMOOTH 
!               PIWLATTICE 
!               MODPIWLATT
!               BANEAPPROX 
!               NAGAITSEV
!               INTERPOLAT
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! betaS         : BETA* (METRES)
! theta         : HALF CROSSING ANGLE (RAD)
! lumimax       : MAXIMUM LUMINOSITY FOR LEVELLING
! betaSMin      : KEEPS TRACK OF THE ORIGINAL BETA STAR INPUT AT THE RESPECTIVE IP, TO KEEP TRACK DURIING BETA-LEVELLING
!
! *** CHARACTER DIMENSION 10 ***************************************************************************************
! leng          : LENGTH OF ELEMENT IN TWISS FILE
! betx          : HORIZONTAL OPTICAL BETA AT EACH ELEMENT IN TWISS FILE
! bety          : VERTICAL OPTICAL BETA AT EACH ELEMENT IN TWISS FILE
! dispx         : HORIZONTAL DISPERSION AT EACH ELEMENT IN TWISS FILE
! alfx          : VERTICAL OPTICAL ALFA AT EACH ELEMENT IN TWISS FILE - REMEMBER ALFA = BETA'
! alfy          : HORIZONTAL OPTICAL ALFA AT EACH ELEMENT IN TWISS FILE - REMEMBER ALFA = BETA'
! dispxP        :
! dispy         : VERTICAL DISPERSION AT EACH ELEMENT IN TWISS FILE
! dispyP        :
! bendAng       : BENDING ANGLE GENERATED BY EACH ELEMENT
! K1L           : SEXTUPOLE K1 NORMAL
! K1S           : SEXTUPOLE K1 SKEW
!
! *** CHARACTER DIMENSION 120 ***************************************************************************************
! Ap            : IBS INTERPOLAT GROWTH RATES
! Ax            : IBS INTERPOLAT GROWTH RATES
! Ay            : IBS INTERPOLAT GROWTH RATES
! ------------------------------------------------------------------------------------------------------------------
! ******************* DESCRIPTION VARIABLES **** END ***************************************************************
! ------------------------------------------------------------------------------------------------------------------

      INTEGER(I4B), PARAMETER :: nb=6000000,nE=13500,nIbs=120

      REAL(SP),     PARAMETER :: PI    = 3.1415926536_sp
      REAL(SP),     PARAMETER :: CLIGHT= 2.998e8_sp

      INTEGER(I1B)            :: nIPs,ibsSwitch,collimationSwitch,writeAllCoordSwitch,transvCoordMethod
      INTEGER(I1B)            :: raddampSwitch,writeMountSwitch,blowupSwitch,coupleIBS,writeEqTimeSwitch
      INTEGER(I1B)            :: collisionSwitch,angleSwitch,collimAvgSwitch,RFswitch,betatronSwitch,writeIntensitySwitch
      INTEGER(I1B)            :: writeEmittanceSwitch,writeIBSSwitch,writeCollImpactsxbetatron,writeCollImpactsybetatron
      INTEGER(I1B)            :: writeCollImpactsxmomentum
      INTEGER(I1B)            :: levellingSwitch,longCoordMethod,writeScreenSummarySwitch,writeLuminositySwitch
      INTEGER(I4B)            :: kturns,nturns,nMacro,nwrite,nresamp,nperturn, &
                                 nturnon,nharm,nharm2,nBins,nextra,iseed,nwrite0,&
                                 ntdelay,longIntBins,nbunches,&
                                 nElem,xbins,ybins,zbins,nAperCuts
      
      INTEGER(I1B), DIMENSION(10) :: levIPSwitch


      CHARACTER(len=200)            :: filebeam1,filebeam2
      CHARACTER(len=2)              :: collRoutine
      CHARACTER(len=10)             :: ibsMethod
      CHARACTER(len=6)              :: radMethod
      CHARACTER(len=5)              :: emitMethod,blowupMethod

      REAL(SP) :: gammat,circ,gamma0,vrf,trf,omegarf,tcoeff,betax,betay, &
                  pcoeff,frf,tunex,tuney,xinject,vrf2,beta,coulombLog,k2L,k2Lskew,&
                  vrev,trev,chromx,chromy,taupart,aatom,qatom,power,&
                  gain0,gain1,fimped0,fimped1,fracmixbad,thib,&
                  gains0,gains1,fimpeds0,fimpeds1,phipks,cutoffAmpl,&
                  fracibstot,fracibsx,ampbtf,phibtf,harmbtf,&
                  alint,dqmin,phipk,sigI,pxKickFac,pyKickFac,&
                  tradlong,tradperp,siglong,sigperp,timeRatio,&
                  xmin,xmax,ymin,ymax,zmin,zmax,nSigCutBeta,nSigCutMom,dispxMom,&
                  pnumInterp,eqTime,rho0,bunchLenPrecis,refEmxy,refMomSpread,betaxMom,&
                  xcut,ycut !tlob (STILL USED ?)

      REAL(SP), DIMENSION(10)   ::  betaS,theta,IPmult,betaSMin,lumimax
      REAL(SP), DIMENSION(nE)   ::  betx,bety,dispx,alfx,alfy,dispxP,&
                                    leng,dispy,dispyP,bendAng,K1L,K1S

      REAL(SP), DIMENSION(nIbs,nIbs,nIbs) :: Ap,Ax,Ay

      REAL(SP), DIMENSION(1502,2) :: gTab
! ---------------------------------------------------------------------------------------------------------------
! ***************************************************************************************************************
!
! CODE DESCRIPTION: SAVE
! **********************
!
! THE SAVE COMMAND IS NECESSARY TO KEEP THE VALUES OF THE VARIABLES
! ***************************************************************************************************************
! ---------------------------------------------------------------------------------------------------------------

SAVE

END MODULE commondata
