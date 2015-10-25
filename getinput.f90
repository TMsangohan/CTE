subroutine getinput(pnumber1,pnumber2,emix1,emiy1,emix2,emiy2,tauhat1,tauhat2,rmsDelta1,rmsBunchLen1,rmsDelta2,rmsBunchLen2)
    include 'collider_time_evolution_common.f90'
    real :: pnumber1,pnumber2,emix1,emiy1,emix2,emiy2,tauhat1,tauhat2,&
           eta,tunesexact,rcoeff,tunes,rmsDelta1,rmsBunchLen1,rmsDelta2,&
           rmsBunchLen2,v1,s,dgamma_hat,rIon,CalphaE3C,I2,I4x,I4y,I3,&
           sigEoE0,Cq,alfax,alfay,gamax,gamay,Dx,Dy,DxP,DyP,&
           Hx,Hy,I5x,I5y,Jx,Jy,ki,rhoi

    integer :: k,nfft,m,i,j,IOstatus
    character(4) :: col1,col2,col3,col4,col5,col6,col7,col8,col9,&
          col10,col11,col12,col13,col14
    character(150) :: twissFile,ibsIntFile,gTabFile
!
	! opening input file with simulation settings
    open(unit=10,file='collider_time_evolution.in',status='unknown')

    read(10,*) col1
    write(6,*) col1

! switches to turn on various processes
    read(10,*) RFswitch       !(set to 1 to activate synchrotron motion)
    read(10,*) betatronSwitch !(set to 1 to activate betatron motion)
    read(10,*) raddampSwitch  !(set to 1 to activate radiation damping and quantum excitation)
    read(10,*) IBSswitch      !(set to 1 to activate ibs)
    read(10,*) collimationSwitch !(set to 1 to activate losses on aperture cuts)
    read(10,*) blowupSwitch !(set to 1 to activate artificial blowup - ADT loss maps etc)
    read(10,*) collisionSwitch !(set to 1 to activate luminosity)
    read(10,*) col1

! nturns=number of turns
! nMacro=number of macro particles per bunch
! timeRatio=real machine turns per simulation turn
! nwrite=number of turns between writes of output
! writeAllCoordSwitch = set to 1 to activate a dump of all particle coorinates on every turn where output is written
! writeMountSwitch = set to 1 to activate mountain range output of distribution shape in transv and long planes
! iseed=random number seed
    read(10,*)nturns,nMacro,timeRatio,nwrite,writeAllCoordSwitch,writeMountSwitch,iseed
    write(6,*) '    nturns=',nturns,', # macro part.=',nMacro

! gammat=transition gamma
! circ=ring circumference
! gamma0=gamma of beam
    read(10,*)gammat,circ,gamma0

! vrf=voltage of first RF system (volt)
! nharm=harmonic number of first RF system
! vrf2=voltage of second RF system
! nharm2=harmonic number of first RF system
    read(10,*) vrf,nharm,vrf2,nharm2

! tunex=betatron tune hor.
! tuney=betatron tune ver.
! chromx=chromaticity hor.
! chromy=chromaticity ver.
! dqmin=linear coupling term between planes
! k2L = thin sextupole strength
! k2Lskew = thin skew sextupole strength
    read(10,*)tunex,tuney,chromx,chromy,dqmin,k2L,k2Lskew

! smooth approximation of lattice functions
      betax=circ/(2*pi*tunex)
      betay=circ/(2*pi*tuney)
      write(6,*)
      write(6,*) 'reference beta functions x,y (m) = ',betax,betay
      write(6,*)

! aatom=mass number of ion species
! qatom=charge number of ion species
! [-thib/2:thib/2] is the interval used for losses. particles with abs(t)>thib/2 are considered lost.
!  For single RF, thib=RF period. If thib=0.0, an infinite boundary is assumed and no losses occur longitudinally.
      read(10,*) aatom,qatom,thib
! optics twissfile - used for ibs with ibsMethod='piwLattice', 'modPiwLat', or 'nagaitsev' and for radiation damping with radMethod='lattic'
! files with long. starting conditions, used only if longCoordMethod=1
 89   format(A200)
      read(10,89) twissFile 
      read(10,*) col1

! emixy1= geom. 1 sigma starting emittance bunch 1 (meters)
! emix1, emiy1 = geom. 1 sigma starting emittance bunch 1 (meters)
! pnumber1=starting bunch population bunch1
      read(10,*) emix1,emiy1,pnumber1

! emixy2= geom. 1 sigma starting emittance bunch 2 (metres)
! emix2, emiy2 = geom. 1 sigma starting emittance bunch 1 (metres)
! pnumber1=starting bunch population bunch1
      read(10,*) emix2,emiy2,pnumber2

! longCoordMethod=switch for long. starting cond. Possible values:
!     0: parabolic with smokering,
!     1: read from file
!     2: generate bi-Gaussian, given energy spread and bunch length. Can only be matched in small-angle approximation
!     3: "pseudo-Gaussian",  distribution for hamiltonian rho(h)=1/H*exp(-h/H), so rho(pt,t)~exp(-C1*pt^2-C2*(1-cos(C3*t)). at 2nd order in t this is a bi-gaussian
!        This creates an exactly matched phase-space which does not filament (see Lee eq. 3.29 - to check in accelerator handbook)
! transvCoordMethod : same as for longCoordMethod, but only options 1 or 2
      read(10,*) longCoordMethod, transvCoordMethod

! rmsBunchLen1=rms of bunch length for beam 1 (in metres, used only with longCoordMethod=2 and longCoordMethod=3)
      read(10,*)  rmsBunchLen1,rmsBunchLen2

! rmsDelta1=rms of deltaP/P0 for beam1 (only used with longCoordMethod=2)
      read(10,*) rmsDelta1,rmsDelta2

! tauhat1=half bucket length (seconds) for beam 1 of the RF system with shortest wavelength
! tauhat2=half bucket length (seconds) for beam 2 of the RF system with shortest wavelength
! bunchLenPrecis=precision in sampling of bunch length [upper bound on abs(sampled/wanted -1)]. ex 0.01 gives 1% accuracy
! note: for single RF, tauhat=thib/2
! power=bunch shape parameter, only used with smoke ring
! alint = smoke ring parm, only used with smoke ring
      read(10,*) tauhat1,tauhat2,bunchLenPrecis,power,alint


      read(10,89) filebeam1
      read(10,89) filebeam2
      read(10,*) col1

! radMethod=method to use for radiation damping times. Possible values:
!     manual : damping times and eq. emittances are written in input file
!     lattic : damping times and eq. emittances calculated from radiation integral over lattice. Twiss file required
!     approx : damping times and eq. emittances approximated. Smooth lattice assumed, I4/I2=0. Same formula as in ODE model. rho0 required
      read(10,*) radMethod

! tradlong=longitudinal radiation damping time (seconds). Used if radMethod=manual
! tradperp=transverse radiation damping time. Used if radMethod=manual
! siglong=equilibrium sigma of pt from radiation damping and quantum excitation. Used if radMethod=manual
! sigperp=eqilibrium transverse beam size from radiation damping and quantum excitation (metres). Used if radMethod=manual
      read(10,*) tradlong,tradperp,siglong,sigperp

! rho0=dipole bending raduus in metres. Used if radMethod=approx
      read(10,*) rho0
      read(10,*) col1

! ibsMethod=switch for ibs, possible values: piwiSmooth, piwLattice, modPiwLatt, baneApprox, nagaitsev or interpolat
      read(10,*) ibsMethod

! coupleIBS=switch for transverse IBS coupling (0 gives separate growth rates, 1 gives same alfa in x and y)
! coulombLog used in nagaitsev IBS routine
      read(10,*) coupleIBS,coulombLog

! fracibstot=fraction of IBS strength that should be used
! nBins=number of bins used for long. ibs and mountain range plots
      read(10,*) fracibstot, nBins

      read(10,89) ibsIntFile ! file with interpolated values for IBS, used with ibsMethod='interpolat'
      read(10,89) gTabFile ! file with tabulated values of Bane's g-function. Used only with ibsMethod='baneApprox'
      read(10,*) col1

! refEmxy = reference geometric transverse emittance used for the collimation/aperture cuts in betatron space (metres)
! cutoffAmpl = amplitude in refence sigma =sqrt(beta*refEmxy) at which initial transverse distr is cut (dimensionless)
! collimAvgSwitch = set to 1 to collimate the max amplitude over all phases sqrt(x^2_p^2) and to 0 to collimate only x, thus accounting for the betatron phase
! emitMethod= method used to calculate transverse emittance. Can be:
      ! stdev : just calculate it from the standard deviation of the coordinates
      ! stcal : calculate standard devation, then assumue Gaussian with cut tails. Calculate backwards from the cut the standard dev. of the uncut Gaussian
      ! exfit : calculate betatron actions for all part., fit exponential to get emittance
      ! if tails are not cut, stdev and stcal should be equivalent.
! nSigCutBeta= betatron cut in sigma (applied both in x and y. no skew implemented---this complicates emittance calculation) (dimensionless)
! nSigCutMom= number of betatron sigma at which the momentum collimator is placed
! betaxMom= horizontal beta function at momentum collimator (metres)
! dispxMom= horizontal dispersion
      read(10,*) refEmxy,cutoffAmpl,collimAvgSwitch,emitMethod
      read(10,*) nSigCutBeta,nSigCutMom,betaxMom,dispxMom
      read(10,*) col1

      xcut=cutoffAmpl*sqrt(betax*refEmxy) ! cutoff in metres of the transverse distributions x=sqrt(beta refemxy) for matched distribution
      ycut=cutoffAmpl*sqrt(betay*refEmxy)

! pxKickFac,pyKickFac = kicks given by blowup method to px, py on every turn. unit m (norm. coord.)
      read(10,*) pxKickFac,pyKickFac,blowupMethod
      read(10,*) col1

! collRoutine=choose collision model. possible values: 6a and 1d (see beginning of this file for explanation)
      read(10,*) collRoutine

! nIPs=no. of IPs with different parameters (IPs with identical parameters are given with )
! sigI=total interaction cross section in collisions (barn)
      read(10,*) nIPs,sigI

! nbunches=number of bunches
! longIntBins=number of bins for numeric integration of hourglass effect
! angleSwitch=switch to change crossing angle plane on odd turns (1 alternating angle, 0 not alternating)
      read(10,*) nbunches,longIntBins,angleSwitch

! betaS=beta* (metres)
! theta=half crossing angle (rad)
! IPmult=number of IPs with this beta* and crossing angle
      do m=1,nIPs               !loop over all IPs, read beta*, crossing angle/2 and multiplicity
         read(10,*) betaS(m),theta(m),IPmult(m)
      enddo

! calculate trev, beta etc
      beta= sqrt(1-1/gamma0**2)
      vrev = clight*beta
      trev = circ/vrev
      trf = trev/nharm ! tref \approx 2.5E-9
      frf = 1/trf !frf = 400 MHz
      omegarf = 2*pi/trf
      eta = 1/gammat**2 - 1/gamma0**2 ! slipfactor see Lee (2.162)

! calculate equivalent simulation time
      eqTime=nturns*trev*timeRatio/3600.
      write(6,*)'    equivalent time (hours) = ', eqtime
      write(6,*)
      open(unit=34,file='eqtime.out')
      write(34,*)' equivalent time (hours) = ', eqtime

! write summary on screen
      if (RFswitch.eq.1) then
         write(6,*) '    RF motion is:         ON'
      else
         write(6,*) '    RF motion is:         OFF'
      endif
      write(6,*)
      if (betatronSwitch.eq.1) then
         write(6,*) '    betatron motion is:   ON'
         write(6,*)'        -->reference beta functions:'
         write(6,*) '          circ/(2*pi*tunex) = ',circ/(2*pi*tunex),',  circ/(2*pi*tuney) = ',circ/(2*pi*tuney)

      else
         write(6,*) '    betatron motion is:   OFF'
      endif
      write(6,*)
      if (collimationSwitch.eq.1) then
         write(6,*) '    Collimation is:        ON'
         write(6,*) '       --> number of aperture limitations:',nAperCuts
      else
         write(6,*) '    Collimation is:       OFF'
      endif

      if (blowupSwitch.eq.1) then
         write(6,*) '    Blowup is:             ON'
      else
         write(6,*) '    Blowup is:            OFF'
      endif


      write(6,*)
      if (collisionSwitch.eq.1) then
         write(6,*) '    Collisions are:       ON'
         write(6,*) '       --> using collision routine:',collRoutine
         do m=1,nIPs
            write(6,*)'       -->',IPmult(m),' IPs with beta*=',betaS(m), 'm, phi/2=',theta(m)
         enddo

      else
         write(6,*) '    Collisions are:       OFF'
      endif

      close(10)
      write(6,*)
! makesure nresamp is a power of 2 ??? where is this nresamp defined the first time ???
      nfft=2
      do k=1,1000
         nfft = 2*nfft
         if(nfft.ge.nresamp)go to 22
      enddo
 22   continue
!      write(6,*)'nresamp,nfft',nresamp,nfft
      nresamp=nfft

! calculate trev, beta etc
      beta= sqrt(1-1/gamma0**2)
      vrev = 2.998e8*beta
      trev = circ/vrev
      trf = trev/nharm
      frf = 1/trf
      omegarf = 2*pi/trf
      eta = 1/gammat**2 - 1/gamma0**2
!      write(6,*)'eta = ',eta
! longitudinal eom, t(k) = arrival time with synchronous point at pi rf radians
!                  pt(k) = gamma-gamma0
! d pt(k)/dn = e*qatom*vrf*omegarf*(t(k)-trf/2)/(aatom*pmass*clight**2)
! d  t(k)/dn = trev*eta*pt(k)/(beta0**2*gamma0)
! coefficients  d pt(k)/dn = pcoeff*(t(k)-trf/2) , d  t(k)/dn = tcoeff*pt(k)
      if (vrf.eq.0) then
         v1=vrf2
      else
         v1=vrf2
      endif
      pcoeff = qatom*v1*omegarf/(938.e6*aatom)
      tcoeff = trev*eta/(beta**2*gamma0)
      rcoeff = pcoeff/tcoeff
      if(rcoeff.gt.0)then
         write(6,*)'wrong sign of voltage for sign of eta'
         stop
      endif
!      write(6,*)' tuney,mode',tuney,mode
! get the frequency of the slowest betatron sideband
!      write(6,*)'lowest sideband (MHz) = ',1.e-6*(tuney-mode)/trev
! given gap volts get (gamma-gamma0) corresponding to amplitude in arrival time variation
      dgamma_hat = tauhat1*sqrt(-rcoeff)
!tauhat1=half bucket length (seconds) for beam 1 of the RF system with shortest wavelength
!      write(6,*)' amp of de/e = ',dgamma_hat/gamma0
!      write(6,*)'max change in tau per turn = ',dgamma_hat*tcoeff
! synchrotron tune
      tunes = sqrt(-pcoeff*tcoeff)/(2*pi)
!      write(6,*) tunes,pcoeff,vrf,omegarf,qatom
!      stop
      tunesexact = acos(1+pcoeff*tcoeff/2)/(2*pi)
!      write(6,*)' synchrotron frequency = ',tunes/trev
!      write(6,*)' revolution period = ',trev

!      write(6,*)' exact synchrotron tune = ',tunesexact
! numerical parms
!      write(6,*)'tbin/sigtsmooth=', 2.506628*(thib-tlob)/
!     : (taupart*nresamp)
      call system('rm intensity.out')
      call system('rm luminosity.out')
      call system('rm ibs.out')
      call system('rm mountain*.out')
      call system('rm emittance.out')
      call system('rm eqtime.out')
      call system('rm collimator_impacts_x_betatron.out')
      call system('rm collimator_impacts_x_momentum.out')
      call system('rm collimator_impacts_y_betatron.out')



!      open(unit=32,file='tran.xvz',status='unknown')
!      open(unit=33,file='csmon.out',status='unknown')
!      open(unit=52,file='nLostBurnoff.out',status='unknown')

!      open(unit=34,file='tran.cs',status='unknown')
!      open(unit=44,file='tran.sig',status='unknown')

!      call getimped !! used for cooling - uncomment later
! get fiducial normalizations
!      write(6,*)' dpwake for 1 meter offset with step wake = 1'
!      write(6,*)(circ/(2*pi*tuney))*1.6e-19*qatom*pnumber1/
!     :          (beta*gamma0*938.e6*aatom/qatom)

! normalization for longitudinal profiles
      nwrite0 = nwrite
! equivalent time in machine of simulation, rescaling with number of sim turns per real turns


!     read mad twiss file with lattice (for some methods of ibs and radiation damping only). Assume fixed format of tfs table for simplicity
      if (((raddampSwitch.eq.1).and.(radMethod.ne."manual")).or.&
           ((ibsSwitch.eq.1).and.((ibsMethod.eq.'piwLattice').or.&
           (ibsMethod.eq.'modPiwLatt').or.(ibsMethod.eq.'baneApprox').or.(ibsMethod.eq.'nagaitsev')))) then
         open(unit=99,file=trim(twissFile))
         col1="||"
!        read until header of twiss file

         do i=1,45
            read(99,*) col1
         enddo
         write(6,*) '    Reading twiss file ',trim(twissFile)
         write(6,*)
         read(99,*) col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14 ! read column headers

!        check that we have the right columns in twiss file
         if((col2=="NAME").and.(col3=="S   ").and.(col4=="L   ").and.(col5=="BETX").and.(col6=="BETY")&
              .and.(col7=="ALFX").and.(col8=="ALFY").and.(col9=="DX  ").and.(col10=="DPX ")&
              .and.(col11=="DY  ").and.(col12=="DPY ").and.(col13=="ANGL").and.(col14=="K1L ")) then
         else
            write(6,*) "Error in format of twiss file ", trim(twissFile)
            write(6,*) "Columns must be NAME|S   |L   |BETX|BETY|ALFX|ALFY|DX  |DPX |DY  |DPY |ANGLE|K1L"
            stop
         endif
         read(99,*) col1        ! read one more line to get to the actual table

!     read in lattice
         nElem=1                ! count the number of elements in lattice
         IOstatus=0             ! flag to check when end of file is reached

         do while(IOstatus.eq.0) ! read elements until end of file
            read(99,*,iostat=IOstatus) &
                 col1,S,leng(nElem),betx(nElem),bety(nElem),alfx(nElem),alfy(nElem),dispx(nElem),dispxP(nElem),dispy(nElem),&
                 dispyP(nElem),bendAng(nElem),K1L(nElem)
            if(leng(nElem)/=0) nElem=nElem+1
         enddo
         nElem=nElem-1
         close(99)
      endif



      if (ibsSwitch.eq.1) then
         write(6,*) '    IBS is:               ON'
         write(6,*) '       -->using method: ',ibsMethod

         if (ibsMethod.eq.'baneApprox') then
!           read in tabulated values of Bane's g-function
            open(unit=99,file=gTabFile)
            do i=1,1501
               read(99,*) gTab(i,1),gTab(i,2)
            enddo
            close(99)
         endif

         if (ibsMethod.eq.'interpolat') then
!           read in tabulated function values for ibs
            open(unit=10,file=trim(ibsIntFile))
            write(6,*) '       -->reading tabulated IBS values from ',trim(ibsIntFile)
            read(10,*) xMin,xMax,xBins
            read(10,*) yMin,yMax,yBins
            read(10,*) zMin,zMax,zBins
            read(10,*) pnumInterp
            do i=1,xBins
               do j=1,yBins
                  do k=1,zBins
                     read(10,*) Ap(i,j,k),Ax(i,j,k),Ay(i,j,k)
                  enddo
               enddo
            enddo
            close(10)
         endif

      else
         write(6,*) '    IBS is:               OFF'
      endif

      write(6,*)
      if (raddampSwitch.eq.1) then
         write(6,*) '    Radiation damping is: ON'
         write(6,*) '       -->using method: ',radMethod
      else
         write(6,*) '    Radiation damping is: OFF'
      endif


      if((raddampSwitch.eq.1).and.(radMethod.ne."manual")) then
!        calcualte radiation damping times and equilibrium emittances if not given manually
!        Reference: Chao, Tigner: Handbook of Accelerator physics and engineering, (1998) page 186

         if (radMethod.eq."approx") then
!        approx: smooth lattice, assume I4_j<<I2 for all j. see handbook p186 and put I4_j=0, I2=2*pi/rho0 and c_alpha=rIon/(3*c^5*mIon^3)
!        additional div by 2 compared to handbook to get damping time for emittance instead of amplitude. See  paper derivation.

            alfax=0.
            alfay=0.
            gamax=(1.+alfax**2)/betax
            gamay=(1.+alfay**2)/betay
            Dx=circ/(2*pi*gammat**2)
            Dy=0.
            DxP=0.1             ! should find an approximation formula. However not very important
            DyP=0.
            Hx=(betax*DxP+2*alfax*Dx*DxP+gamax*Dx)
            Hy=(betay*DyP+2*alfay*Dy*DyP+gamay*Dy)

!     define smooth approximation of radiation integrals
            I2=2*pi/rho0
            I3=2*pi/rho0**2
            I4x=0.0
            I4y=0.0
            I5x=Hx*2*pi/rho0**2
            I5y=Hy*2*pi/rho0**2

         elseif(radMethod.eq."lattic") then
!     calculate radiation integrals over lattice using twiss file
            I2=0.
            I3=0.
            I4x=0.
            I4y=0.
            I5x=0.
            I5y=0.

            do i=1,nElem
               if(bendAng(i).ne.0.) then
                  rhoi=leng(i)/bendAng(i)
                  ki=K1L(i)/leng(i)
                  I2=I2+1./rhoi**2*leng(i)
                  I3=I3+1./rhoi**3*leng(i)
                  I4x=I4x+(dispx(i)/rhoi**3*(1.+2.*rhoi**2*ki)+2.*dispx(i)*ki/rhoi)*leng(i)
                  I4y=I4y+(dispy(i)/rhoi**3*(1.-2.*rhoi**2*ki)-2.*dispy(i)*ki/rhoi)*leng(i)
                  gamax=(1.+alfx(i)**2)/betx(i)
                  gamay=(1.+alfy(i)**2)/bety(i)
                  Hx=betx(i)*dispxP(i)**2 + 2*alfx(i)*dispx(i)*dispxP(i)+gamax*dispx(i)**2
                  Hy=bety(i)*dispyP(i)**2 + 2*alfy(i)*dispy(i)*dispyP(i)+gamay*dispy(i)**2
                  I5x=I5x+Hx*2*pi/rho0**2*leng(i)
                  I5y=I5y+Hy*2*pi/rho0**2*leng(i)
               endif
            enddo
            write(6,*) "          calculating radiation integrals: "
            write(6,*) "          I2=",I2,"I3=",I3,"I4x=",I4x,"I4y=",I4y,"I5x=",I5x,"I5y=",I5y
         endif

         rIon=(qatom**2/aatom)*1.54e-18
         CalphaE3C=rIon*gamma0**3/(3*trev)  ! C_alpha*E^3/C in handbook formulas p 186

         tradperp=1./(CalphaE3C*I2*(1.-I4x/I2))/2. ! eq 9, but div by 2 to get rise time for emittance and not amplitude
         tradlong=1./(CalphaE3C*I2*(2.+(I4x+I4y)/I2))/2.

         Cq=55./(32.*sqrt(3.0))*1.0546e-34*clight / (938.e6*aatom*1.60218e-19)
         sigEoE0=Cq*gamma0**2*I3/(2*I2+I4x+I4y)  ! eq 18
         siglong=gamma0*sigEoE0 ! sigma of pt
         Jx=1.-I4x/I2
         Jy=1.-I4y/I2
         sigperp=Cq*gamma0**2*I5x/(Jx*I2)
         write(6,*) '          calculating radiation damping times and equilibrium beam sizes:'
         write(6,*) '          Tx=',tradperp/3600,'h,  Ts=',tradlong/3600,'h,  sigPt=',siglong,'ex=',sigperp,'m'
      endif

      write(6,*)
      write(6,*)'... collider_time_evolution.in read'
      write(6,*)

! open output files, write headers
      open(unit=80,file='luminosity.out')
      write(80,*) 'sim.turn     t(hours)        L(cm^-2 s^-1)  reduction factor (hourglass,crossing angle) '
      write(80,*) '----------------------------------------------------------------------------------------'
      close(80)

      open(unit=80,file='intensity.out')
!     write(80,*)' sim.turn     t(hours)        N1_macro    N1_real        NlostLum1     Sum   NlostDebunch1  Sum   N2_macro&
      write(80,*)' sim.turn     t(hours)        N1_macro    N1_real        NlostLum1     Sum   NlostDebunch1  Sum   NLostBet&
           &1   Sum   NlostMom1     Sum     N2_macro&
           &     N2_real     NlostLum2     Sum   NlostDebunch2  Sum&
           &NLostBet2   Sum   NlostMom2     Sum'
!      write(80,*) 'sim.turn     t(hours)        N1_macro    N1_real     &
!     &    N2_macro     N2_real'
      write(80,*) ' ---------------------------------------------------------------------------------------------------------&
           &-----------------------------------------------------------------------------&
           &---------------------------------------------------------'
      close(80)

      open(unit=80,file='emittance.out')
      write(80,*) 'sim.turn     t(hours)         ex1(m)           ey1(m)&
     &           el1(eV/s/charge) sig_T1           sig_dP/P_2      ex2(m&
     &)           ey2(m)           el2(eV/s/charge) sig_T1           sig&
     &_dP/P_2'
      write(80,*) '-----------------------------------------------------&
     &-----------------------------------------------------------------'
      close(80)

      open(unit=80,file='ibs.out')
      write(80,*) '    sim.turn     t(hours)      Tp(hours)       Tx(hou&
     &rs)        Ty(hours)'
      write(80,*) '-----------------------------------------------------&
     &-----------------------------------------------------------------'
      close(80)


      open(unit=80, file='collimator_impacts_x_betatron.out')

      write(80,*)'sim.turn     beam         t(hours)     impact_par(m)    impact_par(sigma)'

      write(80,*) '-------------------------------------------------------------------------'
      close(80)

      open(unit=80, file='collimator_impacts_y_betatron.out')

      write(80,*)'sim.turn     beam         t(hours)     impact_par(m)    impact_par(sigma)'

      write(80,*) '-------------------------------------------------------------------------'
      close(80)

      open(unit=80, file='collimator_impacts_x_momentum.out')

      write(80,*)'sim.turn     beam         t(hours)     impact_par(m)    impact_par(sigma)'

      write(80,*) '-------------------------------------------------------------------------'
      close(80)

      end