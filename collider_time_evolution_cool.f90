
!compile the programm with
!    gfortran -fdefault-real-8 -fdefault-double-8 -fno-align-commons -mcmodel=medium -o collider_time_evolution_cool collider_time_evolution_cool.f90 

! status: implemented in output to screen and file. but stcal gives wrong results for a 3 sigma cut! Correct though for 2 sigma. this is weird! Can we understand it?

! next steps:
! merge with Tom's Nagaitsev version
! run standard test case (collision ions) to verify that there is no difference
! to model tails, should we introduce form factor also in transverse plane?




      program collider_time_evolution_cool
!     Authors: Roderik Bruce, Mike Blaskiewicz and Tom Mertens
!**************************************************************
!**************************************************************
! Program to track two bunches of macro-particles in time in a collider.
! Subroutines act on the bunches on a turn-by-turn basis. One simulation turn
! can correspond to any chosen number of machine turns.
! The following processes are taken into account:
! --collisions. user can choose between 2 collision routines:
!     --1d: very slow, integrates interaction probability for every particle
!           by sorting particles in opposing beam in discrete bins. 
!           No assumptions on the shape of the beam distribution
!     --6a: fast routine, assumes Gaussian transverse distribution and 
!           calcualtes interaction probability from transverse distribution
!           analytically and uses global reduction factor (hourglass and crossing angle) 
!           for all particles. No assumptions on longitudinal distribution.     
! --IBS. rise time calculated using a standard method and modulated to account for non-Gaussian longitudinal profiles
!   user can choose between the following methods:
!      --Nagaitsev full lattice
!      --Bane full lattice
!      --smooth lattice Piwinski
!      --full lattice Piwinski
!      --full lattice modified Piwinski
!      --interpolation from tabulated risetimes in external file at given points in emittance-space
! --betatron motion
! --synchrotron motion (particles outside RF bucket are lost)
! --radiation damping and quantum excitation
! --transverse aperture cut from collimation


! Input parameters are given in the input file collider_time_evolution.in
! Common blocks are declared in collider_time_evolution_common.f90
! compile the executable with:
! g95 -g -o collider_time_evolution collider_time_evolution.f
!     (as no special libraries are used, also other compilers like gfortran should work)

!     list of changes to implement
!----------------------------------
!     - betatron noise from feedback gives emittance blowup
!     - give actual energy kick in radiation damping algorithm, non-zero synch. phase       
!     - inelasic beam gas
!     - elastic beam gas: both particle loss and emittance blowup
!     - RF noise



!     change log (Roderik Bruce)
!-----------------------
!    *******Version 1.8 (Michaela Schaumann)
!         --> Bug in equations for radiation integrals corrected to equations from handbook, scond edition (2013), p. 220
!         --> had to add K1S (skew quadrupole strength) to MAD-X lattice file for calculation of I4x
!         --> had to add additional variables to common.f90 
!             -> small effects: negligible change in radiation damping times, I5 changed to smaller values, reducing equilibrium emittance and quantum excitation.
!
!    ****version 1.7 (Michaela Schaumann)
!	--> Putting back in the Stochastic cooling method from Mike's original program
!	    - To do that the relation equivalent time = pnumber/nmacro * nturns * trev/3600 has to be fulfilled to treat the stochastic cooling correctly.
!	    - Therefore, change in input parameters: 
!			-> add stochCoolSwitch to switch cooling on and off
!			-> timeRatio is exchanged by eqTime: now timeRatio is calculated in program
!			-> vlaue for nturn will be overwritten if stochastic cooling is on, but still be valid if stochCoolSwitch=0
!			-> nMacro will be the number of macro particles for beam1, nMacro2 will be calculated to fullfill relation above to have the same eqTime for BOTH beams
!	   - Adding additional variables for SC in before main loop and in common file
!	   - in main loop:
!			-> ntruns -> nturns+nextra (where nextra is the number of extra turns without cooling switched on)
!			-> for cooling the particle coordinates have to be saved when measured at the pickup - break up RF and betatron motion for fractions of a turn and stop at the 3 kickers to save coordinates.
!	   - addparticles: add nMacro to argument list, since nMacro1 and nMaro2 are now different
!	   - stransChrom: separate dqmin coupling to external subroutine, since stransChrom is called after each part of the betatron motion, but the coupling has to be done only once per turn
!	   - writemomentsshort: (was not called in the program before), extention to produce additional writeout and call if writeStochCoolSwitch =1
!	--> Adding an rfupdate with fmix =0 right after the collision took place to get rid of large amplitude particles
!   
!   ***version1.6 (Michaela Schaumann)
!       --> Add beta* leveling in collision6a routine
!       --> the input value for betaS is saved in betaSMin
!       --> extend collision6a variable list by a buffer variable hglassfacOld to save hour glass factors from previous turn
!       --> add if-statement to calculate the potential luminosity with betaSMin and check if this is above lumimax  (the levelling value of the luminosity), 
!           if yes, change value of betaS such that a luminosity of lumimax is produced in real luminosity calculation afterwards
!           
!
!   ***version 1.5
!       -->Adding transverse aperture cuts in betatron and momentum (subroutine collimation)
!       -->Adding transverse cut on initial distribution. This complicates emittance calculation!
!       -->Therefore, introducing several methods to calcualte the emittance once per turn for all subroutines (calcEmitWCut and fitemittance).
!       -->switch to turn on or off saving of all particle coordinates on turns when output is written
!       -->possibility to turn off longitudinal losses when particles leave bucket. Losses will then only take place when the physical aperture is hit due to dispersion or betatron action
!
!    ***version 1.4 (Tom Mertens)
!        -->Added Nagaitsev IBS method
!
!    ***version 1.3.1 (Tom Mertens)
!       -->Adapted to give initial emittances in x and y plane for both beams, hence the possibility to have non round initial beams
!        -> changes are made in the following modules :
!             - getinput
!             - addparticles
!             - keepTransvProf
!             - mountTransv
!
!    ***version 1.3
!      --> calculating radiation damping times and equilibrium emittances within code. 
!          User can now choose between calculation by code (lattic, with full rad. integrals, or approx) or writing manually the times in the input file
!      --> small restructuring of input file. path to file with tabulated values of Bane's g-function given as input parameter.
!      --> fixing bug with longitudinal phase space mismatch for large bunches: phase space density is now a function of the exact hamiltonian and not in small angle approximation
!
!    ***version 1.2
!      --adding Bane ibs method: fast and taking into account lattice. at collision gives similar result to modified Piwinski, but not good @ injection.
!        benchmark: LHC, 20e3 turns, 1e5 particles = 178 min. This is about 10% faster than Piwinski smooth.
!      --adding an additional output file ibs.out with the IBS rise times at each turn. 

!    ***version 1.1
!     --restructuring of input file, removing unused input parameters
!     --restructuring of output. using more logic naming conventions of output files, removing unused output. 
!     --writing output at the beginning of each turn instead of end
!     --cleaning up screen output

!    ***version 1.0.8
!     --changed order in input file: put the IP parameters last, 
!     for larger flexibility when automatically launching multiple stores. Can have multiple IPs in that are not used.


!    ***version 1.0.7
!     taken out calculation of IBS growth rate into separate routine. Different IBS methods possible:
!            --piwinski in smooth lattice approximation (same as original code)
!            --piwinski with lattice read in from madx twiss file (slow but more accurate)
!            --modified Piwinski with lattice
!            --interpolation from external table


!     ***version 1.0.6
!     in definition of q Piwinski's q, line 915 in ibslong:
!     -corrected missing factor sqrt(2) 
!     -corrected definition of d, so that the smallest beam size is taken
!     -improved screen output

!     ***version 1.0.5
!     cleaned up code. added comments and documentation in the beginning of file.
!     made notation consistent: x1,y1 are the transverse coords of bunch1, 
!     x2,y2 are transverse coords of bunch2 everywhere in code.
!     
!     ***version 1.0.4
!     in order not to introduce assymetries, the crossing angle in the the detailed routine 
!     collision1d is optionally in the hor. plane on every 2nd turn and in the vertical
!     plane on the other turn. Better would be to have a specific angle for each IP, 
!     but this would increase the CPU time even more.
!     removed multiple output of transv emittance
!     changed the ending of the allcoord* files to .out, so that all output files end with that

!     ***version 1.0.3
!     changed the way emittance is dumped to file
!     corrected bug in dumping the vertical beam profile to file-zeroing out avgliney
!     for debugging, dumping emittances from within the collision routine

!     ***version 1.0.2
!     exact expression for vertical coordinates and ds2/dz2 in collision1d 
!     using the collision probabilities to calculate the reduction factor

!     ***version 1.0.1
!     corrected bug in Moller factor in collision1d



      include 'collider_time_evolution_common.f90'
      integer :: np1,np2,np10,np20,nLostLum1,nLostLum2,nLostLumSum1,nLostLumSum2
      integer :: nLostMom1,nLostBeta1,nLostMomSum1,nLostBetaSum1,nLostMom2,nLostBeta2,nLostMomSum2,nLostBetaSum2
      integer :: iwrite,kcheck,nLostDebunch1,nLostDebunch2,nLostDebunchSum1,nLostDebunchSum2
      real :: fmix1,fmix2,lumi,redfac, fdum
      
      real :: x1(nb),px1(nb),t1(nb),pt1(nb),y1(nb),py1(nb),xpavg1(nb)  ! variables for beam1
      real :: ykeep1(7,nb),pnumber1,emixy1,tauhat1,avgline1(nb),avglinex1(nb),avgliney1(nb),rmsDelta1,rmsBunchLen1,&
              emix1,emiy1     
!added emix1, emiy1 to have unround beams. These emittances are starting conditions
      real :: ex1,ex2,ey1,ey2 ! these emittances are dynamically updated throughout the run

      real :: x2(nb),px2(nb),t2(nb),pt2(nb),y2(nb),py2(nb),xpavg2(nb)  ! variables for beam2
      real :: ykeep2(7,nb),pnumber2,emixy2,tauhat2,avgline2(nb),avglinex2(nb),avgliney2(nb),rmsDelta2,rmsBunchLen2,&
              emix2,emiy2
!added emix2, emiy2 to have unround beams

!**** MSCHAUMANN *******
! Add variables needed for Stochastic cooling
      real :: denlono1(nb),denlonoo1(nb),dpwakenotchx1(nb),dpwakenotchy1(nb),dpytot1(nb),dpxtot1(nb),coherecs1,coherecss1,&
              wakeyy1(nb),wakeyx1(nb),wakes1(nb),dpalmer1(nb),dpwakex1(nb),dpwakey1(nb),dpwakes1(nb)
      real :: denlono2(nb),denlonoo2(nb),dpwakenotchx2(nb),dpwakenotchy2(nb),dpytot2(nb),dpxtot2(nb),coherecs2,coherecss2,&
              wakeyy2(nb),wakeyx2(nb),wakes2(nb),dpalmer2(nb),dpwakex2(nb),dpwakey2(nb),dpwakes2(nb)


! just for testin cpu usage      
      real :: start, finish 

! dpxtot remembers horizontal pickup measurement
! dpytot remembers vertical pickup measurement
! dpalmer remembers longitudinal pickup measurement
! denlonoo = dpalmer from previous turn
! denlono = denlonoo from previous turn = dpalmer from two truns ago (only used if ntdelay = 1)

!************* MSCHAUMA ***************************
! Needed for beta* leveling
      real :: hglassfacOld(10)
      integer :: kkk
      do kkk=1, 10
       hglassfacOld(kkk) =1
      enddo
!*****************************

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!***** DUMMY Value for now **** will be given in input file later
!nresamp = 4000
!ntdelay = 0
!snrinv = 0.7e-6
!nturnon = 1
!harmbtf = 5.33e5
!ampbtf = 0.
!phibtf = 0.
!nextra = 0.
!writeStochCoolSwitch = 0
!fimpedx0 = 5e9
!fimpedx1 = 20e9
!fimpedy0 = 5e9
!fimpedy1 = 20e9
!fimpeds0 = 5e9
!fimpeds1 = 20e9 
!gainx0 = -4e6
!gainx1 = -4e6
!gainy0 = 4e6
!gainy1 = 4e6
!fracx = 0.675355
!fracy = 0.667456
!fracs = 0.66666
!coffsetx =0
!coffsety =0
!taupart = 8.E-12


      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '*************************************'
      write(6,*) '*      COLLIDER TIME EVOLUTION      *'
      WRITE(6,*) '*************************************'
      WRITE(6,*)
      write(6,*) 'processing input file ... '

!changed to get emix and emiy for both beams as input
      call getinput(pnumber1,pnumber2,emix1,emiy1,emix2,emiy2,tauhat1,tauhat2,rmsDelta1,rmsBunchLen1,rmsDelta2,rmsBunchLen2)

!************ MSCHAUMANN ***************
      if(stochCoolSwitch.eq.1) then
        call getimped(tlob1,wakeyy1,wakeyx1,wakes1) 
        call getimped(tlob2,wakeyy2,wakeyx2,wakes2)
      
        !zeror out memory arrays for stochastic cooling
        call zerooutarrays(denlono1,denlonoo1,dpwakenotchx1,dpwakenotchy1,dpytot1,dpxtot1,dpalmer1)
        call zerooutarrays(denlono2,denlonoo2,dpwakenotchx2,dpwakenotchy2,dpytot2,dpxtot2,dpalmer2)
     endif
!************************************

      call addparticles(np1,nMacro1,x1,px1,t1,pt1,y1,py1,np10,emix1,emiy1,tauhat1,rmsDelta1,rmsBunchLen1)      
      call addparticles(np2,nMacro2,x2,px2,t2,pt2,y2,py2,np20,emix2,emiy2,tauhat2,rmsDelta2,rmsBunchLen2)

! zero out sums
      call zerokeep(ykeep1,avgline1,avglinex1,avgliney1,np1)
      call zerokeep(ykeep2,avgline2,avglinex2,avgliney2,np2)
      nLostLum1=0
      nLostLumSum1=0
      nLostLum2=0
      nLostLumSum2=0
      nLostDebunch1=0
      nLostDebunch2=0
      nLostDebunchSum1=0
      nLostDebunchSum2=0
      nLostMom1=0
      nLostMomSum1=0
      nLostBeta1=0
      nLostBetaSum1=0
      nLostMom2=0
      nLostMomSum2=0
      nLostBeta2=0
      nLostBetaSum2=0


      write(6,*) 
      write(6,*) 'Starting main loop ...'
      write(6,*)
      write(6,*) '     turn      np1        np2     ex1(m)        ey1(m)    &
     &     sigt1(s)        FWHM_1(s)       ex2(m)         ey2(m)   &
     &     sigt2(s)        FWHM_2(s)     '
      write(6,*) '______________________________________________________&
     &__________________________________________________________________&
     &_______________________________________________'


!     calculate transverse emittances
      call getemittance(x1,y1,px1,py1,np1,x2,y2,px2,py2,np2,ex1,ex2,ey1,ey2)

! initial emittances
      ! write(6,*) 
      ! write(6,*) 'initial emittances: '
      ! write(6,*) 'ex1=',ex1
      ! write(6,*) 'ey1=',ey1
      ! write(6,*) 'ex2=',ex2
      ! write(6,*) 'ey2=',ey2
      ! write(6,*) 


! write first turn
      kturns=0
      iwrite=1
      if (writeMountSwitch.eq.1) then
	 call mountainr(np1,y1,py1,t1,pt1,x1,px1,1,np10,pnumber1,avgline1)
         call mountainr(np2,y2,py2,t2,pt2,x2,px2,2,np20,pnumber2,avgline2)
         call mountTransv(avglinex1,avgliney1,emix1,emiy1,1)
         call mountTransv(avglinex2,avgliney2,emix2,emiy2,2)
      endif
      call writeEmi(ex1,ey1,t1,pt1,ex2,ey2,t2,pt2,np1,np2,tauhat1,tauhat2)
      call writeNb(np1,np2,np10,np20,pnumber1,pnumber2,nLostLum1,nLostLum2,nLostLumSum1,nLostLumSum2,&
           nLostDebunch1,nLostDebunch2,nLostDebunchSum1,nLostDebunchSum2,&
           nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1,nLostBeta2,nLostBetaSum2,nLostMom2,nLostMomSum2)

      if (writeAllCoordSwitch.eq.1) then   !     write all starting conditions to file
         call writecoord(np1,x1,px1,y1,py1,t1,pt1,0,1)
         call writecoord(np2,x2,px2,y2,py2,t2,pt2,0,2)
      endif

      !     Main loop over all turns
!********************************************************************************
      do kturns = 0,nturns+nextra
!     logic for writes. Write output every nwrite turns
!------------------------------------------------------
       ! test cpu usage	
!	call cpu_time(start)
 
         kcheck = kturns/nwrite
         kcheck = kcheck*nwrite
         if(kcheck.eq.kturns)then
            iwrite = 1
         else
            iwrite = 0
         endif
!     always write the last turn
         if(kturns.eq.nturns+nextra)then
            write(6,*)' writing last turn'
            nwrite=-1
            iwrite=1
         endif

!     calculate transverse emittances
         call getemittance(x1,y1,px1,py1,np1,x2,y2,px2,py2,np2,ex1,ex2,ey1,ey2)
 
!     Step through different physics processes:
!------------------------------------------------
!
!********** MSCHAUMANN ************
	if(stochCoolSwitch.eq.1) then

! ***     Start at Kicker, fracx is fraction of turn between Pickup and Kicker
          if(fracy.lt.fracx)then
!           go in beam direction: Kicker to pickup distance = 1-fracx
            fmix1 = 1-fracx

            if(RFswitch.eq.1) then
              call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fmix1,nLostDebunch1,iwrite)
              call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fmix1,nLostDebunch2,iwrite)
            endif
            if(betatronSwitch.eq.1) then
              call sptransChrom(np1,x1,px1,y1,py1,pt1,fmix1)
              call sptransChrom(np2,x2,px2,y2,py2,pt2,fmix1)
            endif

!           use dpxtot for horizontal kick 
!           filling dpxtot at horizontal pickup (dpalmer, denlono, denlonoo are not touched)
            call setupres('tran',np1,x1,t1,tlob1,coffsetx,dpalmer1,denlono1,denlonoo1,dpxtot1)
            call setupres('tran',np2,x2,t2,tlob2,coffsetx,dpalmer2,denlono1,denlonoo2,dpxtot2)

!           the beam goes to the next device in the tunnel: vertical pickup
            fmix1 = fracx-fracy
            if(RFswitch.eq.1) then
              call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fmix1,nLostDebunch1,iwrite)
              call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fmix1,nLostDebunch2,iwrite)
            endif
            if(betatronSwitch.eq.1) then
              call sptransChrom(np1,x1,px1,y1,py1,pt1,fmix1)
              call sptransChrom(np2,x2,px2,y2,py2,pt2,fmix1)
            endif
!           filling dpytot at vertical pickup (dpalmer, denlono, denlonoo are not touched)
            call setupres('tran',np1,y1,t1,tlob1,coffsety,dpalmer1,denlono1,denlonoo1,dpytot1)
            call setupres('tran',np2,y2,t2,tlob2,coffsety,dpalmer2,denlono1,denlonoo2,dpytot2)

!           go to longitudinal pickup, which should always be the closest one to the kicker (from beam point of view)
!            fmix1 = fracx-2*fracs
            fmix1 = fracy-fracs
            if(RFswitch.eq.1) then 
              call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fmix1,nLostDebunch1,iwrite)
              call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fmix1,nLostDebunch2,iwrite)
            endif
            if(betatronSwitch.eq.1) then
              call sptransChrom(np1,x1,px1,y1,py1,pt1,fmix1)
              call sptransChrom(np2,x2,px2,y2,py2,pt2,fmix1)
            endif
!           filling dpalmer, denlono, denlonoo at longitudinal pickup (dpxtot is not touched)
            call setupres('long',np1,x1,t1,tlob1,coffsetx,dpalmer1,denlono1,denlonoo1,dpxtot1)
            call setupres('long',np2,x2,t2,tlob2,coffsetx,dpalmer2,denlono1,denlonoo2,dpxtot2)

!           finish the turn
            fmix1 = fracs
            if(RFswitch.eq.1) then
              call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fmix1,nLostDebunch1,iwrite)
              call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fmix1,nLostDebunch2,iwrite)
            endif
            if(betatronSwitch.eq.1) then
              call sptransChrom(np1,x1,px1,y1,py1,pt1,fmix1)
              call sptransChrom(np2,x2,px2,y2,py2,pt2,fmix1)
            endif

! ***      particles are at kicker again

! in case the vertical pickup comes before the horizontal one start at kicker and first go to vertical pickup
          else
            fmix1 = 1-fracy
            if(RFswitch.eq.1) then
              call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fmix1,nLostDebunch1,iwrite)
              call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fmix1,nLostDebunch2,iwrite)
            endif
            if(betatronSwitch.eq.1) then
              call sptransChrom(np1,x1,px1,y1,py1,pt1,fmix1)
              call sptransChrom(np2,x2,px2,y2,py2,pt2,fmix1)
            endif
            call setupres('tran',np1,y1,t1,tlob1,coffsety,dpalmer1,denlono1,denlonoo1,dpytot1)
            call setupres('tran',np2,y2,t2,tlob2,coffsety,dpalmer2,denlono1,denlonoo2,dpytot2)

            fmix1 = fracy-fracx
            if(RFswitch.eq.1) then
              call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fmix1,nLostDebunch1,iwrite)
              call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fmix1,nLostDebunch2,iwrite)
            endif
            if(betatronSwitch.eq.1) then
              call sptransChrom(np1,x1,px1,y1,py1,pt1,fmix1)
              call sptransChrom(np2,x2,px2,y2,py2,pt2,fmix1)
            endif
            call setupres('tran',np1,x1,t1,tlob1,coffsetx,dpalmer1,denlono1,denlonoo1,dpxtot1)
            call setupres('tran',np2,x2,t2,tlob2,coffsetx,dpalmer2,denlono1,denlonoo2,dpxtot2)

!            fmix1 = fracx-2*fracs
            fmix1 = fracx-fracs
            if(RFswitch.eq.1) then
              call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fmix1,nLostDebunch1,iwrite)
              call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fmix1,nLostDebunch2,iwrite)
            endif
            if(betatronSwitch.eq.1) then
              call sptransChrom(np1,x1,px1,y1,py1,pt1,fmix1)
              call sptransChrom(np2,x2,px2,y2,py2,pt2,fmix1)
            endif
            call setupres('long',np1,x1,t1,tlob1,coffsetx,dpalmer1,denlono1,denlonoo1,dpxtot1)
            call setupres('long',np2,x2,t2,tlob2,coffsetx,dpalmer2,denlono1,denlonoo2,dpxtot2)

            fmix1 = fracs
            if(RFswitch.eq.1) then
              call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fmix1,nLostDebunch1,iwrite)
              call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fmix1,nLostDebunch2,iwrite)
            endif
            if(betatronSwitch.eq.1) then
              call sptransChrom(np1,x1,px1,y1,py1,pt1,fmix1)
              call sptransChrom(np2,x2,px2,y2,py2,pt2,fmix1)
            endif

! ***     particles are at kicker 
          endif
!         endif for choise between x or y pickup comes first

! update routines are done  dpytot is for vertical kick
! dpwakenotchx is routine for memory of x
! cooling kicks only after rf is stable

! apply kick
          call resupdate(pnumber1,np1,np10,t1,pt1,px1,py1,emix1,emiy1,tlob1,dpxtot1,dpytot1,dpalmer1,denlono1,denlonoo1,&
                          coherecs1,coherecss1,wakeyy1,wakeyx1,wakes1,dpwakex1,dpwakey1,dpwakes1,iwrite)
          call resupdate(pnumber2,np2,np20,t2,pt2,px2,py2,emix2,emiy2,tlob2,dpxtot2,dpytot2,dpalmer2,denlono2,denlonoo2,&
                          coherecs2,coherecss2,wakeyy2,wakeyx2,wakes2,dpwakex2,dpwakey2,dpwakes2,iwrite)
! update coupling
          call sptransCoupling(np1,x1,px1,y1,px1)
          call sptransCoupling(np2,x2,px2,y2,px2)

!***********************


        else !for stochastic cooling switched off
!         synchrotron motion
          fmix1= 1
          if(RFswitch.eq.1) then
!             write(6,*) 'RF active'
             call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fmix1,nLostDebunch1,iwrite)
             call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fmix1,nLostDebunch2,iwrite)
          endif
         
!         betatron motion
          if(betatronSwitch.eq.1) then
!             write(6,*) 'betatron active'
             call sptransChrom(np1,x1,px1,y1,py1,pt1,fmix1)
             call sptransChrom(np2,x2,px2,y2,py2,pt2,fmix1)

!            ********* MSCHAUMANN **********
!            seperate the coupling between x and y to be able to apply it only once per turn
             call sptransCoupling(np1,x1,px1,y1,px1)
             call sptransCoupling(np2,x2,px2,y2,px2)
          endif

        endif !stochastic cooling switch


! all other processes have to be called ony once per turn, also in case of stochastic cooling

         call keepTransvProf(avglinex1,avgliney1,x1,y1,np1,emix1,emiy1)
         call keepTransvProf(avglinex2,avgliney2,x2,y2,np2,emix2,emiy2)

!     radiation damping and quantum excitation
         if(raddampSwitch.eq.1) then
!            write(6,*) 'raddamp active'
            call raddamp(np1,x1,px1,y1,py1,pt1)
            call raddamp(np2,x2,px2,y2,py2,pt2)
         endif

!     ibs. normalize equivalent time to beam 1
         if(IBSswitch.eq.1) then
!            write(6,*) 'IBS active'
           call ibslong(np1,y1,py1,t1,pt1,x1,px1,ex1,ey1,np10,pnumber1,pnumber1,avgline1,iwrite)
           call ibslong(np2,y2,py2,t2,pt2,x2,px2,ex2,ey2,np20,pnumber2,pnumber1,avgline2,iwrite)
        endif

!     collimation/aperture
        if(collimationSwitch.eq.1) then
           call collimation(np1,y1,py1,t1,pt1,x1,px1,np10,pnumber1,nLostMom1,nLostBeta1)
           call collimation(np2,y2,py2,t2,pt2,x2,px2,np20,pnumber2,nLostMom2,nLostBeta2)
        endif
       
!       collision between beam 1 and 2     
        if(collisionSwitch.eq.1) then
!          call writecoord(np1,x1,px1,y1,py1,t1,pt1,0,1)
!          call writecoord(np2,x2,px2,y2,py2,t2,pt2,0,2)
           if(collRoutine.eq.'6a') then
              call collision6a(np1,x1,px1,t1,pt1,y1,py1,np10,pnumber1,ex1,ey1,&
                               np2,x2,px2,t2,pt2,y2,py2,np20,pnumber2,ex2,ey2,&
                               nLostLum1,nLostLum2,iwrite,lumi,redfac,hglassfacOld)
!	    ***** MSCHAUMANN 19/06/2013 *****
!	    add an rfupdate right after the collison took place to get rid of high amplitude particles
!	    hope this helps with the strage behaviour of the bunch length for large initial values, help in Mikes programm
!            fdum = 0
!            call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fdum,nLostDebunch1,iwrite)
!            call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fdum,nLostDebunch2,iwrite)
           elseif(collRoutine.eq.'1d') then
              call collision1d(np1,x1,px1,t1,pt1,y1,py1,np10,pnumber1,&
     &             np2,x2,px2,t2,pt2,y2,py2,np20,pnumber2,nLostLum1,nLostLum2,iwrite,lumi,redfac,emixy1)
!	     **** MSCHAUMANN also for the second routine ******
!             fdum = 0
!             call rfupdate(np1,x1,px1,t1,pt1,y1,py1,xpavg1,fdum,nLostDebunch1,iwrite)
!             call rfupdate(np2,x2,px2,t2,pt2,y2,py2,xpavg2,fdum,nLostDebunch2,iwrite)
           else
              write(6,*) 'unkonwn collision routine'
              stop
           endif
!          call writecoord(np1,x1,px1,y1,py1,t1,pt1,1,1)
!          call writecoord(np2,x2,px2,y2,py2,t2,pt2,1,2)
        endif
              
        call keepturn(np1,y1,py1,ykeep1)
        call keepturn(np2,y2,py2,ykeep2)

!     first turn already written
        if (kturns.eq.0) then 
           iwrite=0
           call writeLumi(lumi,redfac)
        endif
!     write output if desired
        if (iwrite.eq.1) then
	   if(writeStochCoolSwitch.eq.1) then
!            call writemomentsshort(np1,y1,py1,t1,pt1,x1,px1,1)
!            call writemomentsshort(np2,y2,py2,t2,pt2,x2,px2,2)
             call writemomentsshort(np1,y1,py1,t1,pt1,x1,px1,1,coherecs1,coherecss1,tlob1,pnumber1,np10,&
                                    dpalmer1,dpxtot1,dpytot1,dpwakex1,dpwakey1,dpwakes1,ykeep1,xpavg1)
             call writemomentsshort(np2,y2,py2,t2,pt2,x2,px2,2,coherecs2,coherecss2,tlob2,pnumber2,np20,&
                                    dpalmer2,dpxtot2,dpytot2,dpwakex2,dpwakey2,dpwakes2,ykeep2,xpavg2)
           endif
           write(6,*)
           if (writeAllCoordSwitch.eq.1) then   !     write all starting conditions to file
              call writecoord(np1,x1,px1,y1,py1,t1,pt1,kturns,1) 
              call writecoord(np2,x2,px2,y2,py2,t2,pt2,kturns,2) 
           endif
           if (writeMountSwitch.eq.1) then
              call mountainr(np1,y1,py1,t1,pt1,x1,px1,1,np10,pnumber1,avgline1)
              call mountainr(np2,y2,py2,t2,pt2,x2,px2,2,np20,pnumber2,avgline2)
              call mountTransv(avglinex1,avgliney1,emix1,emiy1,1)
              call mountTransv(avglinex2,avgliney2,emix2,emiy2,2)
           endif
           call writeEmi(ex1,ey1,t1,pt1,ex2,ey2,t2,pt2,np1,np2,tauhat1,tauhat2)
           call writeNb(np1,np2,np10,np20,pnumber1,pnumber2,nLostLum1,nLostLum2,nLostLumSum1,nLostLumSum2,&
                nLostDebunch1,nLostDebunch2,nLostDebunchSum1,nLostDebunchSum2,&
           nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1,nLostBeta2,nLostBetaSum2,nLostMom2,nLostMomSum2)
           if(collisionSwitch.eq.1) call writeLumi(lumi,redfac)


        endif
! test cpu usage	
!	call cpu_time(finish)
!        write(6,*) ' cpu time used = ', finish-start

     enddo
     write(6,*) '_____________________________________________________________________________&
          &__________________________________________'
     write(6,*)
     write(6,*) 'COLLIDER TIME EVOLUTION: run finished normally'
     write(6,*) '**********************************************'
     stop
   end program collider_time_evolution_cool
      

!*********************************************************************************
!*********************************************************************************

      subroutine getinput(pnumber1,pnumber2,emix1,emiy1,emix2,emiy2,tauhat1,tauhat2,rmsDelta1,rmsBunchLen1,rmsDelta2,rmsBunchLen2)
      include 'collider_time_evolution_common.f90'
      real :: pnumber1,pnumber2,emix1,emiy1,emix2,emiy2,tauhat1,tauhat2,&
           eta,tunesexact,rcoeff,tunes,rmsDelta1,rmsBunchLen1,rmsDelta2,&
           rmsBunchLen2,v1,s,dgamma_hat,rIon,CalphaE3C,I2,I4x,I4y,I3,&
           sigEoE0,Cq,alfax,alfay,gamax,gamay,Dx,Dy,DxP,DyP,&
           Hx,Hy,I5x,I5y,Jx,Jy,ki,kis,rhoi

      integer :: k,nfft,m,i,j,IOstatus
      character(4) :: col1,col2,col3,col4,col5,col6,col7,col8,col9,&
          col10,col11,col12,col13,col14,col15
      character(150) :: twissFile,ibsIntFile,gTabFile
! 
      open(unit=10,file='collider_time_evolution_cool.in',status='unknown')
      read(10,*) col1
      
! switches to turn on various processes
      read(10,*) RFswitch       !(set to 1 to activate synchrotron motion)
      read(10,*) betatronSwitch !(set to 1 to activate betatron motion)
      read(10,*) raddampSwitch  !(set to 1 to activate radiation damping and quantum excitation)
      read(10,*) IBSswitch      !(set to 1 to activate ibs)
      read(10,*) collimationSwitch !(set to 1 to activate losses on aperture cuts)
      read(10,*) collisionSwitch !(set to 1 to activate luminosity)
      read(10,*) levellingSwitch !(set to 1 to activate levelling the luminosity)
      read(10,*) stochCoolSwitch !(set to 1 to activate stochastic cooling)
      read(10,*) col1

! **** MSCHAUMANN ****
! to re-implement stochstic cooling the input parameter timeRatio is replaced by eqTime (equivalent time [hour], real time in the machine).
! If stochCoolSwitch = 1, the nturns and nMacro is overwritten 


! nturns=number of turns 
! nMacro=number of macro particles per bunch
! timeRatio=real machine turns per simulation turn
! nwrite=number of turns between writes of output
! writeAllCoordSwitch = set to 1 to activate a dump of all particle coorinates on every turn where output is written
! writeMountSwitch = set to 1 to activate mountain range output of distribution shape in transv and long planes
! iseed=random number seed

      read(10,*)nturns,nMacro1,eqTime,nwrite,writeAllCoordSwitch,writeMountSwitch,iseed
!      read(10,*)nturns,nMacro,timeRatio,nwrite,writeAllCoordSwitch,writeMountSwitch,iseed
!      write(6,*) '    nturns=',nturns,', # macro part.=',nMacro

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
      read(10,*)tunex,tuney,chromx,chromy,dqmin

! smooth approximation of lattice functions
      betax=circ/(2*pi*tunex)
      betay=circ/(2*pi*tuney)

! aatom=mass number of ion species
! qatom=charge number of ion species
! [-thib/2:thib/2] is the interval used for losses. particles with abs(t)>thib/2 are considered lost. 
!                 For single RF, thib=RF period. If thib=0.0, an infinite boundary is assumed and no losses occur longitudinally.
      read(10,*) aatom,qatom,thib
      read(10,89) twissFile ! optics twissfile - used for ibs with ibsMethod='piwLattice' or 'modPiwLat' and for rad. damping with radMethod='lattic'
      read(10,*) col1

! emixy1= geom. 1 sigma starting emittance bunch 1 (metres)
! emix1, emiy1 = geom. 1 sigma starting emittance bunch 1 (metres)
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
!        This creates an exactly matched phase-space which does not filament
      read(10,*) longCoordMethod

! rmsBunchLen1=rms of bunch length for beam 1 (in metres, used only with longCoordMethod=2 and longCoordMethod=3 and Stochastic cooling)
      read(10,*)  rmsBunchLen1,rmsBunchLen2
      tlob1 = 4*rmsBunchLen1/clight
      tlob2 = 4*rmsBunchLen2/clight

! rmsDelta1=rms of deltaP/P0 for beam1 (only used with longCoordMethod=2)
      read(10,*) rmsDelta1,rmsDelta2

! tauhat1=half bucket length (seconds) for beam 1 of the RF system with shortest wavelength
! tauhat2=half bucket length (seconds) for beam 2 of the RF system with shortest wavelength
! bunchLenPrecis=precision in sampling of bunch length [upper bound on abs(sampled/wanted -1)]. ex 0.01 gives 1% accuracy
! note: for single RF, tauhat=thib/2 
! power=bunch shape parameter, only used with smoke ring
! alint = smoke ring parm, only used with smoke ring
      read(10,*) tauhat1,tauhat2,bunchLenPrecis,power,alint

! files with long. starting conditions, used only if longCoordMethod=1
 89   format(A120)
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

!************ MSCHAUMANN ********************
! Stochastic Cooling Parameters
! fimped = lower and upper frequencies 
! gain = lower and upper impedances
! fraci = the fraction of a turn delay for cooling in i=x,y,s(=longitudinal)
! fracs smaller than fracx and fracy
      read(10,*) fimpedx0, fimpedx1, gainx0, gainx1, fracx
      read(10,*) fimpedy0, fimpedy1, gainy0, gainy1, fracy
      read(10,*) fimpeds0, fimpeds1, gains0, gains1, fracs

! nextra = extra turns with cooling off
! ntdelay = 0 for 1 turn delay,  1 for cascaded 1 turn delays
! nturnon=number of turns to turn on nonlinear RF and longitudinal wakes
      read(10,*) nextra, ntdelay, nturnon

! amplitude, phase and carrier harmonic for btf kick
      read(10,*) ampbtf, phibtf, harmbtf

! nresamp = number of samples used for ffts
! coffset = closed orbit offset in x and y
! snrinv = ratio of rms noise to rms signal (inverse signal-to-noise)
! taupart = equivalent length of a macro-particle
! vmaxreal = max voltage allowed on longitudinal kicker, requested voltages above the given one are set to vmaxreal (cutoff/saturation value)
      read(10,*) nresamp, coffsetx, coffsety, snrinv,taupart
      read(10,*) vmaxreal
      read(10,*) writeStochCoolSwitch
      read(10,*) col1
!***********************

! refEmxy = reference geometric transverse emittance used for the collimation/aperture cuts in betatron space (metres)
! cutoffAmpl = amplitude in refeence sigma =sqrt(beta*refEmxy) at which initial transverse distr is cut (dimensionless)
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

      xcut=cutoffAmpl*sqrt(betax*refEmxy) ! cutoff in metres of the transverse distributions
      ycut=cutoffAmpl*sqrt(betay*refEmxy)


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
! lumimax = maximum luminosity for levelling
! levIPSwitch = defines if the IP is levelled (set to 1) or not
! ********* MSCHAUMA **************
! add lumimax and levIPSwitch to input parameters for beta* leveling per IP
      do m=1,nIPs               !loop over all IPs, read beta*, crossing angle/2 and multiplicity
         read(10,*) betaS(m),theta(m),IPmult(m),lumimax(m),levIPSwitch(m)
         betaSMin(m) = betaS(m)
      enddo

! calculate trev, beta etc
      beta= sqrt(1-1/gamma0**2)
      vrev = clight*beta
      trev = circ/vrev
      trf = trev/nharm
      frf = 1/trf
      omegarf = 2*pi/trf
      eta = 1/gammat**2 - 1/gamma0**2

! calculate equivalent simulation time
! **** MSCHAUMA **** eqTime is now an input parameter out of which timeRatio is calculated
! only overwrite nturns if stochastic cooling is switched on

	if (stochCoolSwitch.eq.1) then 
           nturns = NINT(eqTime*nMacro1/pnumber1/trev*3600)
           nMacro2 = NINT(nturns*trev*pnumber2/eqTime/3600)
           write(6,*) '      Stochastic Cooling is :         ON'
           write(6,*) '  ====>>>> nturns have been overwriten! <<<<==== '
           write(6,*) '  ====>>>>      nMacro1 != nMacro2      <<<<==== '
           write(6,*)
           write(6,*) ' lower and upper frequencies for x cooling: ', fimpedx0, ',  ' , fimpedx1
           write(6,*) ' lower and upper frequencies for y cooling: ', fimpedy0, ',  ' , fimpedy1
           write(6,*) ' lower and upper frequencies for s cooling: ', fimpeds0, ',  ' , fimpeds1
           write(6,*)
           write(6,*) ' transverse delays in oscillation:  ', fracx*tunex, ',  ' , fracy*tuney
           write(6,*)
        else
           nMacro2 = nMacro1
           write(6,*) '      Stochastic Cooling is :         OFF'
           write(6,*)
        endif
	
	timeRatio = FLOAT(NINT(eqTime/trev/nturns*3600))

      write(6,*) '    nturns=',nturns,', # macro part. B1=',nMacro1,', # macro part. B2=',nMacro2
      write(6,*)
      write(6,*) '    timeRatio=',timeRatio

! ****      eqTime=nturns*trev*timeRatio/3600.  
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

!!*********** MSCHAUMA ****************************
      write(6,*)
      if (levellingSwitch.eq.1) then 
         write(6,*) '    Levelling is:       ON'
         do m=1,nIPs  
            write(6,*) '       --> IP ', m, '  is levelled: ' , levIPSwitch(m), ' levelling value for luminosity:',lumimax(m)  
         enddo

      else
         write(6,*) '    Levelling is:       OFF'
      endif
!**********************************************


      close(10)
      write(6,*)
! makesure nresamp is a power of 2
      nfft=2
      do k=1,1000
         nfft = 2*nfft
         if(nfft.ge.nresamp)go to 22
      enddo
 22   continue
      write(6,*)'nresamp,nfft',nresamp,nfft
      nresamp=nfft
      
! calculate trev, beta etc
!      beta= sqrt(1-1/gamma0**2)
!      vrev = 2.998e8*beta
!      trev = circ/vrev
!      trf = trev/nharm
!      frf = 1/trf
!      omegarf = 2*pi/trf
!      eta = 1/gammat**2 - 1/gamma0**2
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
      call system('rm csmon*.out')
      call system('rm fort.*')
      call system('rm res_den.out')
      call system('rm z.out')
      call system('rm ws.out')
      call system('rm kickerVoltage*.out')

!      open(unit=32,file='tran.xvz',status='unknown')
!      open(unit=33,file='csmon.out',status='unknown')
!      open(unit=52,file='nLostBurnoff.out',status='unknown')

!      open(unit=34,file='tran.cs',status='unknown')
!      open(unit=44,file='tran.sig',status='unknown')

!     call getimped  !! used for cooling - call outside this function
! get fiducial normalizations
!      write(6,*)' dpwake for 1 meter offset with step wake = 1'
!      write(6,*)(circ/(2*pi*tuney))*1.6e-19*qatom*pnumber1/
!     :          (beta*gamma0*938.e6*aatom/qatom)

! normalization for longitudinal profiles
      nwrite0 = nwrite
! equivalent time in machine of simulation, rescaling with number of sim turns per real turns



!     read mad twiss file with lattice (for some methods of ibs and radiation damping only). Assume fixed format of tfs table for simplicity
      if (((raddampSwitch.eq.1).and.(radMethod.eq."lattic")).or.&
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
         read(99,*) col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12,col13,col14,col15 ! read column headers
            
!        check that we have the right columns in twiss file
         if((col2=="NAME").and.(col3=="S   ").and.(col4=="L   ").and.(col5=="BETX").and.(col6=="BETY")&
              .and.(col7=="ALFX").and.(col8=="ALFY").and.(col9=="DX  ").and.(col10=="DPX ")&
              .and.(col11=="DY  ").and.(col12=="DPY ").and.(col13=="ANGL ").and.(col14=="K1L ").and.(col15=="K1S ")) then
         else
            write(6,*) "Error in format of twiss file ", trim(twissFile)
            write(6,*) "Columns must be NAME|S   |L   |BETX|BETY|ALFX|ALFY|DX  |DPX |DY  |DPY |ANGLE|K1L|K1S"
            stop
         endif
         read(99,*) col1        ! read one more line to get to the actual table
         
!     read in lattice
         nElem=1                ! count the number of elements in lattice
         IOstatus=0             ! flag to check when end of file is reached
         
         do while(IOstatus.eq.0) ! read elements until end of file
            read(99,*,iostat=IOstatus) &
                 col1,S,leng(nElem),betx(nElem),bety(nElem),alfx(nElem),alfy(nElem),dispx(nElem),dispxP(nElem),dispy(nElem),&
                 dispyP(nElem),bendAng(nElem),K1L(nElem),K1S(nElem)
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
            Hx=(betax*DxP**2+2*alfax*Dx*DxP+gamax*Dx**2)
            Hy=(betay*DyP**2+2*alfay*Dy*DyP+gamay*Dy**2)

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
!                  I4x=I4x+(dispx(i)/rhoi**3*(1.+2.*rhoi**2*ki)+2.*dispx(i)*ki/rhoi)*leng(i)  ! corrected to equations in Handbook, Second Edition, p220
!                  I4y=I4y+(dispy(i)/rhoi**3*(1.-2.*rhoi**2*ki)-2.*dispy(i)*ki/rhoi)*leng(i)
                  I4x=I4x+(dispx(i)/rhoi**3*(1.+2.*rhoi**2*ki)+2.*dispy(i)*K1S(i)/rhoi)*leng(i)
                  I4y=0  ! in case of no vertical bending
                  gamax=(1.+alfx(i)**2)/betx(i)
                  gamay=(1.+alfy(i)**2)/bety(i)
                  Hx=betx(i)*dispxP(i)**2 + 2*alfx(i)*dispx(i)*dispxP(i)+gamax*dispx(i)**2
                  Hy=bety(i)*dispyP(i)**2 + 2*alfy(i)*dispy(i)*dispyP(i)+gamay*dispy(i)**2
!                  I5x=I5x+Hx*2*pi/rho0**2*leng(i)
!                  I5y=I5y+Hy*2*pi/rho0**2*leng(i)
                  I5x=I5x+Hx*1./rhoi**3*leng(i)
                  I5y=I5y+Hy*1./rhoi**3*leng(i)
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
         sigEoE0=sqrt(Cq*gamma0**2*I3/(2*I2+I4x+I4y))  ! eq 18
         siglong=gamma0*sigEoE0 ! sigma of pt
         Jx=1.-I4x/I2
         Jy=1.-I4y/I2
         sigperp=Cq*gamma0**2*I5x/(Jx*I2)   !  eq 22, transverse equilibrium emittance 
         write(6,*) '          calculating radiation damping times and equilibrium beam sizes:'
         write(6,*) '          Tx=',tradperp/3600,'h,  Ts=',tradlong/3600,'h,  sigPt=',siglong,'ex=',sigperp,'m'
      endif

      write(6,*)
      write(6,*)'... collider_time_evolution_cool.in read'
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
!      write(80,*) 'sim.turn     t(hours)         ex1(m)           ey1(m)&
!     &           el1(eV/s/charge) sig_T1           sig_dP/P_2      ex2(m&
!     &)           ey2(m)           el2(eV/s/charge) sig_T1           sig&
!     &_dP/P_2'

      write(80,*) 'sim.turn     t(hours)         ex1(m)           ey1(m)&
     &           el1(eV/s/charge) sig_T1(s)          FWHM_1(s)      sig_&
     &dP/P_1         ex2(m)           ey2(m)         el2(eV/s/charge) &
     &sig_T2(s)       FWHM_2(s)        sig_dP/P_2'
      write(80,*) '-----------------------------------------------------&
     &-----------------------------------------------------------------'
      close(80)

      open(unit=80,file='ibs.out')
      write(80,*) '    sim.turn     t(hours)      Tp(hours)       Tx(hou&
     &rs)        Ty(hours)'
      write(80,*) '-----------------------------------------------------&
     &-----------------------------------------------------------------'
      close(80)


      end
      

!*********************************************************************************
!*********************************************************************************

      REAL FUNCTION RAN3(IDUM)
!         IMPLICIT REAL*4(M)
!         PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      common/ran3stuff/ MA(55),inext,inextp
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      RETURN
      END

      
!*********************************************************************************
!*********************************************************************************

      subroutine sampleLongMatched(np,nMacro,t,pt,ham1sig,sigs,tauhat,ptmax,v00,omega0,hammax)
      include 'collider_time_evolution_common.f90'

      real :: t(nb),pt(nb),tk,ptk,ham1sig,sigs,tauhat,p1,p2,ptmax,prob,v00,test,omega0,ran3,ham,hammax
      integer :: np, nMacro
      
      sigs=0.

      do np=1,nMacro
55       continue
         tk = tauhat*(2*ran3(iseed)-1)
         ptk = ptmax*(2*ran3(iseed)-1)
         
         p1 = nharm*tk*omega0   !
         p2 = nharm2*tk*omega0

!     exact hamiltonian:
         ham=0.5*ptk*ptk*tcoeff+(cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)

!     small angle approximation:
!$$$  ham=0.5*ptk*ptk*tcoeff - p1**2/2*vrf/(nharm*v00*omega0) - 
!$$$  &         p2**2/2*vrf2/(nharm2*v00*omega0)

         
         if(ham.gt.hammax)go to 55 ! restart sampling if outside bucket
         
         prob=exp(-ham/ham1sig)
         test = ran3(iseed)
         if(prob.lt.test) go to 55
         sigs=sigs+tk**2
!         write(88,*) tk,ptk,ham
         t(np)=tk
         pt(np)=ptk
      enddo
      sigs=clight*sqrt(sigs/nMacro)
      return
      end

!*********************************************************************************
!*********************************************************************************


      subroutine addparticles(np,nMacro,x,px,t,pt,y,py,np0,emix,emiy,tauhat,rmsDelta,rmsBunchLen)
      include 'collider_time_evolution_common.f90'
      external ran3
      real :: x(nb),px(nb),t(nb),pt(nb),y(nb),py(nb) ! local beam variables
      real :: emix,emiy,tauhat,dum1,dum2,dum3,dum4,rmsDelta,rmsBunchLen,testemit

      real :: sigv,v00,omega0,ham,hammax,hamzero,amp,facc,p1,p2,pcoeff2,&
           phik,pcoeff3,phizero,phik0,prob,ptk,ptmax,radharm2,r1,r2,vk,&
           vzero,tk,ran3,test,ampx,ampy,ampt,ampPt,&
           epsl,sigs0,sigs1,ham1sig0,ham1sig1,fPrime
      integer :: i
      logical readfile1
      save readfile1
      data readfile1 /.true./

      integer :: np,np0,k,j,nMacro

! first get the zero of the voltage
      phizero=real(nharm)/real(nharm2)*3*pi/4
      vzero=abs(vrf)
      do k=-1000,1000
         phik=real(nharm)/real(nharm2)*(pi+k*pi/2000)
         vk = vrf*sin(phik)+vrf2*sin(nharm2*phik/nharm)
!         write(68,*)k,phik,vk
         if(abs(vk).lt.vzero)then
            vzero = abs(vk)
            phizero=phik
         endif
      enddo
      phik0=phizero
      do k=-1000,1000
         phik = phik0+ k/100000.
         vk = vrf*sin(phik)+vrf2*sin(nharm2*phik/nharm)
!         write(68,*)k,phik,vk
         if(abs(vk).lt.vzero)then
            vzero = abs(vk)
            phizero=phik
         endif
      enddo



!      write(6,*)'vzero,phizero',vzero,phizero

       np=0
       pcoeff2=pcoeff/omegarf
       pcoeff3 = pcoeff2/omegarf
       radharm2 = vrf2/vrf
! voltage for unit change in gamma
       v00 = aatom*938.e6/qatom
       sigv=0
       omega0 = 2*pi/trev
       p1 = nharm*tauhat*omega0
       p2 = nharm2*tauhat*omega0
  !     write(6,*) 'p1,p2',p1,p2
       hammax = (cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0) 

       p2=phizero*nharm2/nharm
       hamzero = (cos(phizero)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)
 
!--------------------------------------------
 !      write(6,*)' hammax = ',hammax
 !      write(6,*) 'ptmax=',ptmax
 !      write(6,*) 'ham0=',hamzero
 !      write(6,*) 'nharm=',nharm2
 !      write(6,*) 'gamma0=',gamma0
 !      write(6,*) 'v00=',v00
 !      write(6,*) 'omega0=',omega0
 !      write(6,*) 'vrf=',vrf
 !      write(6,*) 'tcoeff=',tcoeff
 !      write(6,*) 'rmsBunchLength=',rmsBunchLen
 !--------------------------------------------

       ptmax = sqrt(2*hammax/tcoeff)

! get longitudinal coordinates
! if switch is set to
! 0: parabolic distribution with smoke ring 
! 1: read from file
! 2: gaussian
       write(6,*) 'generating longitudinal coordinates'
       if (longCoordMethod==0) then ! smoke ring dist.
          do np = 1,nMacro
             
 43          continue
             tk = tauhat*(2*ran3(iseed)-1)
             ptk = ptmax*(2*ran3(iseed)-1)
             
             p1 = nharm*tk*omega0
             p2 = nharm2*tk*omega0
             ham=0.5*ptk*ptk*tcoeff+(cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)
             if(ham.gt.hammax)go to 43
             prob = 1-(ham/hammax)
             prob = prob**power
             if((abs(p1).gt.phizero).and.(ham.le.hamzero))then
                prob = prob*(ham/hamzero)**alint
             endif
          
             test = ran3(iseed)
             if(prob.lt.test)go to 43
             pt(np) = ptk
             t(np)=tk
          enddo
       elseif (longCoordMethod==1) then ! read from file
          if (readfile1) then ! logical variable to keep track of beam 1 and beam 2
             open(unit=10,file=trim(filebeam1))
             readfile1=.false.
             write(6,*) 'reading starting coordinates for beam1 from ', filebeam1
          else
             open(unit=10,file=trim(filebeam2))
             write(6,*) 'reading starting coordinates for beam2 from ', filebeam2
          endif
		write(6,*) 'Enter Loop over particles'
          do np = 1,nMacro
             read(10,*) dum1,dum2,dum3,dum4,t(np),pt(np)
		!write(6,*) np ,'   ', t(np), '    ' , pt(np)
          enddo

       elseif (longCoordMethod==2) then ! generate bi-gaussian (possibly unmatched, can only be matched in small angle approximation)
56          continue
             ampPt=gamma0*rmsDelta
             ampt=rmsBunchLen/clight
!             write(6,*) ampt,ampPt
!     r1 and r2 are uniform on (-1,1)
             r1 = 2*ran3(iseed)-1
             r2 = 2*ran3(iseed)-1
             amp = r1*r1+r2*r2
             if(amp.ge.1) go to 56
             facc = sqrt(-2.*log(amp)/amp)
             tk = ampt*r1*facc
             ptk = ampPt*r2*facc
             p1 = nharm*tk*omega0
             p2 = nharm2*tk*omega0
             ham=0.5*ptk*ptk*tcoeff+(cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)
             if(ham.gt.hammax)go to 56 ! restart sampling if outside bucket
             if (abs(tk).ge.tauhat) go to 56 ! restart sampling if t is too large
             t(np)=tk
             pt(np)=ptk

       elseif (longCoordMethod==3) then ! generate matched 'pseudo-Gaussian' phase space
!         Assuming distribution function for hamiltonian rho(h)=1/H*exp(-h/H), so rho(pt,t)~exp(-C1*pt^2-C2*(1-cos(C3*t)). at 2nd order in t this is a bi-gaussian
!         average H can not be calculated analytically for the exact hamiltonian, only for small angle approximation.
!         Use first this value to calculate bunch length, then use Newton's method to calculate new H and iterate until convergence

!         get average H in small osc. approx. see derivation in check_phase_space_matching.nb
          ham1sig0=-nharm2*omega0*vrf2/(clight**2*v00)*rmsBunchLen**2
          call sampleLongMatched(np,nMacro,t,pt,ham1sig0,sigs0,tauhat,ptmax,v00,omega0,hammax) 
          ham1sig1=ham1sig0*(rmsBunchLen/sigs0)**2
          call sampleLongMatched(np,nMacro,t,pt,ham1sig1,sigs1,tauhat,ptmax,v00,omega0,hammax)
!         use Newton's method and iterate until we get the right bunch length
!         consider sigs=f(h) as function of ham1sig. we want to find zero of f(h)-sRms
!         iterate h1=h0 - (f(h0)-sRms)/f'(h0)
!         approximate f'=(f(h1)-f(h0))/(h1-h0)
          do while (abs(sigs1/rmsBunchLen-1.).ge.bunchLenPrecis)
             write(6,*) '       iterating: bunch length = ',sigs1,' m.'
             fPrime=(sigs1-sigs0)/(ham1sig1-ham1sig0)
             ham1sig0=ham1sig1
             sigs0=sigs1
             ham1sig1=ham1sig0-(sigs0-rmsBunchLen)/fPrime
             call sampleLongMatched(np,nMacro,t,pt,ham1sig1,sigs1,tauhat,ptmax,v00,omega0,hammax)
!             write(*,*) ham1sig0, ham1sig1, sigs0, sigs1
          enddo
          write(6,*) '       iterating: bunch length = ',sigs1,' m.'
       else
          write(6,*) 'unknown method for long. coord.'
       endif
       sigv = sigv + pt(np)**2


! transverse coordinates:
! determine transverse RMS amplitude (assuming same emi. in x and y)
       ampx = sqrt(betax*emix)
       ampy = sqrt(betay*emiy)
!  **** MSCHAUMA XXX **** this has to be changed back to ampy = sqrt(betay*emiy) just for this test to compare with Mike's program
!       ampy = ampx 	

     
       write(6,*)
       write(6,*) 'generating transverse coordinates'

       

       do np = 1,nMacro
44        continue          
! r1 and r2 are uniform on (-1,1)
          r1 = 2*ran3(iseed)-1
          r2 = 2*ran3(iseed)-1
          amp = r1*r1+r2*r2
          if((amp.ge.1).or.(amp.lt. 3.e-6))go to 44
          facc = sqrt(-2.*log(amp)/amp)
! inject with no angle error
!        x(np) = xinject + ampx*r1*facc
! transverse kick will be given on turn nturnon
          x(np) = ampx*r1*facc
! px has same amplitude as x ie px = beta_L*x'
          px(np) = ampx*r2*facc
! reject if x>initial cut off (given in sigmas with the reference emittance)
          if (sqrt(x(np)**2+px(np)**2).ge.xcut) go to 44
!          if ((x(np)**2+px(np)**2)/(2*betax*refEmxy).ge.cutoffAmpl) go to 44

! other transverse variable
45        continue
! r1 and r2 are uniform on (-1,1)
          r1 = 2*ran3(iseed)-1
          r2 = 2*ran3(iseed)-1
          amp = r1*r1+r2*r2
          if((amp.ge.1).or.(amp.lt. 3.e-6))go to 45
          facc = sqrt(-2.*log(amp)/amp)
          y(np)= ampy*r1*facc
          py(np)= ampy*r2*facc
! reject if x>initial cut off (given in sigmas with the reference emittance)
          if (sqrt(y(np)**2+py(np)**2).ge.ycut) go to 45
!          if ((y(np)**2+py(np)**2)/(2*betay*refEmxy).ge.cutoffAmpl) go to 45
! reject if skew amplitude > cutoff
!          if ((y(np)**2+py(np)**2)/(2*betay*refEmxy)+(x(np)**2+px(np)**2)/(2*betax*refEmxy).ge.cutoffAmpl) go to 44
       enddo
       close(10)
       sigv = sqrt(sigv/np)
!       write(6,*)' rms of pt = ',sigv
! initial number of macro-particles for normalization of current etc
       np=np-1
       np0 = np

       call fitemittance(x,px,np,betax,cutoffAmpl,testemit)
       write(6,*) 'emitx from fit = ',testemit
       call fitemittance(y,py,np,betay,cutoffAmpl,testemit)
       write(6,*) 'emity from fit = ',testemit

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


!*********************************************************************************
!*********************************************************************************

     subroutine emitStDev(x,px,np,bx,emit)
       include 'collider_time_evolution_common.f90'
       real :: x(nb),px(nb),bx,emit
       integer :: i,np
  
       emit=0.0
       do i=1,np
          emit=emit+x(i)**2+px(i)**2
       enddo
       emit=emit/bx/np/2.

       return
     end subroutine emitStDev

!*********************************************************************************
!*********************************************************************************

subroutine fitEmittance(x,px,np,bx,Jmax,emit)
  include 'collider_time_evolution_common.f90'
  real :: x(nb),px(nb),binVals(1000),Jx,Jmax,bx,sigma1,sigma0,getNewtonSigma,bin1,emit,binwidth
  integer :: i,np,bin,maxbin
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


!*********************************************************************************
!*********************************************************************************

real function getNewtonSigma(sigmaOld,binVals,Jmax,binwidth)
  include 'collider_time_evolution_common.f90'
! for fitting sigma in the function g(x)=binVals(1)*Exp[-x/sigma] with Newton's method to the binned betatron actions, in order to calculate the emittance with a least square fit
! The sum of the squares is f(sigma)=sum_i[ (g(xi)-binCount(i))^2 ]. Find sigma that minimizes f
! a minimum is found for f'(sigma)=0
! with newton's method, we thus iterate sigma_(n+1) = sigma_n - f'(sigma_n)/f''(sigma_n)
! the fit only takes bins up to Jmax into account. Therefore we will get the correct emittance even if the tail of the Gaussian is cut.  
  integer :: maxBin,i
  real :: sigmaOld,binVals(500),Jmax,f,fprime,fbis,A,AExpXiSig,xi,yi,binwidth
  
  maxbin=nint(Jmax/binwidth) 
  fprime=0.
  fbis=0.
  A=binVals(1) ! amplitude of f
  do i=1,maxBin ! sum over all bins to calculate f'(sigmaOld) and f''(sigmaOld)
     xi = real(i)*binwidth-binwidth/2. ! get the center of the bin
     yi = binVals(i)
     AExpXiSig = A*Exp(-xi/sigmaOld)
     fprime=fprime + 2*AExpXiSig*xi*( AExpXiSig-yi )/sigmaOld**2
     fbis=fbis + (2*A*xi*AExpXiSig**2*(-(1/AExpXiSig*yi*(xi - 2*sigmaOld)) + 2*A*(xi - sigmaOld)))/sigmaOld**4
  enddo
!  write(6,*) "fi=",fprime,"fii=",fbis,"sigmanew=", sigmaOld - fprime/fbis
  getNewtonSigma = sigmaOld - fprime/fbis
  return
end function getNewtonSigma



!*********************************************************************************
!*********************************************************************************

subroutine calcEmitWCut(x,px,np,cutx,bx,emit)
  include 'collider_time_evolution_common.f90'
! routine to calculate the emittance if the distribution is Gaussian with cut tails at cutx
! just taking the standard deviation then gives a too small value
! calculate backwards from the standard deviation with a cut (see cut_Gaussian_tails_fit_emittance.nb for derivation)

  real :: x(nb),px(nb),emit,cutx,bx,sigObtained,sigma0,sigma1,func,funcP
  integer :: i,np

  ! write(6,*)
  ! write(6,*) 'emittance calc'
  ! write(6,*) '***************************'
  sigObtained=0.0
  do i=1,np
     sigObtained=sigObtained+x(i)**2+px(i)**2
  enddo
  sigObtained=sqrt(sigObtained/np/2.0) ! divide by 2 since we are including both x and px in the sum. this is now raw sigma^2
!  write(6,*) 'sigObtained=',sigObtained
!  write(6,*) 'obtained sigma from distribution = ',sigObtained

! now solve for original sigma without tails with Newton's method
! use obtained sigma as starting value 

  sigma0=sigObtained
  sigma1=-1.0
  do i=1,100 ! max 100 iterations
     func = cutx**2/(2-2*Exp(cutx**2/(2*sigma0**2))) + sigma0**2 - sigObtained**2
     funcP = 2*sigma0 - cutx**4 / ( 8*sigma0**3 *sinh(cutx**2/(4*sigma0**2)) )
     sigma1=sigma0-func/funcP
!     write(6,*) 'sigma0 = ',sigma0,'sigma1 = ',sigma1
     if (abs((sigma1-sigma0)/sigma0).le.0.005) exit
     sigma0=sigma1
  enddo
  if (i.eq.100) then
     write(6,*) 'emittance calculation not converging - STOP'
     stop
  endif
  emit=sigma1**2/bx
  ! write(6,*) 'emit = ',emit
  ! write(6,*) '**********************'
  ! write(6,*) 
  return
end subroutine calcEmitWCut

!*********************************************************************************
!*********************************************************************************

subroutine getemittance(x1,y1,px1,py1,np1,x2,y2,px2,py2,np2,ex1,ex2,ey1,ey2)
  include 'collider_time_evolution_common.f90'
  real :: x1(nb),px1(nb),y1(nb),py1(nb),x2(nb),px2(nb),y2(nb),py2(nb)
  real :: ex1,ey1,ex2,ey2
  integer :: np1,np2
  !     calculate emittance
  if (emitMethod.eq."stdev") then
     call emitStDev(x1,px1,np1,betax,ex1)
     call emitStDev(y1,py1,np1,betay,ey1)
     call emitStDev(x2,px2,np2,betax,ex2)
     call emitStDev(y2,py2,np2,betay,ey2)
  elseif (emitMethod.eq."stcal") then
     call calcEmitWCut(x1,px1,np1,xcut,betax,ex1)
     call calcEmitWCut(y1,py1,np1,ycut,betay,ey1)
     call calcEmitWCut(x2,px2,np2,xcut,betax,ex2)
     call calcEmitWCut(y2,py2,np2,ycut,betay,ey2)
  elseif (emitMethod.eq."exfit") then
     call fitemittance(x1,px1,np1,betax,cutoffAmpl,ex1)
     call fitemittance(y1,py1,np1,betay,cutoffAmpl,ey1)
     call fitemittance(x2,px2,np2,betax,cutoffAmpl,ex2)
     call fitemittance(y2,py2,np2,betay,cutoffAmpl,ey2)
  else
     write(6,*) 'unknown emittance routine. Stop.'
     stop
  endif
  ! write(6,*) 'emittance routine = ',emitMethod
  ! write(6,*) ex1
  return
end subroutine getemittance


!*********************************************************************************
!*********************************************************************************
     
      subroutine writecoord(np,x,px,y,py,t,pt,num,beamnr)
!     write 6D cordinates of all particles in distribution to file
!     for 50 000 particles, the file size is 4.3 MB

      include 'collider_time_evolution_common.f90'
      real :: x(nb),px(nb),t(nb),pt(nb),y(nb),py(nb)
      integer :: np,num,k,beamnr
      character(50) :: fname
      character(20) :: cnum
      character(5) :: beam
    
!     BEAM 1
!     construct filename
      write(cnum,'(i0)') num ! converting number to string
      write(beam,'(i0)') beamnr
      fname=trim('allcoord-B'//trim(beam)//'-'//trim(cnum)//'.out')
!     write(6,*) 'constructing filename: ',fname     
      
      open(55,file=fname)
      
      do k = 1,np
         write(55,*) y(k),py(k),x(k),px(k),t(k),pt(k)
      enddo

      close(55)
 
      return
      end


!*********************************************************************************
!*********************************************************************************

      subroutine zerokeep(ykeep,avgline,avglinex,avgliney,np)
      include 'collider_time_evolution_common.f90'
      integer :: k,np
      real :: avgline(nb),avglinex(nb),avgliney(nb),ykeep(7,nb)

      do k=1,np
         ykeep(1,k)=0
      enddo

      do k=1,nBins
         avgline(k)=0
         avglinex(k)=0
         avgliney(k)=0
      enddo
      return
      end

!*********************************************************************************
!*********************************************************************************

      subroutine rfupdate(np,x,px,t,pt,y,py,xpavg,fmix,nLostDebunch,iwrite)
      include 'collider_time_evolution_common.f90'
      real :: fmix,omega0,v00,volt,ptk,tk,p1,p2,ham
      real :: t(nb),pt(nb),xpavg(nb),thib2,x(nb),y(nb),px(nb),py(nb),pttot
      integer :: koo,np,k,iwrite,nLostDebunch


      omega0 = 2*pi/trev
      v00 = 938.e6*aatom/qatom
      koo = 0
      pttot=0.0
      thib2 = thib/2
      do k=1,np
         tk = t(k)
         ptk = pt(k)
! update the time

         tk = tk + fmix*tcoeff*ptk
         p1 = tk*nharm*omega0
         p2 = tk*nharm2*omega0
         volt = vrf*sin(p1)+vrf2*sin(p2)
         ptk = ptk + fmix*volt/v00
         pttot=pttot+ptk
! particles with abs(tk).gt.thib are lost. if thib==0, an infinite boundary is assumed and all particles are kept.
         if((abs(tk).le.thib2).or.(thib.eq.0.0)) then 
            koo = koo + 1
            t(koo)=tk
            pt(koo)=ptk
            x(koo)=x(k)
            y(koo)=y(k)
            px(koo)=px(k)
            py(koo)=py(k)
            ham=0.5*ptk*ptk*tcoeff+(cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)
! longitudinal hamiltonian
            xpavg(koo)=ham
         endif
      enddo
      nLostDebunch=nLostDebunch+np-koo
      np = koo
! put the normaization in np+1
      p1 = 0.5*thib*nharm*omega0
      p2 = 0.5*thib*nharm2*omega0
      xpavg(np+1)= (cos(p1)-1)*vrf/(nharm*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)
! write the average long. momentum as a check. radiation and RF should cancel, mean should be 0
!      if(iwrite.eq.1) write(16,*) 'mean pt=',pttot/real(np)
      return
      end

!*********************************************************************************
!*********************************************************************************

      subroutine raddamp(np,x,px,y,py,pt)
! does radiation damping and quantum excitation once per turn
      include 'collider_time_evolution_common.f90'
      integer :: np,k,np10,m
      real :: x(nb),px(nb),y(nb),py(nb),pt(nb),coeffdecaylong,coeffexcitelong,coeffgrow,coeffdecay,xrad,prad,ran3
      
! skip if transverse damping time is not positive
      if(tradperp.le.0)return

! damping time should scale down with the ratio of real turns/sim. turns
! exact damping coefficient is exp(-trev/tradlong). expand to first order

! get longitudinal radiation damping and excitation terms 
      coeffdecaylong = 1-trev/tradlong*timeRatio
! excitation uses a uniform deviate on [-1:1]
      coeffexcitelong = siglong*sqrt(3.)*sqrt(2*trev/tradlong*timeRatio)
! tradperp is the damping time for EMITTANCE, therefore need to multiply by 2
! assume same damping in horizontal and vertical plane (I4x,I4y<<I2)
      coeffdecay = 1 - trev/(2*tradperp)*timeRatio
! exact     coeffgrow= sigperp*sqrt(3.)*sqrt(1-coeffdecay**2)
! but trev << tradperp so 
      coeffgrow = sigperp*sqrt(3.)*sqrt(2*trev/(2*tradperp)*timeRatio)
!      write(6,*)coeffdecay,coeffgrow
 
      do k=1,np
! longitudinal
         pt(k) = pt(k)*coeffdecaylong +coeffexcitelong*(2*ran3(iseed)-1)
! transverse
         x(k) = coeffdecay*x(k) + (2*ran3(iseed)-1)*coeffgrow
         px(k) = coeffdecay*px(k)+(2*ran3(iseed)-1)*coeffgrow
         y(k) = coeffdecay*y(k) + (2*ran3(iseed)-1)*coeffgrow
         py(k) = coeffdecay*py(k)+(2*ran3(iseed)-1)*coeffgrow

      enddo
      return
      end


!*********************************************************************************
!*********************************************************************************


      subroutine sptransChrom(np,x,px,y,py,pt,fracturn)
! full turn (fracturn =1) update for x, horizontal
! fro chromx =chromy=0 should give same result as Mikes original routine
      include 'collider_time_evolution_common.f90'
      real :: x(nb),px(nb),y(nb),py(nb),pt(nb),fracturn
      real :: psix,psiy,a11,a12,tpdqmin,xk,pkx,xk1,pxk1,yk,pky,yk1,pyk1,coeffchromx,coeffchromy,psi1,ptk
      integer :: np,k,m

      psix = 2*pi*tunex*fracturn
      psiy = 2*pi*tuney*fracturn
      coeffchromx = fracturn*chromx*2*pi/(gamma0-1/gamma0)
      coeffchromy = fracturn*chromy*2*pi/(gamma0-1/gamma0)

!      tpdqmin = 2*pi*dqmin

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


!  *********** MSCHAUMANN ***************
! now have dqmin part - coupling between x and y
!         px(k)=pxk1+tpdqmin*yk1
!         py(k)=pyk1+tpdqmin*xk1

! seperate the coupling between x and y into an other subroutine sptransCoupling(), for Stochastic Cooling
         px(k) = pxk1
         py(k) = pyk1
         
      enddo
!      write(44,*)psi1,ptk*coeffchrom,ptk
      return
      end

!**********************************************************************************
!**********************************************************************************
!********** MSCHAUMANN *********

     subroutine sptransCoupling(np,x,px,y,py)
!Take out the dqmin part for the x-y skew quad-coupling from sptransChrom() 
!this is only needed once per trun, but the sptransChrom routine is called once after each Pickup and Kicker for a fraction of the turn if Stochastic Cooling is switched on
     include 'collider_time_evolution_common.f90'

     real :: x(nb),px(nb),y(nb),py(nb),tpdqmin
     integer :: k, np

     tpdqmin = 2*pi*dqmin
     do k =1, np
         px(k)=px(k)+tpdqmin*y(k)
         py(k)=py(k)+tpdqmin*x(k)
     enddo
     return
     end

!*********************************************************************************
!*********************************************************************************
      subroutine ibsPiwSmooth(pnumber,epsx,epsy,sigs,dponp,alfap0,alfax0,alfay0)

! uses Piwinski's formulae for ibs, pg 126 in handbook
! smooth lattice approximation

      include 'collider_time_evolution_common.f90' 

      real :: xdisp,sigh2inv,sigh,r0,atop,abot,rmsx,rmsy,&
           ceuler,fmohlp,fmohlx,fmohly,q,pnumber,epsx,epsy,sigs,dponp,&
           alfap0,alfax0,alfay0,ca,a,b,d,fmohl
      integer :: npp

! dispersion in smoooth approximation
      xdisp = circ/(2*pi*gammat**2)

      rmsx=sqrt(epsx*betax)
      rmsy=sqrt(epsy*betay)

      sigh2inv =  1/dponp**2  + (xdisp/rmsx)**2 
      sigh = 1/sqrt(sigh2inv)

! get capitol A
!      write(6,*)sigx,sigy,dponp,beta,sigh
! classical radius
      r0 = (qatom**2/aatom)*1.54e-18
      atop = r0*r0*clight*pnumber
      abot = 64*pi*pi*beta**3*gamma0**4*epsy*epsx*sigs*dponp
! ca is Piwinski's A
      ca = atop/abot
! mohl's a,b,and q
! a is horizontal, x
      a = sigh*betax/(gamma0*rmsx)
! b is vertical
      b = sigh*betay/(gamma0*rmsy)
! log is good enough
      if (rmsx.le.rmsy) then
         d=rmsx
      else
         d=rmsy
      endif
      q = sigh*beta*sqrt(2*d/r0)
! calculate fmohl(a,b,q) with 1000 points
      npp=1000
      ceuler = 0.577215
!      pi = 3.14159265
!      write(6,*)' a, b, q ', a,b,q
      fmohlp = fmohl(a,b,q,npp)
      fmohlx = fmohl(1/a,b/a,q/a,npp)
      fmohly = fmohl(1/b,a/b,q/b,npp)
      alfap0 =ca*fmohlp*(sigh/dponp)**2
      alfax0 =ca*(fmohlx+fmohlp*(xdisp*sigh/rmsx)**2)
      alfay0 =ca*fmohly
      return
      end
     

!*********************************************************************************
!*********************************************************************************
      subroutine ibsPiwLattice(pnumber,epsx,epsy,sigs,dponp,alfap0,alfax0,alfay0)

! uses Piwinski's formulae for ibs, pg 126 in handbook
! taking variations of optical functions in lattice into account

      include 'collider_time_evolution_common.f90' 

      real :: bx,by,xdisp,sigh2inv,sigh,r0,atop,abot,rmsx,rmsy,ceuler,fmohlp,fmohlx,fmohly,q,pnumber,epsx,epsy,sigs,dponp,&
           alfap0,alfax0,alfay0,ca,a,b,d,fmohl
      integer :: npp,i

!     classical radius
      r0 = (qatom**2/aatom)*1.54e-18

!     get capitol A
      atop = r0*r0*clight*pnumber
      abot = 64*pi*pi*beta**3*gamma0**4*epsy*epsx*sigs*dponp
! ca is Piwinski's A
      ca = atop/abot
      
      npp=1000
      ceuler = 0.577215

! zero out average quantities
      alfax0=0
      alfay0=0
      alfap0=0

      do i=1,nElem
         bx = betx(i)
         by = bety(i)
         xdisp = Dispx(i)

         rmsx=sqrt(epsx*bx)
         rmsy=sqrt(epsy*by)

         if (rmsx.le.rmsy) then
            d=rmsx
         else
            d=rmsy
         endif

         sigh2inv =  1/dponp**2  + (xdisp/rmsx)**2 
         sigh = 1/sqrt(sigh2inv)

! mohl's a,b,and q
! a is horizontal
         a = sigh*bx/(gamma0*rmsx)
! b is vertical
         b = sigh*by/(gamma0*rmsy)
         q = sigh*beta*sqrt(2*d/r0)
! calculate fmohl(a,b,q) with 1000 points

!      pi = 3.14159265
!      write(6,*)' a, b, q ', a,b,q
         fmohlp = fmohl(a,b,q,npp)
         fmohlx = fmohl(1/a,b/a,q/a,npp)
         fmohly = fmohl(1/b,a/b,q/b,npp)
         alfap0 =ca*fmohlp*(sigh/dponp)**2*leng(i)+alfap0
         alfax0 =ca*(fmohlx+fmohlp*(xdisp*sigh/rmsx)**2)*leng(i)+alfax0
         alfay0 =ca*fmohly*leng(i)+alfay0
      enddo
      alfap0=alfap0/circ
      alfax0=alfax0/circ
      alfay0=alfay0/circ

!$$$      write(6,*) 'sigH=',sigh,'ca=',ca,'fmohlp=',fmohlp
!$$$      write(6,*) 'alfx=',alfax0,'alfap0=',alfap0
!$$$      write(6,*) 'Tp=',1/alfap0/3600/2
!$$$      write(6,*) 'Tx=',1/alfax0/3600/2
!$$$      stop
      
      return
      end
     
!*********************************************************************************
!*********************************************************************************
      subroutine modPiwLattice(pnumber,epsx,epsy,sigs,dponp,alfap0,alfax0,alfay0)

! uses Piwinski's formulae for ibs, pg 126 in handbook
! taking variations of optical functions in lattice into account

      include 'collider_time_evolution_common.f90' 

      real :: bx,by,xdisp,sigh2inv,sigh,r0,atop,abot,rmsx,rmsy,ceuler,fmohlp,fmohlx,fmohly,q,pnumber,epsx,epsy,sigs,dponp,&
           alfap0,alfax0,alfay0,ca,a,b,d,fmohl,H,xdispP,bxP
      integer :: npp,i

!     classical radius
      r0 = (qatom**2/aatom)*1.54e-18

!     get capitol A
      atop = r0*r0*clight*pnumber
      abot = 64*pi*pi*beta**3*gamma0**4*epsy*epsx*sigs*dponp
! ca is Piwinski's A
      ca = atop/abot
      
      npp=1000
      ceuler = 0.577215

! zero out average quantities
      alfax0=0
      alfay0=0
      alfap0=0

      do i=1,nElem
         bx = betx(i)
         by = bety(i)
         xdisp = Dispx(i)
         xdispP = dispxP(i)
         bxP=-2*alfx(i)

         H=(xdisp**2+(bx*xdispP-0.5*xdisp*bxP)**2)/bx         

         rmsx=sqrt(epsx*bx)
         rmsy=sqrt(epsy*by)

         if (rmsx.le.rmsy) then
            d=rmsx
         else
            d=rmsy
         endif
         
         sigh2inv =  1/dponp**2 + H/epsx
         sigh = 1/sqrt(sigh2inv)

! mohl's a,b,and q
! a is horizontal
         a = sigh*bx/(gamma0*rmsx)
! b is vertical
         b = sigh*by/(gamma0*rmsy)
         q = sigh*beta*sqrt(2*d/r0)
! calculate fmohl(a,b,q) with 1000 points

         fmohlp = fmohl(a,b,q,npp)
         fmohlx = fmohl(1/a,b/a,q/a,npp)
         fmohly = fmohl(1/b,a/b,q/b,npp)
         alfap0 =ca*fmohlp*(sigh/dponp)**2*leng(i)+alfap0
         alfax0 =ca*(fmohlx+fmohlp*H*sigh**2/epsx)*leng(i)+alfax0
         alfay0 =ca*fmohly*leng(i)+alfay0

      enddo
      alfap0=alfap0/circ
      alfax0=alfax0/circ
      alfay0=alfay0/circ

!$$$      write(6,*) 'a=',a,'b=',b,'q=',q
!$$$      write(6,*) 'sigH=',sigh,'ca=',ca,'fmohlx=',fmohlx,'rmsx=',rmsx
!$$$      write(6,*) 'H=',H,'alfx=',alfax0,'alfap0=',alfap0
!$$$      write(6,*) 'Tp=',1/alfap0/3600/2
!$$$      write(6,*) 'Tx=',1/alfax0/3600/2
!$$$      stop
      
      return
      end

!*********************************************************************************
!*********************************************************************************
      subroutine ibsBane(nparti,ex,ey,sigs,sigp,alfap0,alfax0,alfay0)

!     high-energy IBS approximation. Reference: K.L.F Bane, A Simplified Model of Intrabeam Scattering, SLAC-PUB-9226 (2002)

      include 'collider_time_evolution_common.f90' 

      external gBane
      
      real :: ex,ey,r0,sigp,sigs,coulLog,nparti !input to be defined
      real :: bx,by,dxx,dyy,betxp,betyp ! twiss functions
      real :: Hx,Hy,sigH,a,b,sigx,sigy,d,Tp,HxAVG,HyAVG,brAVG,Tx,Ty,gBane,alfap0,alfax0,alfay0 ! parameters to be calculated

      integer :: i,j

!     classical radius
      r0 = (qatom**2/aatom)*1.54e-18

!     zero out average quantities
      brAVG=0  ! br is the quantity in brackets < > to be averaged in eq. 11
      HxAVG=0
      HyAVG=0

!     loop over all elements in lattice, calculate average quantities required by eqs 11,14 in slacpub9226 

      do i=1,nElem
         bx=betx(i)
         by=bety(i)
         dxx=Dispx(i)
         dyy=Dispy(i)
         betxp=-2*alfx(i)
         betyp=-2*alfy(i)
       
         sigx=sqrt(ex*bx+sigp**2*dxx)
         sigy=sqrt(ey*by+sigp**2*dyy)

         Hx=(Dxx**2+(bx*dispxP(i)-0.5*betxp*Dxx)**2)/bx
         Hy=(Dyy**2+(by*dispyP(i)-0.5*betyp*Dyy)**2)/by

         sigH=1/sqrt(1/sigp**2+Hx/ex+Hy/ey)
         a=sigH/gamma0*sqrt(bx/ex)
         b=sigH/gamma0*sqrt(by/ey)
      
         if (sigx<sigy) then
            d=sigx
         else
            d=sigy
         endif
         coulLog=log(d*sigH**2/(4*r0*a**2))

!     calcualte average over lattice of quantities required
         brAVG=brAVG+sigH*gBane(a/b)*(bx*by)**(-0.25)*coulLog*leng(i)  ! the factor to be averaged over in Eq.11, slacpub9226
         HxAVG=HxAVG+Hx*leng(i)
         HyAVG=HyAVG+Hy*leng(i)
      enddo

      brAVG=brAVG/circ
      HxAVG=HxAVG/circ
      HyAVG=HyAVG/circ

      Tp=1/(r0**2*clight*nparti/(16*gamma0**3*(ex*ey)**0.75*sigs*sigp**3)*brAVG)
      
      Tx=Tp*ex/(HxAvg*sigp**2)
      Ty=Tp*ey/(HyAvg*sigp**2)
      
      alfap0=1/Tp
      alfax0=1/Tx
      alfay0=1/Ty
      
      return
      end

!*********************************************************************************
      real function gBane(alpha)
      include 'collider_time_evolution_common.f90'
!     calculates value of Bane's g-function through linear interpolation of tabulated values in file
      integer :: nLow
      real :: alpha
      

      if((alpha<0.02).or.(alpha>=15.0)) then
         write(*,*) 'gBane: called with argument ',alpha,' outside range. allowed range is 0.02<arg<10'//CHAR(13)//CHAR(10)
         stop
      endif
!     find the closest lower tabulated value for alpha--nLow is the bin in the list where interpolation starts
      nLow=int(100*alpha+1)
      gBane=(gTab(nLow+1,2)-gTab(nLow,2))*(alpha-gTab(nLow,1))/(gTab(nLow+1,1)-gTab(nLow,1))+gTab(nLow,2)
      
      return
      end

!*********************************************************************************
!*********************************************************************************

      subroutine nagaitsevIBSlattice(pnumber,epsx,epsy,dponp,sigs,alfap0,alfax0,alfay0)
!     Authors: Roderik Bruce, Tom Mertens
!     calcualtes IBS growth rates using the Nagaitsev formulation of Bjorken-Mtingwa. Reference: PRSTAB 8, 064403 (2005)
      include 'collider_time_evolution_common.f90'
!      implicit none
   
      external rd_s

      integer :: i
     
      REAL(8) :: rd_s
      real :: pnumber,epsx,epsy,sigs,dponp,S,rIon,b1,alfap,alfax,alfay,r0 !
      real ::betaxP,phi,axx,ayy,sigmax,sigmay,as,a1,a2,lambda1,lambda2,lambda3 !
      real :: R1,R2,R3,sp,sx,sxp,alfap0,alfax0,alfay0,BMlog,alfapp,alfaxx,alfayy !

      rIon=(qatom**2/aatom)*1.54e-18
      alfap0=0.
      alfax0=0.
      alfay0=0.
 
!      write(*,*) nElem,rIon

      do i=1,nElem
         phi=dispxP(i)+(alfx(i)*(dispx(i)/betx(i)))
!         write(6,*) betx(i)
!         write(6,*) phi,nElem
         axx=betx(i)/epsx
         ayy=bety(i)/epsy
         sigmax=sqrt(dispx(i)**2*dponp**2+epsx*betx(i))
         sigmay=sqrt(epsy*bety(i))       
         as=axx*(dispx(i)**2/betx(i)**2+phi**2)+(1/(dponp**2))

         a1=0.5d0*(axx+gamma0**2*as)
         a2=0.5d0*(axx-gamma0**2*as)

         b1 = sqrt(a2**2 + gamma0**2*axx**2*phi**2)

         lambda1=ayy
         lambda2 = a1 + b1
         lambda3 = a1 - b1

         R1=(1/lambda1)*rd_s((1.d0/lambda2),(1.d0/lambda3),(1.d0/lambda1))
         R2=(1/lambda2)*rd_s((1.d0/lambda3),(1.d0/lambda1),(1.d0/lambda2))
         R3=3*sqrt((lambda1*lambda2)/lambda3)-(lambda1/lambda3)*R1-(lambda2/lambda3)*R2

         sp=(gamma0**2/2.d0) * ( 2.d0*R1 - R2*( 1.d0 - 3.d0*a2/b1 ) &
                                         - R3*( 1.d0 + 3.d0*a2/b1 ))

         sx=0.5d0 * (2.d0*R1 - R2*(1.d0 + 3.d0*a2/b1) &
                              -R3*(1.d0 - 3.d0*a2/b1))

         sxp=(3.d0*gamma0**2*phi**2*axx)/b1*(R3-R2)

         alfapp = sp/(sigmax*sigmay)

         alfaxx=(betx(i)/(sigmax*sigmay)) * (sx+sxp+sp*(dispx(i)**2/betx(i)**2+phi**2))

         alfayy=(bety(i)/(sigmax*sigmay)) * (-2.d0*R1+R2+R3)
         
         ! if(leng(i).gt.0.0) then
         !    write(44,*) 1/alphap,1/alphax,betax(i),betay(i),sigmax,sigmay,ax,ay,as,a1,a2,lambda1,lambda2,lambda3,R1,R2,R3,sp,sx,&
         !         sxp,alphap,alphax
         !    ! write(6,*) sp,sx,sxp,betax(i)
         !    ! stop
         ! endif

         alfap0 = alfapp*leng(i)/circ + alfap0
         alfax0 = alfaxx*leng(i)/circ + alfax0
         alfay0 = alfayy*leng(i)/circ + alfay0

         !                open(unit=80,file='carlsson.out',access='append')
         !                write(80,*) alphap0,alphax0,alphay0
         !                close(80)

      enddo
! extra division by 2 because later in program is multiplied with, comes from other models
! that simulate sigmas and not emittances
      alfap0 = alfap0 / dponp**2 * (pnumber*rIon**2*clight*coulombLog)/(12.d0*pi*beta**3*gamma0**5*sigs)/2
      alfax0 = alfax0 / epsx *     (pnumber*rIon**2*clight*coulombLog)/(12.d0*pi*beta**3*gamma0**5*sigs)/2
      alfay0 = alfay0 / epsy *     (pnumber*rIon**2*clight*coulombLog)/(12.d0*pi*beta**3*gamma0**5*sigs)/2

      ! write(6,*) 'pnumber=',pnumber,'epsx=',epsx,'epsy=',epsy,'dponp=',dponp,'sigs=',sigs
      ! write(6,*) 'alfax0=',alfax0,'alfay0=',alfay0,'alfap0=',alfap0
      ! write(6,*) 'Tp=',1/alfap0/3600/2.
      ! write(6,*) 'Tx=',1/alfax0/3600/2.
      ! write(6,*) 'Ty=',1/alfay0/3600/2.
      
      return
    end subroutine nagaitsevIBSlattice

!**********************************************************************************
!**********************************************************************************

FUNCTION rd_s(x,y,z)
!USE nrtype; USE nrutil, ONLY : assert
!IMPLICIT NONE
REAL(8), INTENT(IN) :: x,y,z
REAL(8) :: rd_s
REAL(8), PARAMETER :: ERRTOL=0.05,TINY=1.0e-25,BIG=4.5e21,&
C1=3.0/14.0,C2=1.0/6.0,C3=9.0/22.0,&
C4=3.0/26.0,C5=0.25*C3,C6=1.5*C4
!Computes Carlson elliptic integral of the second kind, RD(x, y, z). x and y must be
!nonnegative, and at most one can be zero. z must be positive. TINY must be at least twice
!the negative 2/3 power of the machine overflow limit. BIG must be at most 0.1×ERRTOL
!times the negative 2/3 power of the machine underflow limit.
REAL(4) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,&
ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
integer :: iter
! call assert(min(x,y) >= 0.0, min(x+y,z) >= TINY, max(x,y,z) <= BIG,'rd_s args')
xt=x
yt=y
zt=z
sum=0.0
fac=1.0
iter=0
do
   iter=iter+1
   sqrtx=sqrt(xt)
   sqrty=sqrt(yt)
   sqrtz=sqrt(zt)
   alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
   sum=sum+fac/(sqrtz*(zt+alamb))
   fac=0.25*fac
   xt=0.25*(xt+alamb)
   yt=0.25*(yt+alamb)
   zt=0.25*(zt+alamb)
   ave=0.2*(xt+yt+3.0*zt)
   delx=(ave-xt)/ave
   dely=(ave-yt)/ave
   delz=(ave-zt)/ave
   if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
end do
ea=delx*dely
eb=delz*delz
ec=ea-eb
ed=ea-6.0*eb
ee=ed+ec+ec
rd_s=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
END FUNCTION rd_s

     
!*********************************************************************************
!*********************************************************************************

      subroutine ibsInterpolat(pnumber,epsx,epsy,sigs,dponp,alfap0,alfax0,alfay0)
      include 'collider_time_evolution_common.f90'
!     trilinear interpolation of ibs growth rates Ai(x,y,z) , for i=x,y,z, from tabulated values in external file 
!     reading A as a 3D array A(i,j,k)
!     reference: see for example http://en.wikipedia.org/wiki/Trilinear_interpolation

      integer :: i,j,k
      real :: a000,a001,a010,a100,a110,a101,a011,a111,a00,a01,a10,a11,a0,a1,alfap0,alfax0,alfay0,invDeltaX,invDeltaZ,invDeltaY,&
           DeltaX,DeltaY,DeltaZ,x0,y0,z0,pnumber,epsx,epsy,sigs,dponp

      real :: x,y,z ! point where function value should be interpolated

!     rewrite input arguments as x,y,z for use in trilinear program
      x=(epsx+epsy)/2 ! average transverse emittance
      y=sigs
      z=dponp

!     if arguments outside tabulated range, exit
      if ((x>xmax).or.(x<xmin).or.(y>ymax).or.(y<ymin).or.(z>zmax).or.(z<zmin)) then
         write(*,*) 'input argument out of bounds'
         write(*,*) 'x=',x,' xmin= ',xmin,' xmax=',xmax
         write(*,*) 'y=',y,' ymin= ',ymin,' ymax=',ymax
         write(*,*) 'z=',z,' zmin= ',zmin,' zmax=',zmax
         stop
      endif

!     find correct bins in 000 corner
      i=int(1+(x-xMin)/(xMax-xMin)*(xBins-1))
      j=int(1+(y-ymin)/(ymax-ymin)*(yBins-1))
      k=int(1+(z-zmin)/(zmax-zmin)*(zBins-1))

!     calculate bin width in all 3 directions. 
      DeltaX=(xMax-xMin)/(xBins-1)
      DeltaY=(yMax-yMin)/(yBins-1)
      DeltaZ=(zMax-zMin)/(zBins-1)

      invDeltaX=1/DeltaX
      invDeltaY=1/DeltaY
      invDeltaZ=1/DeltaZ

!     calculate the lower values of (x,y,z) in the corner of the cube
      x0=xmin+DeltaX*(i-1)
      y0=ymin+DeltaY*(j-1)
      z0=zmin+DeltaZ*(k-1)

!     interpolate alfap
!------------------------------------------------------
!     get tabulated function values in corners of cube
      a000=Ap(i,j,k)
      a100=Ap(i+1,j,k)
      a010=Ap(i,j+1,k)
      a001=Ap(i,j,k+1)
      a110=Ap(i+1,j+1,k)
      a101=Ap(i+1,j,k+1)
      a011=Ap(i,j+1,k+1)
      a111=Ap(i+1,j+1,k+1)

!     find 4 new T-values by interpolation along x (middle of each side along x in cube)
      a00=a000+(a100-a000)*invDeltaX*(x-x0)
      a10=a010+(a110-a010)*invDeltaX*(x-x0)
      a01=a001+(a101-a001)*invDeltaX*(x-x0)
      a11=a011+(a111-a011)*invDeltaX*(x-x0)

!     find 2 new T-values by interpolation along y (middle points of the two planes with different y-values)
      a0=a00+(a10-a00)*invDeltaY*(y-y0)
      a1=a01+(a11-a01)*invDeltaY*(y-y0)

!     find final T-value by interpolation along z to the point inside cube
      alfap0=a0+(a1-a0)*invDeltaZ*(z-z0)

!     interpolate alfax
!------------------------------------------------------
!     get tabulated function values in corners of cube
      a000=Ax(i,j,k)
      a100=Ax(i+1,j,k)
      a010=Ax(i,j+1,k)
      a001=Ax(i,j,k+1)
      a110=Ax(i+1,j+1,k)
      a101=Ax(i+1,j,k+1)
      a011=Ax(i,j+1,k+1)
      a111=Ax(i+1,j+1,k+1)

!     find 4 new T-values by interpolation along x (middle of each side along x in cube)
      a00=a000+(a100-a000)*invDeltaX*(x-x0)
      a10=a010+(a110-a010)*invDeltaX*(x-x0)
      a01=a001+(a101-a001)*invDeltaX*(x-x0)
      a11=a011+(a111-a011)*invDeltaX*(x-x0)

!     find 2 new T-values by interpolation along y (middle points of the two planes with different y-values)
      a0=a00+(a10-a00)*invDeltaY*(y-y0)
      a1=a01+(a11-a01)*invDeltaY*(y-y0)

!     find final T-value by interpolation along z to the point inside cube
      alfax0=a0+(a1-a0)*invDeltaZ*(z-z0)

!     interpolate alfay
!------------------------------------------------------
!     get tabulated function values in corners of cube
      a000=Ay(i,j,k)
      a100=Ay(i+1,j,k)
      a010=Ay(i,j+1,k)
      a001=Ay(i,j,k+1)
      a110=Ay(i+1,j+1,k)
      a101=Ay(i+1,j,k+1)
      a011=Ay(i,j+1,k+1)
      a111=Ay(i+1,j+1,k+1)

!     find 4 new T-values by interpolation along x (middle of each side along x in cube)
      a00=a000+(a100-a000)*invDeltaX*(x-x0)
      a10=a010+(a110-a010)*invDeltaX*(x-x0)
      a01=a001+(a101-a001)*invDeltaX*(x-x0)
      a11=a011+(a111-a011)*invDeltaX*(x-x0)

!     find 2 new T-values by interpolation along y (middle points of the two planes with different y-values)
      a0=a00+(a10-a00)*invDeltaY*(y-y0)
      a1=a01+(a11-a01)*invDeltaY*(y-y0)

!     find final T-value by interpolation along z to the point inside cube
      alfay0=a0+(a1-a0)*invDeltaZ*(z-z0)

!     rescale the bunch population from the fixed value used in tabulation
      alfap0=alfap0*pnumber/pnumInterp
      alfax0=alfax0*pnumber/pnumInterp
      alfay0=alfay0*pnumber/pnumInterp

      return
      end



!*********************************************************************************
!*********************************************************************************

      subroutine ibslong(np,y,py,t,pt,x,px,epsx,epsy,np0,pnumber,pnnorm,avgline,iwrite)
!     pnumber is the total number of particles in the bunch
!     pnnorm is the bunch population that IBS should be normalized to.
!     In the original code, the IBS strength is a factor pnumber/np0 higher, 
!     but we need the same normalization factor for BOTH bunches, otherwise 
!     their equivalent time in the machine will not match!

      include 'collider_time_evolution_common.f90'
      external ran3
      
      real :: rmsy,rmsx,rmsg,rmst,denlon2(nb),thib2,dtsamp2,tk,gk,yk,&
           pk,dponp,xdisp,sigh2inv,sigh,ro,atop,&
           epsy,epsx,abot,sigs,ca,a,b,q,ceuler,fmohl,fmohlp,fmohlx,&
           fmohly,alfap0,alfax0,alfay0,alfap,alfay,alfax,coeffs,&
           coeffy,coeffx,coeffmult,denlonn,denlon2k,dy,dg,dx,clum,&
           grv1,grv2,ci,pnumber,pnnorm,r0,alfaAvg,d,a1,a2

!     this routine copied as is from Mike, except small changes: 
!     --using timeRatio parameter to scale the strength instead of pnumber1/np10
!     --changed notation x2->x and x->y for consistency with the rest of the code
!     --implemented option of coupling that averages the growth rates over x and y

      real :: y(nb),py(nb),t(nb),pt(nb),x(nb),px(nb),avgline(nb) ! local beam variables
      integer :: np,np0,k,nk,npp,kkk,iwrite

      rmsg=0
      rmst=0
! zero out the array for line density
      do k=1,nBins
         denlon2(k)=0
      enddo
      dtsamp2 = thib/nBins
      thib2 = thib/2
      do k=1,np
!     t=0 is stable fixed point
         tk = t(k)
         gk = pt(k)
         rmst = rmst + tk*tk
         rmsg = rmsg + gk*gk
         tk = (tk+thib2)/dtsamp2
! just use nearest integer
         nk = nint(tk)
         if(nk.le.0)nk=1
         if(nk.gt.nBins)nk=nBins
         denlon2(nk)=denlon2(nk)+1
      enddo
! add into avgline(k)
      do k=1,nBins
         avgline(k) = avgline(k)+denlon2(k)
      enddo
      
! nothing for no ibs
      if(fracibstot.le.0)return

      rmst = sqrt(rmst/np)
      rmsg = sqrt(rmsg/np)
      dponp = rmsg/(beta*beta*gamma0)
      sigs = rmst*clight*beta
     
      rmsx = sqrt(epsx*betax)
      rmsy = sqrt(epsy*betay)

      if (ibsMethod.eq.'piwiSmooth') then
         call ibsPiwSmooth(pnumber,epsx,epsy,sigs,dponp,alfap0,alfax0,alfay0)
      elseif (ibsMethod.eq.'piwLattice') then
         call ibsPiwLattice(pnumber,epsx,epsy,sigs,dponp,alfap0,alfax0,alfay0)
      elseif (ibsMethod.eq.'modPiwLatt') then
         call modPiwLattice(pnumber,epsx,epsy,sigs,dponp,alfap0,alfax0,alfay0)
      elseif (ibsMethod.eq.'baneApprox') then
         call ibsBane(pnumber,epsx,epsy,sigs,dponp,alfap0,alfax0,alfay0)         
      elseif (ibsMethod.eq.'interpolat') then
         call ibsInterpolat(pnumber,epsx,epsy,sigs,dponp,alfap0,alfax0,alfay0)
      elseif (ibsMethod.eq.'nagaitsev') then
         call nagaitsevIBSlattice(pnumber,epsx,epsy,dponp,sigs,alfap0,alfax0,alfay0)
      else
         write(*,*) 'Stop - unknown ibs method: ',ibsMethod
         stop
      endif

!      write(*,*) alfap0,alfax0,alfay0
      alfap0=2*alfap0*float(np)/float(np0)
      alfax0=2*alfax0*float(np)/float(np0)
      alfay0=2*alfay0*float(np)/float(np0)
      

      if (iwrite.eq.1) then
         open(unit=80,file='ibs.out',access='append')
         write(80,'(I10,4E17.7)') kturns,real(kturns)/real(nturns)*eqTime,1/alfap0/3600,1/alfax0/3600,1/alfay0/3600
         close(80)
      endif

!      write(69,*) np,rmsy,rmsx,rmst,rmsg,alfap0,alfax0,alfay0
! alphas are amplitude growth rates for pnumber particles
! and no losses
! correct for small number of macro particles, normalize to pnnorm
! mult by 2 so it a kick for emittance 
! reduce by input fraction
!      write(6,*)ca,sigh,fmohl,dponp
      alfap = fracibstot*timeRatio*alfap0
!      write(6,*)alfap
      alfay = fracibstot*timeRatio*alfay0
      alfax = fracibstot*timeRatio*alfax0
! with rms growth rate for 1 dim of shm need to kick
! sqrt(2.) times harder
      if(alfap.gt.0)then
          coeffs = sqrt(6*alfap*trev)*rmsg
      else
          coeffs=0
      endif
      if(coupleIBS.eq.0) then   ! use separate growth rates in x and y. Original version from Mike
         if(alfay.gt.0)then
            coeffy = sqrt(6*alfay*trev)*rmsy
         else
            coeffy = 0
         endif
         if(alfax.gt.0)then
            coeffx = sqrt(6*alfax*trev)*rmsx
         else
            coeffx = 0
         endif
      else                      ! use average alfa in x and y. corresponds to COOL07 paper and early LHC simulations
         alfaAvg=(alfay+alfax)/2
         if (alfaAvg.gt.0) then
            coeffy=sqrt(6*alfaAvg*trev)*rmsy
            coeffx = sqrt(6*alfaAvg*trev)*rmsx
         else
            coeffy=0
            coeffx=0
         endif
      endif

      coeffmult = sigs*2*sqrt(pi)/(np*dtsamp2*clight)
      denlonn=0
      do k=1,np
          tk = (t(k)+thib2)/dtsamp2
! just use nearest integer
          nk = nint(tk)
          if(nk.le.0)nk=1
          if(nk.gt.nBins)nk=nBins
          
          denlon2k = denlon2(nk)*coeffmult
          denlonn = denlonn+denlon2k
          denlon2k = sqrt(denlon2k)
! that easy??
!         dy = denlon2k*coeffy*(2*ran3(iseed)-1)
!         dg = denlon2k*coeffs*(2*ran3(iseed)-1)
         call get2gaussrv(iseed,grv1,grv2)
         dy = denlon2k*coeffy*grv1
         dg = denlon2k*coeffs*grv2
         call get2gaussrv(iseed,grv1,grv2)
         dx = denlon2k*coeffx*grv1
  
         py(k)=py(k)+dy
         px(k)=px(k)+dx
         pt(k)=pt(k)+dg
      enddo
! get the central luminosity density
      clum=0
      do kkk=1,nBins
         clum = clum + denlon2(kkk)**2
      enddo
      ci = pnumber/np
      clum = clum*ci*ci/(clight*dtsamp2*rmsy**2)

      return
 100  format(i8,1p12e14.4)
      end


!*********************************************************************************
!*********************************************************************************

      subroutine get2gaussrv(iseed,grv1,grv2)
! output is two zero mean gaussian rvs with rms = 1/sqrt(3.)
! get two uniform deviates in the unit circle
      real :: grv1,grv2,r1,r2,facc,amp
      integer :: iseed

 44     continue
! r1 and r2 are uniform on (-1,1)
           r1 = 2*ran3(iseed)-1
           r2 = 2*ran3(iseed)-1
           amp = r1*r1+r2*r2
! just keep the circle and cut at 6 sigma
         if((amp.ge.1).or.(amp.lt. 1.e-8))go to 44
         facc = sqrt(-2.*log(amp)/amp)
         grv1 = r1*facc/1.73205         
         grv2 = r2*facc/1.73205
      return
      end


!*********************************************************************************
!*********************************************************************************

      function fmohl(a,b,q,np)
! direct crib from page 126 in handbook
      ceuler = 0.577215
      pi = 3.14159265
      sum = 0
      du = 1/float(np)
      do k=0,np
         u = k*du
         cp = sqrt(a*a+(1-a*a)*u*u)
         cq = sqrt(b*b+(1-b*b)*u*u)
         dsum = 2*log(q*(1/cp+1/cq)/2) - ceuler
         dsum = dsum*(1-3*u*u)/(cp*cq)
         if(k.eq.0)dsum=dsum/2
         if(k.eq.np)dsum=dsum/2
         sum = sum + dsum
      enddo
      fmohl = 8*pi*du*sum
!      write(6,*)'a,b,q,fmohl',a,b,q,fmohl
      return
      end

!*********************************************************************************
!*********************************************************************************

      subroutine keepturn(np,y,py,ykeep)
      include 'collider_time_evolution_common.f90'
      integer :: k,np
      real :: y(nb),py(nb),ykeep(7,nb),yk,pk

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


      subroutine writemomentsshort(np,x,px,t,pt,x2,px2,beamnr,coherecs,coherecss,tlob,&
                    pnumber,np0,dpalmer,dpxtot,dpytot,dpwakex,dpwakey,dpwakes,xkeep,xpavg)
      include 'collider_time_evolution_common.f90'
      integer :: kdex(0:100),iperm(nb)
      real :: xavg(100),pavg(100),csavg(100),zavg(100),torder(nb),sigx(100),sigpx(100),denlonavg(100),denlonsig(100),coherecss
          
      real :: x(nb),px(nb),t(nb),pt(nb),x2(nb),px2(nb)
      
      real :: csfull,csfull2,siggama,siggama2,sigtime,sigtime2,xk,pxk,pxk2,tk,dtsamp,cperbin,denk,xdenk,pxdenk,coherecsden,pxoffk,&
           amp_res,xoffk,xk2,coherecs

! xkeep and xtune only for output     
      real :: pnumber,dpalmer(nb),dpxtot(nb),dpytot(nb),dpwakex(nb),dpwakey(nb),dpwakes(nb),xtune(nb),xkeep(7,nb),tlob, xpavg(nb)

!*MS* variables to track kicker voltage evolution
      real:: coeffvoltage, peakdpwakex, peakdpwakey, peakdpwakes, rmsdpwakex, rmsdpwakey, rmsdpwakes, &
             meandpwakex, meandpwakey, meandpwakes

      integer :: k,k0,np,np0,beamnr,npCen,npC2,n
      character(5) :: beam
      character(50) :: fname

!*MS*
      peakdpwakex = 0
      peakdpwakey = 0
      peakdpwakes = 0
      rmsdpwakex = 0
      rmsdpwakey = 0
      rmsdpwakes = 0
      meandpwakex = 0
      meandpwakey = 0
      meandpwakes = 0
      coeffvoltage = aatom*938.e6/qatom *sqrt(np0/pnumber)

      csfull=0
      csfull2=0
      siggama=0
      siggama2=0
      sigtime=0
      sigtime2=0
      npCen=0
      npC2=0
      do k0 = 1,np
         xk = x(k0)
         pxk = px(k0)
         xk2 = x2(k0)
         pxk2 = px2(k0)
         csfull = csfull + xk*xk + pxk*pxk
         csfull2 = csfull2 + xk2*xk2 + pxk2*pxk2
! longitudinal
         siggama = siggama + pt(k0)**2
         tk = t(k0)
         sigtime = sigtime +tk*tk
         sigtime2 = sigtime2 + tk
         if (abs(tk)<2.5e-9) npCen=npCen+1   ! hard coding number manually for check
         if (abs(tk)<thib/4) npC2=npC2+1 ! what is really inside the h=360 bucket

      enddo
      sigtime = sqrt( sigtime/np - (sigtime2/np)**2)
      siggama = sqrt(siggama/np)

!     get the smoothed distributions
!     call getsmoothed(nlag)

      csfull = csfull/np
      csfull2 = csfull2/np

! stop the code if the beam is large enough
      if(csfull.gt.20)then
         write(6,*)'sumcsavg too big'
         stop
      endif

!     count the number of particles in the central bucket

 103  format(i10,1p10e13.5)
 100  format(1p12e15.5)

 104  format('kturns,np,log10(csfull),log10(csfull2),coherecs',2i10,1p10e11.3)
 105  format('kturns,np,csfull,csfull2,coherecs',2i10,1p10e11.3)

!     not calculating and writing the cooling and impedance stuff - 0 instead
!     construct filename
      write(beam,'(i0)') beamnr
      fname=trim('csmon-B'//trim(beam)//'.out')
      open(unit=33,file=fname,access='append')
!      write(33,103)kturns,csfull,csfull2,0.,sigtime,siggama,0.,float(npC2),float(np),float(npCen)
      write(33,103)kturns,csfull,csfull2,coherecs,sigtime,siggama,coherecss,float(npC2),float(np),float(npCen)
      close(33)

!**********   MSCHAUMANN ****** Put back impeadance stuff

! write the densities used for the impedance work
       open(unit=20,file='res_den.out',status='unknown')
!       coherecs=0
       coherecsden=0
! need to smooth the arrays
       dtsamp = tlob/nresamp
! cperbin = const. = intitial charge / nmacro/fft_bin (in units of Ampere]
       cperbin = pnumber*1.602e-19*qatom/(np0*dtsamp)
      !! write(6,*)'nresamp = ',nresamp
       do k0=1,nresamp
          denk = dpalmer(k0)
          xdenk = dpxtot(k0)
          pxdenk = dpytot(k0)
          amp_res = cperbin*denk
! **MS** long. position in bunch [s], position at long pickup, position at hor PU, position at ver PU, initial Charge/(nmacro) per bin * long pos. at PU, ... , ..., ...
! ??? cperbin should change with changing particle density per slice, since amp_res is porportional to the strength of the longitudinal kick (?)  
          write(20,100) (k0-1)*dtsamp,denk,xdenk,pxdenk,amp_res,dpwakex(k0),dpwakey(k0),dpwakes(k0)

          coherecsden = coherecsden + denk

!*MS* get peak kicker voltage for the current turn
         if(peakdpwakex<dpwakex(k0))  peakdpwakex=dpwakex(k0)
         if(peakdpwakey<dpwakey(k0))  peakdpwakey=dpwakey(k0)
         if(peakdpwakes<dpwakes(k0))  peakdpwakes=dpwakes(k0)

!*MS* get sums for mean and standard deviation
         rmsdpwakex = rmsdpwakex + dpwakex(k0)**2
         rmsdpwakey = rmsdpwakey + dpwakey(k0)**2
         rmsdpwakes = rmsdpwakes + dpwakes(k0)**2

         meandpwakex = meandpwakex + dpwakex(k0)
         meandpwakey = meandpwakey + dpwakey(k0)
         meandpwakes = meandpwakes + dpwakes(k0)


       enddo
       close(20)
!       write(6,104)kturns,np,log10(csfull),log10(csfull2),coherecs,coherecss
      !! write(6,105)kturns,np,csfull,csfull2,coherecs,coherecss

!*MS* calculate standard deviation and mean of kicker voltage   
        rmsdpwakex = coeffvoltage*sqrt(rmsdpwakex/nresamp - (meandpwakex/nresamp)**2)
        rmsdpwakey = coeffvoltage*sqrt(rmsdpwakey/nresamp - (meandpwakey/nresamp)**2)
        rmsdpwakes = coeffvoltage*sqrt(rmsdpwakes/nresamp - (meandpwakes/nresamp)**2)
 
        meandpwakex = coeffvoltage*meandpwakex/nresamp
        meandpwakey = coeffvoltage*meandpwakey/nresamp
        meandpwakes = coeffvoltage*meandpwakes/nresamp

!*MS* open file and write kicker voltages
       fname=trim('kickerVoltage-B'//trim(beam)//'.out')
       open(unit=36,file=fname,access='append')
       write(36,103) kturns, coeffvoltage*peakdpwakex, meandpwakex, rmsdpwakex, coeffvoltage*peakdpwakey, &
                     meandpwakey, rmsdpwakey, coeffvoltage*peakdpwakes, meandpwakes, rmsdpwakes 
       close(36)

       coherecs=0
       coherecss=0

! write everything to tran.full if nwrite.lt.0
       if(nwrite.gt.0)return
! order the output
       do k=1,np
          iperm(k)=k
! order with respect to momentum
          torder(k)=pt(k)
       enddo
       call sort2(np,torder,iperm)
       open(unit=20,file='tran_full.out',status='unknown')

! include subroutine tunestuff here
! xtune and ykeep only for output
      do n=1,np
        xtune(n) = xkeep(1,n)/nwrite
      enddo

       do k0=1,np
          k = iperm(k0)
! xtune is actually the action
          write(20,100)x(k),px(k),t(k),pt(k),xtune(k),xpavg(k)
       enddo
       close(20)

!**********   End cooling stuff   **********************


      return
      end

!*********************************************************************************
!*********************************************************************************


      subroutine mountainr(np,x,px,t,pt,x2,px2,beamnr,np0,pnumber,avgline)
! makes a mountain range plot of emittance versus 
! dgamma
      include 'collider_time_evolution_common.f90'

      real :: dt2,currcoeff,peaki,currk,pk,fwhmlumi,peakilumi,pnumber
      real :: x(nb),px(nb),t(nb),pt(nb),x2(nb),px2(nb),avgline(nb)

      integer :: np,np0,beamnr,jdex,jhi,jlo,k
      character(5) :: beam
      character(50) :: fname

      dt2 = thib/nBins
! put into units of current
! old      currcoeff = 1.602e-19*pnumber*qatom/(np0*nwrite0*dt2)
      currcoeff = 1.602e-19*pnumber*qatom/(nwrite0*dt2)
      currcoeff = currcoeff/np0
!      write(6,*)currcoeff,nwrite0,np0
      peaki=0
      jdex=0
      do k=1,nBins
         currk = avgline(k)*currcoeff
         if(currk.gt.peaki)then
            peaki=currk
            jdex=k
         endif
         avgline(k)=currk
      enddo
! get the full width half max of the peak
      do k=jdex,nBins
           pk = avgline(k)
           if(pk.le.peaki/2)go to 33
      enddo
 33   continue
      jhi=k
      do k=jdex,1,-1
           pk = avgline(k)
           if(pk.le.peaki/2)go to 34
      enddo
 34   continue
      jlo=k
      fwhmlumi = dt2*(jhi-jlo+1)
      peakilumi=peaki

!     construct filename
      write(beam,'(i0)') beamnr
      fname=trim('mountain-B'//trim(beam)//'.out')
      open(unit=33,file=fname,access='append')

      do k=1,nBins
         write(33,100)k*dt2,avgline(k),float(kturns)
! zero it out
         avgline(k)=0
      enddo
      write(33,*)
 100  format(1p10e14.4)
      close(33)
      return
      end


!*********************************************************************************
!*********************************************************************************
      subroutine keepTransvProf(avglinex,avgliney,x,y,np,emix,emiy)
! bin the transverse coord., store the  profile in arrays
      include 'collider_time_evolution_common.f90'
      real :: x(nb),y(nb),avglinex(nb),avgliney(nb),xmn,dx,emix,ymn,dy,emiy
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
      subroutine writeNb(np1,np2,np10,np20,pnumber1,pnumber2,nLostLum1,nLostLum2,nLostLumSum1,nLostLumSum2,&
           nLostDebunch1,nLostDebunch2,nLostDebunchSum1,nLostDebunchSum2,&
           nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1,nLostBeta2,nLostBetaSum2,nLostMom2,nLostMomSum2)

      include 'collider_time_evolution_common.f90'

      integer :: np1,np2,np10,np20,nLostLum1,nLostLum2,nLostLumSum1,nLostLumSum2
      integer :: nLostDebunch1,nLostDebunch2,nLostDebunchSum1,nLostDebunchSum2
      integer :: nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1,nLostBeta2,nLostBetaSum2,nLostMom2,nLostMomSum2
      real :: pnumber1,pnumber2
      
      nLostLumSum1=nLostLumSum1+nLostLum1 ! total n.o. lost particles over all turns from lumi
      nLostLumSum2=nLostLumSum2+nLostLum2
      nLostDebunchSum1=nLostDebunchSum1+nLostDebunch1 ! total n.o. lost particles over all turns from debunching
      nLostDebunchSum2=nLostDebunchSum2+nLostDebunch2

      nLostBetaSum1=nLostBetaSum1+nLostBeta1
      nLostMomSum1=nLostMomSum1+nLostMom1
      nLostBetaSum2=nLostBetaSum2+nLostBeta2
      nLostMomSum2=nLostMomSum2+nLostMom2
      
      open(unit=80,file='intensity.out',access='append')
      write(80,'(I10,E17.7,I10,E17.7,9I10,E17.7,8I10)') &
           kturns,real(kturns)/real(nturns)*eqTime, &
           np1,real(np1)/real(np10)*pnumber1, nLostLum1,nLostLumSum1,nLostDebunch1,nLostDebunchSum1, &
           nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1, &
           np2,real(np2)/real(np20)*pnumber2, nLostLum2,nLostLumSum2,nLostDebunch2,nLostDebunchSum2, &
           nLostBeta2,nLostBetaSum2,nLostMom2,nLostMomSum2
      close(80)      
      nLostLum1=0 ! reset intermediate sums
      nLostLum2=0
      nLostDebunch1=0
      nLostDebunch2=0
      nLostBeta1=0
      nLostBeta2=0
      nLostMom1=0
      nLostMom2=0
      return
    end subroutine writeNb

!*********************************************************************************
!*********************************************************************************
subroutine writeLumi(lumi,redfac)
  include 'collider_time_evolution_common.f90'
  
  real :: lumi,redfac
  
  write(6,*) 'luminosity=',lumi,'cm^-2 s^-1, ','geom. reduction factor=',redfac
  open(unit=80,file='luminosity.out',access='append')
  write(80,'(I10,4E17.7)') kturns,real(kturns)/real(nturns)*eqTime,lumi,redfac,betaS(1)
  close(80)
!*********** MSCHAUMA ************ include betaS in output for IP1

  return
end subroutine writeLumi

!*********************************************************************************
!*********************************************************************************

!---------------------------------------------------------------------------------
! by MICHAELA SCHAUMANN:
!
! Calculate the FWHM out of the coordinates to be able to compare correctely with 
! the LHC BQM measurement which gives a logged value of 
! tau = (2/ln(2))^(1/2) FWHM
! which corresponds to 4 std deviations in case of a gaussian beam profile
!
	subroutine calculateFWHM(t,np,fwhm,tauhat)
 	include 'collider_time_evolution_common.f90'
	
	real :: t(nb),fwhm,tauhat,binWidth,histo(50),&
		mLeft,mRight, bLeft,bRight,xLeft,xRight
	integer :: np,k,numberBins,maximum,xHalfMin,xHalfMax,&
		   binIndex(nb), maxIndex,&
		   xHalfMin1,xHalfMax1,xHalfMin2,xHalfMax2

	fwhm = 0.

	!devide phase space in numberBins bins with the width binWidth to create a histogram
	numberBins = 50
	binWidth = 2*tauhat/numberBins

	!set binIndex to zero everywhere to use it as a counting array
	do k=1,np
	binIndex(k) = 0
	enddo
	do k=1,numberBins
	histo(k) = 0
	enddo

	!count how many paritcles are in which bin of the phase space
	!and fill histogram
	!calculate fwhm

	do k=1,np
	   binIndex(k) = int(((t(k)+tauhat)/binWidth )+1)
	   histo(binIndex(k)) = histo(binIndex(k))+1
	enddo

!	!11/01/2013 normalise the histogram to the number of nMacro
!	do k=1, numberBins
!	   histo(k) = histo(k)/np
!	enddo

	!find the maximum value of the histogram
	maximum = -1 !set to a value which can not appear
	maxIndex =-1
	do k=1, numberBins
	  if (maximum.lt.histo(k)) then
	    maximum = histo(k)
	    maxIndex = k
	  endif
	enddo
	maximum = (histo(maxIndex-1) + histo(maxIndex) + histo(maxIndex+1))/3

! ***	15/01/2013
!	Improve the calculation of the FWHM by doing an interpolation between the bins just smaller and larger as max/2
!	Search for the position of maximum/2 - save the bin which is just below and above that value to both sides
	xHalfMin1 = 1
	xHalfMin2 = maxIndex
	xHalfMax2 = maxIndex
	xHalfMax1 = numberBins
	do while (histo(xHalfMin1).lt.maximum/2)
	  xHalfMin1 = xHalfMin1 +1
	enddo
	do while (histo(xHalfMin2).gt.maximum/2)
	  xHalfMin2 = xHalfMin2 -1
	enddo
	do while (histo(xHalfMax2).gt.maximum/2)
	  xHalfMax2 = xHalfMax2+1
	enddo
	do while (histo(xHalfMax1).lt.maximum/2)
	  xHalfMax1 = xHalfMax1-1
	enddo

	!do the interpolation and get the equation of a line between the two points
	if ((xHalfMin1.eq.xHalfMin2).or.&
	    ((histo(xHalfMin1)).eq.(histo(xHalfMin2)))) then
		 xHalfMin1 = xHalfMin1 -1
	endif
	if ((xHalfMax1.eq.xHalfMax2).or.&
	   ((histo(xHalfMax1)).eq.(histo(xHalfMax2)))) then
		 xHalfMax1 = xHalfMax1 +1
	endif
	mLeft = (histo(xHalfMin1) - histo(xHalfMin2))/(xHalfMin1 - xHalfMin2)
	mRight = (histo(xHalfMax1) - histo(xHalfMax2))/(xHalfMax1 - xHalfMax2)
	bLeft = histo(xHalfMin1) - (mLeft*xHalfMin1)
	bRight = histo(xHalfMax1) - (mRight*xHalfMax1)
	!interpolated values are
	xLeft = (maximum/2 - bLeft)/mLeft;
	xRight = (maximum/2 - bRight)/mRight;

	fwhm = (xRight - xLeft)*binwidth;
! ***
! old version 
!	!find x positions of max/2
!	xHalfMin = 1 !smallest value possible
!	xHalfMax = numberBins -1 ! largest value possible
!	do k=1, maxIndex
!	  if (histo(k).le.maximum/2) then
!	    xHalfMin = k
!	  endif
!	enddo
!
!	do k=1, maxIndex
!	   if (histo(maxIndex+k).le.maximum/2) then
!	      xHalfMax = maxIndex+k
!	      exit
!	   endif
!	enddo
!!	do while (histo(xHalfMin).le.maximum/2) 
!!	  xHalfMin = xHalfMin +1
!!	enddo
!!	do while (histo(xHalfMax).le.maximum/2) 
!!	  xHalfMax = xHalfMax -1
!!	enddo
!	fwhm = (xHalfMax-xHalfMin)*binWidth

 	 return
	end subroutine calculateFWHM

!---------------------------------------------------------------------------------

!*********************************************************************************
!*********************************************************************************

      subroutine writeEmi(epsx1,epsy1,t1,pt1,epsx2,epsy2,t2,pt2,np1,np2,tauhat1,tauhat2)
      include 'collider_time_evolution_common.f90'
      real :: x1(nb),px1(nb),y1(nb),py1(nb),x2(nb),px2(nb),y2(nb),&
           py2(nb),t1(nb),t2(nb),pt1(nb),pt2(nb),dponp1,dponp2,&
           sigx1,sigx2,sigy1,sigy2,epsx1,epsy1,epsx2,&
           epsy2,sigt1,sigPt1,sigt2,sigPt2,epsl1,epsl2,&
	   fwhm1,fwhm2,tauhat1,tauhat2
      integer :: np1,np2,k

      sigt1=0
      sigPt1=0
      sigt2=0
      sigPt2=0

!-------by MICHAELA -----------------------
!	Calculate FWHM

!	to be able to use tauhat here I hat to change the argruments given to subroutine writeEmit
!	and add tauhat1 and tauhat2 to the list of arguments

	call calculateFWHM(t1,np1,fwhm1,tauhat1)
	call calculateFWHM(t2,np2,fwhm2,tauhat2)

!------------------------------------------


!     determine long. stand. dev.
!     beam 1
      do k=1,np1
         sigt1=sigt1+t1(k)**2
         sigPt1=sigPt1+pt1(k)**2
      enddo
      sigt1=sqrt(sigt1/np1)
      sigPt1=sqrt(sigPt1/np1)

      dponp1 = sigPt1/(beta*beta*gamma0) ! 1 sigma fractional momentum deviation     
      epsl1=sigt1*sigPt1*0.931494e9*Aatom*pi/Qatom ! long. emittance in eV s/charge at 1 sigma. not exact.

      do k=1,np2
         sigt2=sigt2+t2(k)**2
         sigPt2=sigPt2+pt2(k)**2
      enddo
      sigt2=sqrt(sigt2/np2)
      sigPt2=sqrt(sigPt2/np2)
      
      dponp2 = sigPt2/(beta*beta*gamma0) ! 1 sigma fractional momentum deviation        
      epsl2=sigt2*sigPt2*0.931494e9*Aatom*pi/Qatom ! long. emittance in eV s/charge at 1 sigma. not exact.

!     write to file
!      open(unit=33,file='emit.out',access='append')
!      write(33,*) epsx1,epsy1,epsl1,epsx2,epsy2,epsl2,kturns
!      close(33)

      open(unit=80,file='emittance.out',access='append')
      !write(80,'(I10,11E17.7)')kturns,real(kturns)/real(nturns)*eqTime,epsx1,epsy1,epsl1,sigt1,dponp1,epsx2,epsy2,epsl2,sigt2,dponp2
      write(80,'(I10,13E17.7)')kturns,real(kturns)/real(nturns)*eqTime,epsx1,epsy1,epsl1,sigt1,fwhm1,dponp1,&
                               &epsx2,epsy2,epsl2,sigt2,fwhm2,dponp2
      close(80)  

!     write screen output
!      write(6,*) 'Turn ',kturns
      write(6,'(3I10,8E15.5)') kturns,np1,np2,epsx1,epsy1,sigt1,fwhm1,epsx2,epsy2,sigt2,fwhm2
!      write(6,*) '------------------------------------------'

      return
      end

!*********************************************************************************
!*********************************************************************************

      subroutine mountTransv(avglinex,avgliney,emix,emiy,beamnr)
! stores transverse bunch profiles
      include 'collider_time_evolution_common.f90'

      real :: avglinex(nb),avgliney(nb),dx,dy,emix,emiy,xmn,ymn,sigx,sigy,xc,yc,ex,ey
      integer :: k,beamnr,ntot
      character(5) :: beam
      character(50) :: fname

      xmn=-5*sqrt(betax*emix)
      ymn=-5*sqrt(betay*emiy)

      dx = -2*xmn/nBins
      dy = -2*ymn/nBins

!      sigx=0
!      sigy=0
      ntot=0

!     construct filename
      write(beam,'(i0)') beamnr
      fname=trim('mountain-xy-B'//trim(beam)//'.out')
      open(unit=33,file=fname,access='append')

      do k=1,nBins
         xc=xmn+k*dx
         yc=ymn+k*dy
         write(33,100)xc,yc,avglinex(k),avgliney(k),float(kturns)
!         sigx=sigx+avglinex(k)*xc**2
!         sigy=sigy+avgliney(k)*yc**2
         ntot=ntot+avglinex(k)
!     zero it out
         avglinex(k)=0
         avgliney(k)=0
      enddo
      write(33,*)
 100  format(1p10e14.4)
      close(33)

!      sigx=sqrt(sigx/real(ntot))
!      sigy=sqrt(sigy/real(ntot))
      
      ex=sigx**2/betax
      ey=sigy**2/betay
!      fname=trim('emi-transv-B'//trim(beam)//'.out')
!      open(unit=33,file=fname,access='append')
!      write(33,*) ex,ey,kturns
!      close(33)
      
      return
      end

!*********************************************************************************
!*********************************************************************************

      real function I0(X)
!     
!       =========================================================
!       Compute modified Bessel function I0(x)
!       modified to use real*4 instead of real*8
!       using real*8 causes many strange numerical instabilities

        include 'collider_time_evolution_common.f90'
        REAL EL,X,BI0,X2,R,CA,XR,A(12),B(12),A1(8)
        INTEGER K,K0
!        DIMENSION A(12),B(12),A1(8)
!        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        X2=X*X
        IF (X.EQ.0.0) THEN
           BI0=1.0
           RETURN
        ELSE IF (X.LE.18.0) THEN
           BI0=1.0
           R=1.0
           DO 15 K=1,50
              R=0.25*R*X2/(K*K)
              BI0=BI0+R
              IF (ABS(R/BI0).LT.1.0E-15) GO TO 20
15         CONTINUE
20         CONTINUE
           R=1.0

        ELSE
           DATA A/0.125,7.03125E-2,&
                 7.32421875E-2,1.1215209960938E-1,&
                 2.2710800170898E-1,5.7250142097473E-1,&
                 1.7277275025845E0,6.0740420012735E0,&
                 2.4380529699556E01,1.1001714026925E02,&
                 5.5133589612202E02,3.0380905109224E03/
           DATA B/-0.375E0,-1.171875E-1,&
                 -1.025390625E-1,-1.4419555664063E-1,&
                 -2.7757644653320E-1,-6.7659258842468E-1,&
                 -1.9935317337513E0,-6.8839142681099E0,&
                 -2.7248827311269E01,-1.2159789187654E02,&
                 -6.0384407670507E02,-3.3022722944809E03/
           K0=12
           IF (X.GE.35.0) K0=9
           IF (X.GE.50.0) K0=7
           CA=EXP(X)/SQRT(2.0E0*PI*X)
           BI0=1.0E0
           XR=1.0E0/X
           DO 35 K=1,K0
35            BI0=BI0+A(K)*XR**K
           BI0=CA*BI0

        ENDIF
        I0=BI0
        RETURN
        END



!*********************************************************************************
!*********************************************************************************

    subroutine collimation(np,y,py,t,pt,x,px,np0,pnumber,nLostMom,nLostBeta)
      include 'collider_time_evolution_common.f90'
      real :: y(nb),py(nb),t(nb),pt(nb),x(nb),px(nb),pnumber,Jx,Jy,xDx,xmax2,ymax2,betaRatio
      real :: xCutBeta,yCutBeta,xCutBeta2,yCutBeta2,deltaPoP,xCutMom,xCutMom2,xi,pxi,yi,pyi
      integer :: iwrite,np,np0,i,j,koo,lost,nLostMom,nLostBeta

      ! loop through particles. For each particle, calculate betatron amplitude and momentum amplitude. loop through aperture cuts and check if outside.
      ! amplitude can be calcualted either taking into account both amplitude and phase, or amplitude only.

      koo=0

! cutoff values for betatron collimation
      xCutBeta = nSigCutBeta*sqrt(betax*refEmxy) ! sufficient to use beta function of ring
      xCutBeta2 = xCutBeta**2
      yCutBeta = nSigCutBeta*sqrt(betay*refEmxy)
      yCutBeta2=yCutBeta**2

      betaRatio=sqrt(betaxMom/betax)

! cutoff values for momentum collimation. 
      xCutMom = nSigCutMom * sqrt(refEmxy*betaxMom)

!      write(6,*) 'xcutbet=',xCutBeta,'ycutbet=',yCutBeta,'xcutmom=',xCutMom

      do i=1,np
         lost=0
         deltaPoP = pt(i)/(beta*beta*gamma0)
         if (collimAvgSwitch.eq.1) then ! max amplitude, disregarding phase (if 1 sim turn is many machine turns)
            xmax2 = x(i)**2 + px(i)**2 
            ymax2 = y(i)**2 + py(i)**2 
            xDx = sqrt(xmax2)*betaRatio + abs(dispxMom*deltaPoP)
         else ! take into account phase 
            xmax2 = x(i)**2 
            ymax2 = y(i)**2 
            xDx = abs( x(i)*betaRatio + dispxMom*deltaPoP )
         endif

         ! check hor. betatron cut (checking squares to save sqrt operation)
         if (xmax2.ge.xCutBeta2) then
            lost=1
            nLostBeta=nLostBeta+1
         endif

         ! check ver. betatron cut
         if ((ymax2.ge.yCutBeta2).and.(lost.eq.0)) then ! don't double count the particle if lost in both x and y
            lost=1
            nLostBeta=nLostBeta+1
         endif

         ! check momentum collimator
         if ((xDx.ge.xCutMom).and.(lost.eq.0)) then
            lost=1
            nLostMom=nLostMom+1
         endif
         
         ! if (lost.eq.1) then
         !    write(93,*) i,y(i),py(i),x(i),px(i),pt(i)
         ! endif

         if (lost.eq.0) then ! keep particle
            koo = koo + 1
            t(koo)=t(i)
            pt(koo)=pt(i)
            x(koo)=x(i)
            y(koo)=y(i)
            px(koo)=px(i)
            py(koo)=py(i)
         endif
       
      enddo
      np=koo
      return
    end subroutine collimation

!*********************************************************************************
!*********************************************************************************

!     Collisions with hour glass effect and crossing angle, including phase averaging.
!     assuming beams with Gaussian transverse densities and longitudinal
!     density given by the simulation. Calculating collision probability for a single 
!     particle as the product of the reduction factor (calculated once and approximated
!     to be the same for all particles) and the transverse density in 
!     action variable J only, averaged over phase, of the opposing beam and a scaling 
!     factor (containing cross section etc).


subroutine collision6a(np1,x1,px1,t1,pt1,y1,py1,np10,pnumber1,epsx1,epsy1,&
                       np2,x2,px2,t2,pt2,y2,py2,np20,pnumber2,epsx2,epsy2,&
                       nLostLum1,nLostLum2,iwrite,lumi,redfac,hglassfacOld)
!     sigI is the interaction cross section in barn
!     betaS is beta* at the IP in m
!     nlost are the number of lost particles in beam 1 and 2 - used for luminosity calculation

      include 'collider_time_evolution_common.f90'
      external ran3

      real :: x1(nb),px1(nb),t1(nb),pt1(nb),y1(nb),py1(nb),pnumber1 ! variables for beam1
      real :: x2(nb),px2(nb),t2(nb),pt2(nb),y2(nb),py2(nb),pnumber2 ! variables for beam2
      real :: xk,yk,z1,z2,sigx1,sigx2,sigy1,sigy2,dzsamp,probZ,&
           prob,sigm2,ran3,lum,dtsamp,thib2,lumInf,psi,&
           test,probavg1,probavg2,denlon1(nb),denlon2(nb),tk,lumprefac,&
           hglassfac(10),lumProb,lumIP1,probavgIP1,lumprefac1,covxy2,&
           Jx,Jy,epsx1,epsy1,epsx2,epsy2,argx,argy,i0x,i0y,I0,covxy1,lumi,redfac,&
           lumitest,hglassfacOld(10)
      integer :: k,m,k1,k2,nLostLum1,nLostLum2,nlost2,np1,np2,np10,np20,nkt,m2,iwrite

!     central particle for debugging - REMOVE FOR REAL RUNS
!$$$      x1(1)=0.0
!$$$      y1(1)=0.0
!$$$      px1(1)=0.0
!$$$      py1(1)=0.0
!----------------------------------------------------------


      sigm2=sigI*1.e-28   ! convert cross section from barn to m^2

      covxy1=0
      covxy2=0

      probavg1=0
      probavg2=0
      probavgIP1=0

      dtsamp = thib/longIntBins
      dzsamp = clight*dtsamp
      thib2 = thib/2

      do k=1,longIntBins
         denlon1(k)=0
         denlon2(k)=0
      enddo

!     determine transv. stand. dev. and long. distribution
!----------------------------------------------------------
!     beam 1
      do k=1,np1
!         covxy1=covxy1+x1(k)*y1(k)
         tk = (t1(k)+thib2)/dtsamp
         nkt = nint(tk)
         if(nkt.le.0) nkt=1
         if(nkt.gt.longIntBins) nkt=longIntBins
         denlon1(nkt)=denlon1(nkt)+1 ! binning in t, but same binning valid in z with dzsamp=clight*dtsamp
      enddo
      sigx1=sqrt(epsx1*betax)     ! beam size in the ring with average beta
      sigy1=sqrt(epsy1*betay)
      ! covxy1=covxy1/np1
      ! epsx1=sigx1**2/betax      ! geometric emittance
      ! epsy1=sigy1**2/betay

!     same for beam2
      do k=2,np2
!         covxy2=covxy2+x2(k)*y2(k)
         tk = (t2(k)+thib2)/dtsamp
         nkt = nint(tk)
         if(nkt.le.0) nkt=1
         if(nkt.gt.longIntBins) nkt=longIntBins
         denlon2(nkt)=denlon2(nkt)+1 
      enddo
      sigx2=sqrt(epsx2*betax)
      sigy2=sqrt(epsy2*betay)
      ! covxy2=covxy2/np2
      ! epsx2=sigx2**2/betax
      ! epsy2=sigy2**2/betay

!      check if levelling is required, and change betaS accordingly
!------------------------------------------------------------------------------------------------------------
!******************* MSCHAUMA ******************************

      if(levellingSwitch.eq.1) then
         do m2=1,nIPs
           if(levIPSwitch(m2).eq.1) then
             ! claculate potential luminosity
             lumitest=hglassfacOld(m2)*(np1*pnumber1/np10)*(np2*pnumber2/np20)* &
                  nbunches*(clight/circ)/(2*pi*sqrt(betaSMin(m2)**2/(betax*betay)* &
                 (sigx1**2+sigx2**2)*(sigy1**2+sigy2**2)))* 1.e-4             ! scaling factor 1e-4 to get unit cm^-2 s^-1   
 !            write(6,*) 'betaS=',betaS(m2), '   betaSMin=', betaSMin(m2)
             !lumimax is the maximal value of the luminosity allowed - value to which lumi should be levelled        
             if((betaS(m2).ge.betaSMin(m2)).and.(lumitest.gt.lumimax(m2))) then
               betaS(m2) = hglassfacOld(m2)*(np1*pnumber1/np10)*(np2*pnumber2/np20)* &
                      nbunches*(clight/circ)/(2*pi*lumimax(m2)*sqrt(1/(betax*betay)* &
                      (sigx1**2+sigx2**2)*(sigy1**2+sigy2**2)))* 1.e-4 
                 if (iwrite.eq.1) then 
                   write(6,*) 'betaS=',betaS(m2) , '  IP=', m2
                 endif
!             else
!               betaS(m2) = betaSMin
             endif
           endif

         enddo
      endif

!     calculate reduction factors at all IPs. To speed up, use the same factor for all particles
!------------------------------------------------------------------------------------------------------------
!     assume the same beta* in x and y

      do m2=1,nIPs
         hglassfac(m2)=0
      enddo

      do k=1,longIntBins
         z1=((real(k)/longIntBins)*thib-thib2)*clight
         do m=1,longIntBins
            z2=-((real(m)/longIntBins)*thib-thib2)*clight
            do m2=1,nIPs
               hglassfac(m2)=hglassfac(m2)+  &
                    exp(-(2*(z1+z2)**2*betaS(m2)*sin(theta(m2))**2)/ & ! theta is half of the crossing angle
                    ((epsx1+epsx2)*((z1+z2)**2+2*betaS(m2)**2*(1+cos(2*theta(m2))))))* &
                    denlon1(k)/np1*denlon2(m)/np2/(1+((z1+z2)/(2*betaS(m2)*cos(theta(m2))))**2)
            enddo
         enddo
      enddo
 
      lumprefac=0
      lumprefac1=0
      do m2=1,nIPs              ! sum the effective luminosity pre-factor, which is IP-dependent
         lumprefac=lumprefac+IPmult(m2)*hglassfac(m2)/betaS(m2) !not the "real" hour glass factor, contains also 1/beta*
      enddo

      lumprefac=lumprefac*sqrt(betax*betay)*sigm2*timeRatio ! IP-independent pre-factor

      if (iwrite.eq.1) lumprefac1=hglassfac(1)/betaS(1)*sqrt(betax*betay)*sigm2*timeRatio ! prefactor for probability at IP1 only


!     loop through particles, sample interaction according to probability, beam 1
!--------------------------------------------------------------------------------------
!     beam 1

      k1=0 ! counting number of surviving particles
      do k=1,np1
         Jx=(x1(k)**2+px1(k)**2)/(2*betax) ! action variable. x1 given in ring so need to use beta of the ring 
         Jy=(y1(k)**2+py1(k)**2)/(2*betay)
         argx=Jx/(2*epsx2)      ! argument for the density function
         argy=Jy/(2*epsy2)
         prob=lumprefac*np2*exp(-argx)*I0(argx)*exp(-argy)*I0(argy)* & !interaction probability averaged over phase
              (pnumber2/np20)/ & ! scaling factor for the number of real particles in beam 2
              (2*pi*sigx2*sigy2)
         
         if (prob.ge.1) then
            write(6,*) 'prob=',prob
            write(6,*) 'x=',x1(k),'y=',y1(k),'px=',px1(k),'py=',py1(k)
            write(6,*) 'Jx=',Jx,'Jy=',Jy,'argx=',argx,'argy=',argy
            write(6,*) 'turn',kturns
            write(6,*) 'k=',k
            write(6,*) 
         endif            

         if (iwrite.eq.1) then  ! store average probabilities for output
            probavg1=probavg1+prob 
            probavgIP1=probavgIP1 + lumprefac1*np2 * exp(-argx)*I0(argx) * exp(-argy)*I0(argy)*(pnumber2/np20) / (2*pi*sigx2*sigy2)
         endif

         test=ran3(iseed)
                   
         if (test.ge.prob) then ! keep particle
            k1=k1+1             ! count number of surviving macro particles
            t1(k1)=t1(k)
            pt1(k1)=pt1(k)
            x1(k1)=x1(k)
            px1(k1)=px1(k)
            y1(k1)=y1(k)
            py1(k1)=py1(k)            
         endif
      enddo
      nLostLum1=nLostLum1+np1-k1   ! number of lost particles from beam 1  = number of real particles lost from both bunches during a single bunch crossing



!     loop through particles, sample interaction according to probability, beam 2     
!-------------------------------------------------------------------------------------- 
!     beam 2
      k2=0                      ! counting the number of surviving particles
      
      x2(1)=0
      y2(1)=0
      
      do k=1,np2

         Jx=(x2(k)**2+px2(k)**2)/(2*betax) ! action variable. x1 given in ring so need to use beta of the ring 
         Jy=(y2(k)**2+py2(k)**2)/(2*betay)
         argx=Jx/(2*epsx1)        ! argument for the density function
         argy=Jy/(2*epsy1)
         prob=lumprefac*np1*exp(-argx)*I0(argx)*exp(-argy)*I0(argy)* &
              (pnumber1/np10)/ &  ! scaling factor for the number of real particles in beam 2
             (2*pi*sigx1*sigy1)

         probavg2=probavg2+prob ! used for debugging
         test=ran3(iseed)
         
         if (test.ge.prob) then ! keep particle
            k2=k2+1
            t2(k2)=t2(k)
            pt2(k2)=pt2(k)
            x2(k2)=x2(k)
            px2(k2)=px2(k)
            y2(k2)=y2(k)
            py2(k2)=py2(k)       
         endif
      enddo
      nLostLum2=nLostLum2+np2-k2 ! n.o. particles lost from bunch 2 from luminosity
      

!     write data to file: calculated luminosity and reduction factor  for first IP
!------------------------------------------------------------------------------------------------------------
      if (iwrite.eq.1) then

         probavg1=probavg1/np1
         probavg2=probavg2/np2
         probavgIP1=probavgIP1/np1
         
         lum=hglassfac(1)*(np1*pnumber1/np10)*(np2*pnumber2/np20)* &
             nbunches*(clight/circ)/(2*pi*sqrt(betaS(1)**2/(betax*betay)* &
             (sigx1**2+sigx2**2)*(sigy1**2+sigy2**2)))* 1.e-4             ! scaling factor 1e-4 to get unit cm^-2 s^-1

!        calculate luminosity from actual number of particles lost in simulation
         lumInf=nbunches*(clight/circ)*real(nLostLum1)/nwrite/ & 
             (sigI*1.e-24)*  &  ! scaling factor to get unit cm^-2 s^-1
             (pnumber1/np10)/timeRatio ! scal. fact. to account for smaller n.o. macro part. and different time scale

!        calcualte luminosity from average interaction probability  - this is the number of particles from ALL IPs!!
         lumProb=nbunches*(clight/circ)*probavg1*np1/ &
             (sigI*1.e-24)*    &! scaling factor to get unit cm^-2 s^-1
             (pnumber1/np10)/timeRatio ! scal. fact. to account for smaller n.o. macro part. and different time scale

!        calcualte luminosity from average interaction probability at IP1 only
         lumIP1=nbunches*(clight/circ)*probavgIP1*np1/ &
             (sigI*1.e-24)*    & ! scaling factor to get unit cm^-2 s^-1
             (pnumber1/np10)/timeRatio ! scal. fact. to account for smaller n.o. macro part. and different time scale

!        calcualte the tilt angle of the x-y ellipse due to linear betatron coupling
!         psi=atan(2*covxy1/(sigx1**2-sigy1**2))/2
         
!         write(6,*) 'luminosity=',lumIP1,'cm^-2 s^-1, ','geom. reduction factor=',hglassfac(1)
!     &        ,'lum calc.=',lum
!         write(6,*) 'probavg1,probavg2=',probavg1,probavg2         
!        nLostLum1 contains the total number of particles lost from burnoff over nwrite turns
!         open(unit=52,file='nLostBurnoff.out',access='append')
!         write(52,*) kturns,nLostLum1,lumInf,lum,lumProb,lumIP1,
!     &        hglassfac(1),sigx1,sigy1,covxy1,sigx2,sigy2,covxy2,psi 
!         close(52)

     !     open(unit=80,file='luminosity.out',access='append')
     !     write(80,'(I10,2E17.7,I10,2E17.7)') kturns,real(kturns)/real(nturns)*eqTime,lumIP1,&
     ! &        nLostLum1,real(nLostLum1)/real(np10)*pnumber1,hglassfac(1)
     !     close(80)
         
!         nLostLum1=0

         
      endif
      lumi=lumIP1
      redfac=hglassfac(1)
      hglassfacOld = hglassfac
      
      np1=k1
      np2=k2

      return
      end


!*********************************************************************************
!*********************************************************************************

!     Collisions without assumptions on transverse or longitudinal distribution.
!     Including hour glass effect and crossing angle.
!     Transverse and longitudinal density taken from discrete binning of distributions.
!     Calculating collision probability for a single 
!     particle by propagating it along s, using x=x0+s*x0' in the rotated coordinate 
!     system of the opposing bunch. For each bin in s, 
!     the x and y-values are determined and the collision probability is given by the 
!     density of the opposing beam at that particular s.
!     Assuming that the x- and y-distributions are independent of each other,
!     that is f(x,y,s)=f1(x,s)*f2(y,s)*f3(s). 
!     f3(s) gives the longitudinal distribution, while the only s-dependence in f1 and f2 is 
!     a scaling of the width of distribution at s=0 by sqrt[beta(s)/beta(0)] (see paper).

!     This method is very slow, but does not make any assumptions on the beam distributions


      subroutine collision1d(np1,x1,px1,t1,pt1,y1,py1,np10,pnumber1,np2,x2,px2,t2,pt2,y2,py2,np20,pnumber2,nLostLum1,&
           nLostLum2,iwrite,lumi,redfac,emixy)
!     as input argument to make it easy to switch between different collision routines.
!     sigI is the interaction cross section in barn
!     betaSt is beta* at the IP in m
!     nlost are the number of lost particles in beam 1 and 2 - used for luminosity calculation

      include 'collider_time_evolution_common.f90'
      external ran3

      real :: x1(nb),px1(nb),t1(nb),pt1(nb),y1(nb),py1(nb),pnumber1 ! variables for beam1
      real :: x2(nb),px2(nb),t2(nb),pt2(nb),y2(nb),py2(nb),pnumber2 ! variables for beam2

      real :: xk,yk,z1,z2,sigx1,sigx2,sigy1,sigy2,dzsamp,probZ,ddx,&
          prob(10),sigm2,ran3,lum,dtsamp,thib2,lumInf,bx,by,&
          test,probavg1,probavg2,denlon1(nb),denlon2(nb),denx1(nb),&
          denx2(nb),deny1(nb),deny2(nb),tk,lumprefac(10),s,&
          hglassfac(10),lumProb,lumIP1,probavgIP1,lumprefac1,&
          xIP(10),yIP(10),xPrime(10),yPrime(10),ampx,emixy,probTot,&
          kappa,ratX(10),ratY(10),dxIP(10),dyIP(10),x00,x0(10),x0P,&
          Cos2Th(10),Sin2Th(10),CosTh,SinTh,IPprefac(10),z0,&
          epsx1,epsx2,denom(10),y0(10),y0p(10),lumNom,hglassIP1,&
          epsy1,epsy2,lumi,redfac
      integer :: k,m,k1,k2,nLostLum1,nLostLum2,nlost2,np1,np2,np10,np20,nkt,m2,iwrite,nky,nkx

!     central particle for debugging - REMOVE FOR REAL RUNS
!      x1(1)=0.0
!      y1(1)=0.0
!      px1(1)=0.0
!      py1(1)=0.0
!      t1(1)=0.0
!----------------------------------------------------------

!$$$      do m2=1,nIPs
!$$$         betaSX(m2)=betaS(m2)
!$$$         betaSY(m2)=betaS(m2)
!$$$      enddo


!$$$      write(6,*)
!$$$      write(6,*) 'beginning of routine'
!$$$      write(6,*) 'x1=',x1(1),'px1=',px1(1)
!$$$      write(6,*) 'y1=',y1(1),'py1=',py1(1)
!$$$      write(6,*)

!     on odd turns, if angleSwitch==1 => swap x and y, and swap back after the collision
!     this way, the crossing angle switches between hor. and ver. plane in order to 
!     simulate for example the angles in the different planes at IP1 and IP5 in the LHC.
!     this can be used to avoid to introduce assymetries when only collisions are active

      if ((mod(kturns,2).eq.1).and.(angleSwitch.eq.1)) then 
!         write(6,*)'odd turn'
         bx = circ/(2*pi*tuney) ! beta functions in the ring
         by = circ/(2*pi*tunex)
         do k=1,np1
            xk=x1(k)
            x1(k)=y1(k)
            y1(k)=xk

            xk=px1(k)
            px1(k)=py1(k)
            py1(k)=xk
         enddo
         do k=1,np2
            xk=x2(k)
            x2(k)=y2(k)
            y2(k)=xk

            xk=px2(k)
            px2(k)=py2(k)
            py2(k)=xk
         enddo
      else
!         write(6,*) 'even turn'
         by = circ/(2*pi*tuney) ! beta functions in the ring
         bx = circ/(2*pi*tunex)
      endif

!$$$      write(6,*)
!$$$      write(6,*) 'after switch'
!$$$      write(6,*) 'x1=',x1(1),'px1=',px1(1)
!$$$      write(6,*) 'y1=',y1(1),'py1=',py1(1)
!$$$      write(6,*)


      prob(1)=0
      sigm2=sigI*1.e-28         ! convert cross section from barn to m^2


      probavg1=0
      probavg2=0
      probavgIP1=0

!     binning in x and y between +-5*ampx and in t between +-thib2 in longIntBins number of bins
!     determine step sizes
      ampx = sqrt(bx*emixy)
      ddx=(10*ampx)/longIntbins
      dtsamp = thib/longIntBins
      dzsamp = clight*dtsamp
      thib2 = thib/2

      do k=1,longIntBins
         denlon1(k)=0
         denlon2(k)=0
         denx1(k)=0
         denx2(k)=0
         deny1(k)=0
         deny2(k)=0
      enddo

!     determine transv. and long. distribution
!----------------------------------------------------------
!     note that the binning is done in x and px, that is 
!     the coordinates in the ring with the average beta.
!     the interaction probability is calculated using these
!     bins but scaled by a different bin width dxIP reflecting
!     the ratio to the beta at the IP
      
!     beam 1
      do k=1,np1
         xk = (x1(k)+5*ampx)/ddx
         nkx = nint(xk)  ! find nearest integer = bin number
         !nkx=floor(xk)+1
         if(nkx.le.0)  nkx=1
         if(nkx.gt.longIntBins) nkx=longIntBins
         denx1(nkx)=denx1(nkx)+1

         yk = (y1(k)+5*ampx)/ddx
         nky = nint(yk)  ! find nearest integer = bin number
         !nky=floor(yk)+1
         if(nky.le.0)  nky=1
         if(nky.gt.longIntBins) nky=longIntBins
         deny1(nky)=deny1(nky)+1

         tk = (t1(k)+thib2)/dtsamp
         nkt = nint(tk)
         if(nkt.le.0) nkt=1
         if(nkt.gt.longIntBins) nkt=longIntBins
         denlon1(nkt)=denlon1(nkt)+1 ! binning in t, but same binning valid in z with dzsamp=clight*dtsamp
      enddo

!     same for beam2
      do k=2,np2
         xk = (x2(k)+5*ampx)/ddx
         nkx = nint(xk)  ! find nearest integer = bin number
         !nkx=floor(xk)+1
         if(nkx.le.0)  nkx=1
         if(nkx.gt.longIntBins) nkx=longIntBins
         denx2(nkx)=denx2(nkx)+1

         yk = (y2(k)+5*ampx)/ddx
         nky = nint(yk)  ! find nearest integer = bin number
         !nky=floor(yk)+1
         if(nky.le.0)  nky=1
         if(nky.gt.longIntBins) nky=longIntBins
         deny2(nky)=deny2(nky)+1

         tk = (t2(k)+thib2)/dtsamp
         nkt = nint(tk)
         if(nkt.le.0) nkt=1
         if(nkt.gt.longIntBins) nkt=longIntBins
         denlon2(nkt)=denlon2(nkt)+1 
      enddo

!     export densities for debugging
!      do k=1,longIntBins
!         tk=(real(k)/longIntBins)*thib-thib2
!         xk=(real(k)/longIntBins)*10*ampx-5*ampx
!         yk=xk
!         write(10,*) xk,denx1(k)
!         write(11,*) yk,deny1(k)
!         write(12,*) tk,denlon1(k)
!         write(13,*) xk,denx2(k)
!         write(14,*) yk,deny2(k)
!         write(15,*) tk,denlon2(k)
!      enddo


      do m2=1,nIPs
         ratX(m2)=sqrt(betaS(m2)/bx)
         ratY(m2)=sqrt(betaS(m2)/by)
         dxIP(m2)=ddx*ratX(m2)
         dyIP(m2)=ddx*ratY(m2)
         Cos2Th(m2)=cos(2*theta(m2))
         Sin2Th(m2)=sin(2*theta(m2))
         IPprefac(m2)=cos(theta(m2))**2*(pnumber2/np20)*2* & ! common prefactor to P1 for all particles, beam 1
             sigm2*timeRatio/ (dxIP(m2)*dyIP(m2)*real(np2)**2)  *IPmult(m2)       ! multiplicity of IP
      enddo
         


!     loop through particles, sample interaction according to probability, beam 1
!--------------------------------------------------------------------------------------

      k1=0 ! counting number of surviving particles
      
      do k=1,np1
         
         z0=t1(k)*clight        ! long. distance from bunch center of the particle in beam 1
         
         do m2=1,nIPs           ! calcualte transverse angle at each of the IPs
!           x00,x0P=x',z0 are coordinates in system of beam1 at t=0. ratX scales from ring to IP.
!           assuming trajectory is x=x0+x0P(s-z0)=x00+x0P*s and s=c*t+z0 in the long. plane
            x0P=px1(k)/bx/ratX(m2) ! transforming back from normalized coord., assuming alpha=0 at all IPs
            x0(m2)=x1(k)*ratX(m2)
            x00=x0(m2)-z0*x0P

!           transform to the rotated system of the opposing beam
!           trajectory is x=xIP+xPrime*s2 in the rotated system
            denom(m2)=Cos2Th(m2)+x0P*Sin2Th(m2)
            xIP(m2)=x00/denom(m2)
            xPrime(m2)=(x0P*Cos2Th(m2)-Sin2Th(m2))/denom(m2)

!           vertical coordinates 
!           in the rotated system, the trajectory is y=yIP+yPrime*s2
            y0(m2)=y1(k)*ratY(m2)
            y0P(m2)=py1(k)/by/ratY(m2)
            yIP(m2)=y0(m2)-y0P(m2)*(Cos2Th(m2)*z0+Sin2Th(m2)*x0(m2))/denom(m2)
            yPrime(m2)=y0P(m2)/denom(m2)

            prob(m2)=0          ! zero out probability at each IP
         enddo
         do m=1, longIntBins    ! integrating over longitudinal coordinate
            z2=-((real(m)/longIntBins)*thib-thib2)*clight        ! longitudinal distance from center of bunch


            do m2=1,nIPs
               s=(z2*denom(m2)+z0*Cos2Th(m2)+x0(m2)*Sin2Th(m2))/ (1+denom(m2))   ! s is distance in fixed system, s=0 at IP
               kappa=sqrt(1+(s/betaS(m2))**2) ! scaling factor to account for distance from IP
               nkx=nint(((xIP(m2)+s*xPrime(m2))/kappa+5*ampx*ratX(m2))/dxIP(m2))  ! find x-bin at this s-position. 
               nky=nint(((yIP(m2)+s*yPrime(m2))/kappa+5*ampx*ratY(m2))/dyIP(m2))
               if(nky.le.0)  nky=1
               if(nky.gt.longIntBins) nky=longIntBins
               if(nkx.le.0)  nkx=1
               if(nkx.gt.longIntBins) nkx=longIntBins
               prob(m2)=prob(m2)+denx2(nkx)*deny2(nky)*denlon2(m)/(kappa*kappa) ! summing the collision probability at each IP
            enddo
         enddo
         
         probTot=0
         
         do m2=1,nIPs
            probTot=probTot+IPprefac(m2)*prob(m2)/(1+denom(m2))
         enddo      

         if (iwrite.eq.1) then  ! store average probabilities for output
            probavg1=probavg1+probTot ! average over all probability in beam1
            probavgIP1=probavgIP1+IPprefac(1)*prob(1)*denom(1)/(1+denom(1))
         endif

         test=ran3(iseed)
                   
         if (test.ge.probTot) then ! keep particle
            k1=k1+1
            t1(k1)=t1(k)
            pt1(k1)=pt1(k)
            x1(k1)=x1(k)
            px1(k1)=px1(k)
            y1(k1)=y1(k)
            py1(k1)=py1(k)            
         endif
      enddo
      nLostLum1=nLostLum1+np1-k1          ! number of lost particles from beam 1  = number of real particles lost from both bunches during a single bunch crossing

      do m2=1,nIPs
         IPprefac(m2)=cos(theta(m2))**2*(pnumber1/np10)*2* & ! common prefactor to P1 for all particles, beam 2
             sigm2*timeRatio/  &
             (dxIP(m2)*dyIP(m2)*real(np1)**2) &
             *IPmult(m2)       ! multiplicity of IP
      enddo
         

!     loop through particles, sample interaction according to probability, beam 2     
!-------------------------------------------------------------------------------------- 
!     beam 2
      k2=0                      ! counting the number of surviving particles
      
      do k=1,np2
         
         z0=t2(k)*clight        ! long. distance from bunch center of the particle in beam 2
         
         do m2=1,nIPs           ! calcualte transverse angle at each of the IPs
            x0P=px2(k)/bx/ratX(m2) ! transforming back from normalized coord., assuming alpha=0 at all IPs
            x0(m2)=x2(k)*ratX(m2)
            x00=x0(m2)-z0*x0P

            denom(m2)=Cos2Th(m2)+x0P*Sin2Th(m2)
            xIP(m2)=x00/denom(m2)
            xPrime(m2)=(x0P*Cos2Th(m2)-Sin2Th(m2))/denom(m2)

            y0(m2)=y2(k)*ratY(m2)
            y0P(m2)=py2(k)/by/ratY(m2)
            yIP(m2)=y0(m2)-y0P(m2)*(Cos2Th(m2)*z0+Sin2Th(m2)*x0(m2)) /denom(m2)
            yPrime(m2)=y0P(m2)/denom(m2)

            prob(m2)=0
         enddo
         do m=1, longIntBins    ! integrating over longitudinal coordinate




            z1=-((real(m)/longIntBins)*thib-thib2)*clight  ! longitudinal distance from center of bunch

            do m2=1,nIPs
               s=(z1*denom(m2)+z0*Cos2Th(m2)+x0(m2)*Sin2Th(m2))/ (1+denom(m2))   ! s is distance in fixed system, s=0 at IP
               kappa=sqrt(1+(s/betaS(m2))**2) ! scaling factor to account for distance from IP
               nkx=nint(((xIP(m2)+s*xPrime(m2))/kappa+5*ampx*ratX(m2)) /dxIP(m2))  ! find x-bin at this s-position. assuming x=x0+s*x0P at each IP
               nky=nint(((yIP(m2)+s*yPrime(m2))/kappa+5*ampx*ratY(m2)) /dyIP(m2))
               if(nky.le.0)  nky=1
               if(nky.gt.longIntBins) nky=longIntBins
               if(nkx.le.0)  nkx=1
               if(nkx.gt.longIntBins) nkx=longIntBins
               prob(m2)=prob(m2)+denx1(nkx)*deny1(nky)*denlon1(m)/(kappa*kappa) ! summing the collision probability at each IP
            enddo
         enddo
         
         probTot=0

         do m2=1,nIPs
            probTot=probTot+IPprefac(m2)*prob(m2)/(1+denom(m2))
         enddo      

         if (iwrite.eq.1) then  ! store average probabilities for output
            probavg2=probavg2+probTot ! average over all probability in beam1
         endif

         test=ran3(iseed)
                   
         if (test.ge.probTot) then ! keep particle
            k2=k2+1
            t2(k2)=t2(k)
            pt2(k2)=pt2(k)
            x2(k2)=x2(k)
            px2(k2)=px2(k)
            y2(k2)=y2(k)
            py2(k2)=py2(k)            
         endif
      enddo

      nLostLum2=nLostLum2+np2-k2          ! number of lost particles from beam 2

!$$$      write(6,*)
!$$$      write(6,*) 'after collision'
!$$$      write(6,*) 'x1=',x1(1),'px1=',px1(1)
!$$$      write(6,*) 'y1=',y1(1),'py1=',py1(1)
!$$$      write(6,*)

!     switch back x and y on odd turns
!---------------------------------------------------------------------------
      if ((mod(kturns,2).eq.1).and.(angleSwitch.eq.1)) then 
         do k=1,np1
            xk=x1(k)
            x1(k)=y1(k)
            y1(k)=xk

            xk=px1(k)
            px1(k)=py1(k)
            py1(k)=xk
         enddo
         do k=1,np2
            xk=x2(k)
            x2(k)=y2(k)
            y2(k)=xk

            xk=px2(k)
            px2(k)=py2(k)
            py2(k)=xk
         enddo
      endif

!$$$      write(6,*)
!$$$      write(6,*) 'after switched back'
!$$$      write(6,*) 'x1=',x1(1),'px1=',px1(1)
!$$$      write(6,*) 'y1=',y1(1),'py1=',py1(1)
!$$$      write(6,*)




!     write data to file: calculate luminosity and hourglass for first IP
!------------------------------------------------------------------------------------------------------------
      if (iwrite.eq.1) then

         probavg1=probavg1/np1
         probavg2=probavg2/np2
         probavgIP1=probavgIP1/np1

!        calculate the transverse sigmas
         sigx1=0
         sigx2=0
         sigy1=0
         sigy2=0
         do k=1,np1
            sigx1=sigx1+x1(k)**2
            sigy1=sigy1+y1(k)**2
         enddo
         do k=1,np2
            sigx2=sigx2+x2(k)**2
            sigy2=sigy2+y2(k)**2
         enddo

         sigx1=sqrt(sigx1/np1)  ! beam size in the ring with average beta
         sigy1=sqrt(sigy1/np1)
         sigx2=sqrt(sigx2/np2)  
         sigy2=sqrt(sigy2/np2)

         epsx1=sigx1**2/bx   ! geometric emittance
         epsx2=sigx2**2/bx
         epsy1=sigy1**2/by   ! geometric emittance
         epsy2=sigy2**2/by

!         write(17,*) kturns,epsx1,epsy1,epsx2,epsy2

!        calcualte luminosity at IP1, assuming Gaussian transverse distributions
!        first, calculate hourglass factor at IP1

         hglassfac(1)=0

!         write(6,*) 'emittances: ',sigx1,sigx2
!         stop
         CosTh=cos(theta(1))
         SinTh=sin(theta(1))
         do k=1,longIntBins
            z1=((real(k)/longIntBins)*thib-thib2)*clight
            do m=1,longIntBins  ! reduction factor assuming Gaussian distributions
               z2=-((real(m)/longIntBins)*thib-thib2)*clight
               hglassfac(1)=hglassfac(1)+ exp(-(2*(z1+z2)**2*betaS(1)*SinTh**2)/ &! theta is half of the crossing angle
                   ((epsx1+epsx2)*((z1+z2)**2+2*betaS(1)**2*(1+Cos2Th(1)))))*&
                   denlon1(k)/np1*denlon2(m)/np2/ (1+((z1+z2)/(2*betaS(1)*CosTh))**2)
            enddo
         enddo

 
!        luminosity at IP1
         lum=hglassfac(1)*(real(np1)*pnumber1/real(np10))*(real(np2)*pnumber2/real(np20))*&
             nbunches*(clight/circ)/(2*pi*sqrt(betaS(1)**2/(bx*by)*(sigx1**2+sigx2**2)*(sigy1**2+sigy2**2)))*1.e-4 ! 1.e-4 to get unit cm^-2 s^-1

!        calculate luminosity from actual number of particles lost in simulation
         lumInf=nbunches*(clight/circ)*real(nLostLum1)/nwrite/  &
             (sigI*1.e-24)   &  ! scaling factor to get unit cm^-2 s^-1
             *(pnumber1/np10)/timeRatio ! scal. fact. to account for smaller n.o. macro part. and different time scale

!        calcualte luminosity from average interaction probability  - this is the number of particles from ALL IPs!!
         lumProb=nbunches*(clight/circ)*probavg1*np1/ &
             (sigI*1.e-24)  &   ! scaling factor to get unit cm^-2 s^-1
             *(pnumber1/np10)/timeRatio ! scal. fact. to account for smaller n.o. macro part. and different time scale

!        calcualte luminosity from average interaction probability at IP1 only. compensate for multiplicity of IP
         lumIP1=nbunches*(clight/circ)*probavgIP1*np1/IPmult(1)/ &
             (sigI*1.e-24)    & ! scaling factor to get unit cm^-2 s^-1
             *(pnumber1/np10)/timeRatio ! scal. fact. to account for smaller n.o. macro part. and different time scale

!        use lumIP1 to get alternative estimate of the reduction factor
         lumNom=(real(np1)*pnumber1/real(np10))*(real(np2)*pnumber2/real(np20))* &
             nbunches*(clight/circ)/(2*pi*sqrt(betaS(1)**2/(bx*by)*(sigx1**2+sigx2**2)*(sigy1**2+sigy2**2)))*1.e-4
         hglassIP1=lumIP1/lumNom


!$$$         write(6,*) 'collision1d: '
!$$$         write(6,*) 'lum. reduction factor Gauss=',
!$$$     &        hglassfac(1),', lum. red. fac. from prob=',
!$$$     &        hglassIP1,
!$$$     &        ', lum from avg prob=',lumIP1,
!$$$     &        ', lum calc.=',lum
!$$$         write(6,*) 'probavg1,probavg2=',probavg1,probavg2     
!$$$
!$$$         open(unit=52,file='nLostBurnoff.out',access='append')
!$$$         write(52,*) kturns,nLostLum1,lumInf,lum,lumProb,lumIP1,
!$$$     &        hglassfac(1),hglassIP1 ! nLostLum1 contains the total number of particles lost from burnoff over nwrite turns
!$$$         close(52)


         ! write(6,*) 'luminosity=',lumIP1,'cm^-2 s^-1, ', 'geom. reduction factor=',hglassfac(1)

         ! open(unit=80,file='luminosity.out',access='append')
         ! write(80,'(I10,2E17.7,I10,2E17.7)') kturns,real(kturns)/real(nturns)*eqTime,lumIP1,&
         !      nLostLum1,real(nLostLum1)/real(np10)*pnumber1,hglassfac(1)
         ! close(80)
         
         ! nLostLum1=0
         
         
         
      endif
      lumi=lumIP1
      redfac=hglassfac(1)

      np1=k1
      np2=k2

      return
      end


!********************************************************************
!********************************************************************
!   MSCHAUMA    
!********************************************************************
!     Subroutines used for Stochastic Cooling only
!********************************************************************

!*********************************************************************************
!*********************************************************************************


      subroutine zerooutarrays(denlono,denlonoo,dpwakenotch,dpwakenotch2,dpytot,dpxtot,dpalmer)
      include 'collider_time_evolution_common.f90'


      real :: denlono(nb),denlonoo(nb),dpwakenotch(nb),dpwakenotch2(nb),dpytot(nb),dpxtot(nb),dpalmer(nb)
      integer :: k

! declare dpalmer
! zero out memory arrays
      do k=1,nresamp
         denlono(k)=0
         denlonoo(k)=0
         dpwakenotch(k)=0
         dpwakenotch2(k)=0
	 dpytot(k)=0
	 dpxtot(k)=0
         dpalmer(k)=0
      enddo

      return
      end subroutine zerooutarrays


!**********************************************************************
!**********************************************************************

      SUBROUTINE SORT2(N,RA,RB)
      DIMENSION RA(N)
! NEXT DECLARATIONS FOR INTEGER TAG-ALONG ARRAY
! BASED ON QUICKSORT
      INTEGER RB(N),RRB
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB=RB(L)
        ELSE
          RRA=RA(IR)
          RRB=RB(IR)
          RA(IR)=RA(1)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(I)=RRB
      GO TO 10
      END


!**********************************************************************
!**********************************************************************


      SUBROUTINE FOUR1(DATA,NN,ISIGN)
! COMMENTS BY MIKE BLASKIEWICZ
! DOES AN FFT 
!  DATA(2*K-1), K=1,..NN, ARE REAL PARTS OF INPUT ARRAY
!  DATA(2*K)  , K=1,..NN, ARE IMAGINARY PARTS OF INPUT ARRAY
! DOES FFT IN PLACE
! DPHI = ISIGN*(2*PI)/NN
! FT(J)  = SUM_{K=1}^NN  (D(2*K-1) + (0,1)*D(2*K))*EXP((0,1)*DPHI*(K-1)*(J-1))
! ON OUTPUT D(2*K-1) = REAL(FT(K)), D(2*K) = IMAG(FT(K))


      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END
!**********************************************************************
!**********************************************************************
  
     subroutine setupres(plane,np,x,t,tlob,coffset,dpalmer,denlono,denlonoo,dptot)
        include 'collider_time_evolution_common.f90'

	real :: tlob,dpalmer(nb),denlono(nb),denlonoo(nb),t(nb),x(nb),dptot(nb),coffset
	character(4) :: plane
	real :: dtsamp,tlobHalf,tk,sk,mm,fhi,flo,xk
	integer :: k, nlo,nhi,np
	

! calculates line density and dipole driving terms for resonator
! does binning
! also calculates dipole angular density
!  --- tlob is the total (4 sigma) bunch length in ns
!  --- nresamp is the number of samples used for ffts
!  --- dtsamp length of longitudinal bins for ffts
      dtsamp = tlob/nresamp
! zero out the sums
      if(plane.eq.'long') then
        do k=1,nresamp
!         denlonoo is data from previous pickup, 2 turn needs two more arrays
!         or rewrite with only one longitudinal setup call per turn.
           denlono(k)=denlonoo(k)
           denlonoo(k)=dpalmer(k)
           dpalmer(k)=0
        enddo
      elseif(plane.eq.'tran')then
         do k=1,nresamp
            dptot(k)=0
         enddo
      else
	write(6,*)' Unknown plane-settings in routine SETUPRES  '
      endif

!  k = 1 -> 0+dtsamp, k=nresamp -> tlob

! use linear interpolation
! bunch should be between zero and
      tlobHalf = tlob/2 

! dispu is dispersion at pickup
!      coeffpal = dispu/gamma0
      do k=1,np
         tk = t(k)
! 1.le.sk  and sk.le.nresamp
         sk = (tk+tlobHalf)/dtsamp
!*MS* why is the position of the particles in the first and last bin exchanged? 
!*MS* if sk.lt.1 => sk = 0
!*MS* if sk.gt.nresamp => sk = nresamp ? - should be like this!
         if(sk.lt.1)then
!*MS* sk= nresamp
           mm = sk/nresamp
           sk = sk + (1-mm)*nresamp
         endif
         if(sk.gt.nresamp)then
!*MS* sk= 0
           mm = sk/nresamp
           sk = sk - mm*nresamp
         endif

         nlo = int(sk) !*MS* cut comma digits [0,nresamp]
         fhi = sk - nlo !*MS* only comma digits [0,1)
         flo = 1-fhi
         nhi=nlo+1      !*MS* [1, nresamp+1]
         if(nlo.lt.1)nlo=nresamp !*MS* nlo = [1,nresamp], contains all particles in front of bunch in last bin
         if(nhi.gt.nresamp)nhi=1  !*MS*  nhi = [1, nresamp], contains all particles in the end of the bunch in frist bin

         if(plane.eq.'long')then
!           dpalmer is just line density         
            dpalmer(nlo) = dpalmer(nlo) + flo
            dpalmer(nhi) = dpalmer(nhi) + fhi
        endif
        if(plane.eq.'tran') then
!          same offset for both planes
	    xk = x(k)+coffset
            dptot(nlo) = dptot(nlo) + xk*flo
            dptot(nhi) = dptot(nhi) + xk*fhi
	endif



      enddo
      return
      end
! dpalmer, denlonoo, denlono, dptot


 
!**********************************************************************
!**********************************************************************

      subroutine btfwrite(tlob,dpwakex,dpwakey,dpwakes)
      include 'collider_time_evolution_common.f90'

      real :: tlob,dpwakex(nb),dpwakey(nb),dpwakes(nb)
      real :: cosbtfx,sinbtfx,cosbtfy,sinbtfy,cosbtfs,sinbtfs,dtsamp,dphibtf
      integer :: k

! calculate and write the sine and cosine components of the
! btf response
      cosbtfx=0
      sinbtfx=0
      cosbtfy=0
      sinbtfy=0
      cosbtfs=0
      sinbtfs=0
      dtsamp=tlob/nresamp
      dphibtf = harmbtf*dtsamp*2*pi/trev
!      write(6,*)'dphibtf= ',dphibtf
!      write(6,*)trev
!  transverse kick for positive ampbtf
        do k=1,nresamp
         cosbtfx = cosbtfx + cos(k*dphibtf)*dpwakex(k)
         sinbtfx = sinbtfx + sin(k*dphibtf)*dpwakex(k)
         cosbtfy = cosbtfy + cos(k*dphibtf)*dpwakey(k)
         sinbtfy = sinbtfy + sin(k*dphibtf)*dpwakey(k)
         cosbtfs = cosbtfs + cos(k*dphibtf)*dpwakes(k)
         sinbtfs = sinbtfs + sin(k*dphibtf)*dpwakes(k)
        enddo

      write(24,100)float(kturns),cosbtfx,sinbtfx,cosbtfy,sinbtfy,cosbtfs,sinbtfs
 100  format(8e14.4)
      return

      end


!**********************************************************************
!**********************************************************************


     subroutine resupdate(pnumber,np,np0,t,pt,px,py,emix,emiy,tlob,dpxtot,dpytot,dpalmer,denlono,denlonoo,&
                          coherecs,coherecss,wakeyy,wakeyx,wakes,dpwakex,dpwakey,dpwakes,iwrite)
      include 'collider_time_evolution_common.f90'
! does once per turn transverse kick for an impedance
! using an fft approach
!   tk = k*dtsamp, denk = denlon(k),  xdenk = dptot(k)
!   xdenk normalized so that xdenk = xk for small taupart
!   need to calculate the kick to px = betaL*dx/ds
!   <x>*I = dptot*cperbin, I=cperbin*denk
! integral (E+v\cross B)ds = integral (xI*wake)dt = volts
! delta(px) = ((trev*clight)/(2*pi*beta*etotdq*tunex))* integral[(xI)*wake*dt] 
! for array purposes put all particles in [0,tlob)
! stable fixed point of RF is tk=0


      real :: tlob, pnumber, dpxtot(nb),dpytot(nb),dpalmer(nb),denlono(nb),denlonoo(nb),emix,emiy,t(nb),&
              pt(nb),px(nb),py(nb),coherecs, coherecss
      real :: dtsamp, wakefieldx(nb),wakefieldy(nb),wakefields(nb),coeffsnrx,coeffsnry, ampx,ampy,&
              grv1,grv2,wrk,wyrk,wik,wyik,wakeyy(nb),wakeyx(nb),wakes(nb),dpwakex(nb),dpwakey(nb),dpwakes(nb),&
              etotdq, cperbin, dpx_coeff, dpy_coeff,dp_coeff,tk,sk,mm,fhi,flo,dpyk,dpxk,dpsk,omegabtf
      real :: wakesmax, cm, wamplitude,dpcoeff_voltage
      integer :: k, isgn,np0,np,nlo,nhi,iwrite

      dtsamp = tlob/nresamp
      etotdq = aatom*938.e6*gamma0/qatom
       
      cperbin = pnumber*1.602e-19*qatom/(np0*dtsamp)
      dpx_coeff =(cperbin*trev*clight*dtsamp)/(2*pi*beta*etotdq*tunex)
      dpy_coeff =(cperbin*trev*clight*dtsamp)/(2*pi*beta*etotdq*tuney)
      
      ampx=Sqrt(emix*betax)
      ampy=sqrt(emiy*betay)

! first we need to calculate the wakefields
        do k=1,nresamp
         wakefieldy(2*k-1)=dpytot(k)
         wakefieldy(2*k)=0
         wakefieldx(2*k-1)=dpxtot(k)
         wakefieldx(2*k)=0
! filter cooling, denlonoo is from early pickup
         if(ntdelay.eq.0)then
! one turn delay
                 wakefields(2*k-1)=dpalmer(k)-denlonoo(k)
         else
! cascaded one turn delays
                 wakefields(2*k-1)=dpalmer(k)-2*denlonoo(k)+denlono(k)
         endif
         wakefields(2*k)=0
! fill in memory array
!         denlono(k)=denlon(k)
        enddo
! add in the noise to the transverse signals
      coeffsnrx = 1.732*ampx*sqrt(float(np)/float(nresamp))*snrinv
      coeffsnry = 1.732*ampy*sqrt(float(np)/float(nresamp))*snrinv
      do k=1,nresamp
         call get2gaussrv(iseed,grv1,grv2)
         wakefieldy(2*k-1)= wakefieldy(2*k-1)+coeffsnry*grv1
         wakefieldx(2*k-1)= wakefieldx(2*k-1)+coeffsnrx*grv2
      enddo

        

      isgn=1
      call four1(wakefieldy,nresamp,isgn)
      call four1(wakefieldx,nresamp,isgn)
      call four1(wakefields,nresamp,isgn)
! multiply 
      do k=1,nresamp
! vertical
         wrk = wakefieldy(2*k-1)
         wik = wakefieldy(2*k)
         wyrk = wakeyy(2*k-1)
         wyik = wakeyy(2*k)
         wakefieldy(2*k-1)=wrk*wyrk-wik*wyik
         wakefieldy(2*k)=wrk*wyik+wik*wyrk
! horizontal
         wrk = wakefieldx(2*k-1)
         wik = wakefieldx(2*k)
         wyrk = wakeyx(2*k-1)
         wyik = wakeyx(2*k)
         wakefieldx(2*k-1)=wrk*wyrk-wik*wyik
         wakefieldx(2*k)=wrk*wyik+wik*wyrk



! longitudinal
         wrk = wakefields(2*k-1)
         wik = wakefields(2*k)
         wyrk = wakes(2*k-1)
         wyik = wakes(2*k)
         wakefields(2*k-1)=wrk*wyrk-wik*wyik
         wakefields(2*k)=wrk*wyik+wik*wyrk
      enddo

!!!***MS*** get real voltage on the kicker

! need longitudinal kick coefficient    (moved from below) 
      dp_coeff = gamma0*dtsamp*cperbin/etotdq
! calcualte the peak voltage per cavity
! scaled to the actual beam in volts
      dpcoeff_voltage = 2*dp_coeff*(etotdq/gamma0)*sqrt(np0/pnumber)
      wakesmax = vmaxreal/dpcoeff_voltage
! need to modify wakefields so vmaxreal is not exceeded
      do k=1,nresamp
         wamplitude = sqrt(wakefields(2*k)**2 + wakefields(2*k-1)**2)
         if(wamplitude.gt.wakesmax) then
            cm = wakesmax/wamplitude
            wakefields(2*k) = cm*wakefields(2*k)
            wakefields(2*k-1) = cm*wakefields(2*k-1)
         endif
      enddo

! write real voltages to file if writeStochCoolSwitch==1
  101  format(i6,1p20e8.1)   
        if (iwrite.eq.1) then
          if(writeStochCoolSwitch.eq.1) then
             open(unit=40,file='LongitudinalVoltages.out',access='append')
	     write(40,101) kturns, dpcoeff_voltage*sqrt(wakefields(2*5)**2 + wakefields(2*5-1)**2),&
                        dpcoeff_voltage*sqrt(wakefields(2*22)**2 + wakefields(2*22-1)**2) !for k=5 and k=22
           close(40)
         endif
       endif

!!!***MS***


! invert
      isgn=-1
      call four1(wakefieldy,nresamp,isgn)         
      call four1(wakefieldx,nresamp,isgn)         
      call four1(wakefields,nresamp,isgn)    

      do k=1,nresamp
         dpwakey(k)=dpy_coeff*wakefieldy(2*k-1)
         dpwakex(k)=dpx_coeff*wakefieldx(2*k-1)
! disable notch
!         dpwakenotch(k) = (1-enotch)*dpwakenotch(k) + enotch*dpwake(k)
!         dpwakenotch2(k) = (1-enotch)*dpwakenotch2(k) + enotch*dpwake2(k)
!         dpwake(k) = dpwake(k)-dpwakenotch(k)
!         dpwake2(k) = dpwake2(k)-dpwakenotch2(k)
         coherecs = coherecs + dpwakey(k)**2
         dpwakes(k)=dp_coeff*wakefields(2*k-1)
         coherecss = coherecss + dpwakes(k)**2
      enddo
! write the kicks 
      call btfwrite(tlob,dpwakex,dpwakey,dpwakes)
! make sure sums are ok. they are
      write(22,100)float(kturns),coherecss,coherecs,dp_coeff,wakes(2*30)
! return if before time
      if(kturns.lt.nturnon)return
      if(kturns.gt.nturns)return

!***************************************
!*MS* apply kick
      
      do k=1,np
         tk = t(k)
! 1.le.sk  and sk.le.nresamp
         sk = (tk+tlob/2)/dtsamp
         if(sk.lt.1)then
           mm = sk/nresamp
           sk = sk + (1-mm)*nresamp
         endif

         if(sk.gt.nresamp)then
           mm = sk/nresamp
           sk = sk - mm*nresamp
         endif
         nlo = int(sk)
         fhi = sk - nlo
         flo = 1-fhi
         nhi=nlo+1
         if(nlo.lt.1)nlo=nresamp
         if(nhi.gt.nresamp)nhi=1
         dpyk = dpwakey(nlo)*flo + dpwakey(nhi)*fhi
         dpxk = dpwakex(nlo)*flo + dpwakex(nhi)*fhi
         dpsk = dpwakes(nlo)*flo + dpwakes(nhi)*fhi
         py(k) = py(k) + dpyk
         px(k) = px(k) + dpxk
         pt(k) = pt(k) + dpsk
      enddo
 100  format(1p10e14.6)
      if(kturns.ne.nturnon)return
! loop through the btf kick
      omegabtf = 2*pi*harmbtf/trev
! transverse kick for positive amplitude
      if(ampbtf.gt.0)then
        do k=1,np
         tk = t(k)+tlob/2
         py(k)=py(k)+ampbtf*sin(tk*omegabtf+phibtf)
         px(k)=px(k)+ampbtf*sin(tk*omegabtf+phibtf)
        enddo
      else
! longitudinal
        do k=1,np
         tk = t(k)+tlob/2
         pt(k)=pt(k)+ampbtf*sin(tk*omegabtf+phibtf)
        enddo
      endif


      return
      end


!**********************************************************************
!**********************************************************************


      subroutine getimped(tlob,wakeyy,wakeyx,wakes)
      include 'collider_time_evolution_common.f90'
       
       real :: wsmoothedy(nb), wsmootheds(nb),wsmoothedx(nb),&
               wallsum,dtwake,df,fk,fk1,fra,gain,gains,fk0,argexp,&
               wsss,tk, wakeyy(nb),wakeyx(nb),wakes(nb),tlob
       integer :: isgn,k,k0
!      dimension zw0y(nb),zalfy(nb),zw0s(nb),zalfs(nb)
!      external finduct
!      dimension wsmoothed(nb), wsmootheds(nb),wsmoothedd(nb)
! everything in MKS
! all data are in pickup2.in
      dtwake = tlob/nresamp
      !wallsum=0
! get and plot the impedance
      isgn=1
!  wakey is vertical
      df =  1/(dtwake*nresamp)
      do k=1,nresamp
         if(k.le.nresamp/2+1)then
             fk = (k-1)*df
         else
             fk = (-nresamp+k-1)*df

         endif
         fk1 = abs(fk)
         if((fk1.ge.fimpedy0).and.(fk1.le.fimpedy1))then
! calculate
!*MS* if gainy1 = gainy0 => gain = const = gainy1 =gainy0!
           fra = (fk1-fimpedy0)/(fimpedy1-fimpedy0)
           gain = fra*gainy1 + (1-fra)*gainy0
         else
           gain=0
         endif
           wakeyy(2*k-1)=gain*nresamp*df
           wakeyy(2*k)=0
! wakeyx is horizontal
  
         fk1 = abs(fk)
         if((fk1.ge.fimpedx0).and.(fk1.le.fimpedx1))then
! calculate
           fra = (fk1-fimpedx0)/(fimpedx1-fimpedx0)
           gain = fra*gainx1 + (1-fra)*gainx0
         else
           gain=0
         endif
           wakeyx(2*k-1)=gain*nresamp*df
           wakeyx(2*k)=0



! longitudinal
        
         if((fk1.ge.fimpeds0).and.(fk1.le.fimpeds1))then
! calculate
           fra = (fk1-fimpeds0)/(fimpeds1-fimpeds0)
           gains = fra*gains1 + (1-fra)*gains0
         else
           gains=0
         endif
! Palmer is in phase, no phase shift
         wakes(2*k-1)=0
           if(fk.gt.0)then
              wakes(2*k)=gains*nresamp*df
           else
              wakes(2*k)=-gains*nresamp*df
           endif

       enddo

      open(unit=55,file='z.out',status='unknown')

      do k0=1,nresamp
! normalize wakey by nresamp, and smooth with equivalent length of taupart
         if(k0.le.nresamp/2+1)then
             fk0 = (k0-1)*df
         else
             fk0 = (-nresamp+k0-1)*df

         endif
         argexp = pi*(taupart*fk0)**2
!*MS* this .le.20 looks kind of random to me!
         if(argexp.le.20)then
            argexp = 1/exp(argexp)
         else
            argexp=0
         endif
         wakeyy(2*k0-1)=wakeyy(2*k0-1)/nresamp !*MS* =gainy/tlob
         wakeyy(2*k0)=wakeyy(2*k0)/nresamp     !*MS* =0
         wakeyx(2*k0-1)=wakeyx(2*k0-1)/nresamp !*MS* =gainx/tlob
         wakeyx(2*k0)=wakeyx(2*k0)/nresamp     !*MS* =0
         wakes(2*k0-1)=wakes(2*k0-1)/nresamp   !*MS* =0
         wakes(2*k0)=wakes(2*k0)/nresamp       !*MS* =gains/tlob
!*MS* the fk written in file 55 is the fk of k=nresamp (which is negative), since it does not change inside this loop.
!*MS* not much information in fk
!*MS* wakeyy(2*k0-1)/df = gain
!*MS* wakeyy(2*k0)/df = 0
!*MS* argexp completly unimportant: only calculated above but not used anywhere, and taupart is inputparameter, which is also not used elsewhere!
!*MS* -wakes(2*k0-1)/df = 0
!*MS* -wakes(2*k0)/df = -gains (for small k this is negative and for large k this is positive)
         write(55,100)fk,(k0-1)*df,wakeyy(2*k0-1)/df,wakeyy(2*k0)/df,argexp,-wakes(2*k0-1)/df,&
                     -wakes(2*k0)/df,wakeyx(2*k0-1)/df, wakeyx(2*k0)/df
  

!*MS* Everything which comes below is only for output

         wsmoothedy(2*k0)= wakeyy(2*k0)      !*MS* =0
         wsmoothedy(2*k0-1)= wakeyy(2*k0-1)  !*MS* = gain/tlob
         wsmootheds(2*k0)= wakes(2*k0)       !*MS* = gains/tlob
         wsmootheds(2*k0-1)= wakes(2*k0-1)   !*MS* = 0
         wsmoothedx(2*k0)= wakeyx(2*k0)      !*MS* =0
         wsmoothedx(2*k0-1)= wakeyx(2*k0-1)  !*MS* = gain/tlob

       enddo
      close(55)
 100  format(1p10e14.4)
! invert the fft for wsmoothed 
      isgn=-1
      call four1(wsmoothedy,nresamp,isgn)
      call four1(wsmootheds,nresamp,isgn)
      call four1(wsmoothedx,nresamp,isgn)
! write to file ws.out
      open(unit=55,file='ws.out',status='unknown')
      wsss=0
      do k=nresamp/2+2,nresamp
         if(k.le.nresamp/2+1)then
            tk = (k-1)*dtwake
         else
            tk = (k-nresamp-1)*dtwake
         endif
         write(55,100)tk,wsmoothedy(2*k-1),-wsmootheds(2*k-1),wsmoothedx(2*k-1)

         wsss = wsss + wsmootheds(2*k-1)
      enddo

!WHY do you do that twice?
      do k=1,nresamp/2+1
         if(k.le.nresamp/2+1)then
            tk = (k-1)*dtwake
         else
            tk = (k-nresamp-1)*dtwake
         endif
         write(55,100)tk,wsmoothedy(2*k-1),-wsmootheds(2*k-1),wsmoothedx(2*k-1),-wsss

         wsss = wsss + wsmootheds(2*k-1)

      enddo
      close(55)
      return
      end

!**********************************************************************
!**********************************************************************




