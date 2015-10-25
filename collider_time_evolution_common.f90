!     version 1.2

      implicit none
      integer, parameter :: nb=6000000
      integer, parameter :: nE=13500
      integer, parameter :: nIbs=120
      real, parameter :: pi=3.1415926536
      real, parameter ::  clight = 2.998e8

      integer :: nturns,nMacro1,nMacro2,nwrite,nresamp,nperturn,kturns,writeMountSwitch, &
        nturnon,nharm,nharm2,nBins,nextra,iseed,nwrite0,writeAllCoordSwitch,&
        ntdelay,nbunches,nIPs,longIntBins,longCoordMethod,coupleIBS, &
        RFswitch,betatronSwitch,raddampSwitch,ibsSwitch,collimationSwitch, &
        collisionSwitch,angleSwitch,collimAvgSwitch,nElem,xbins,ybins,zbins, nAperCuts, &
        stochCoolSwitch,writeStochCoolSwitch,&
        levellingSwitch,levIPSwitch(10)

      character*120 :: filebeam1,filebeam2
      character*2 :: collRoutine
      character*10 :: ibsMethod
      character*6 :: radMethod
      character*5 :: emitMethod

      real :: gammat,circ,gamma0,vrf,trf,omegarf,tcoeff,betax,betay, &
       pcoeff,frf,tunex,tuney,xinject,vrf2,beta,coulombLog,&
       vrev,trev,chromx,chromy,taupart,aatom,qatom,power,&
       gainx0,gainx1,fimpedx0,fimpedx1,fracmixbad,tlob1,tlob2,thib,&
       gainy0,gainy1,fimpedy0,fimpedy1,&
       gains0,gains1,fimpeds0,fimpeds1,phipks,cutoffAmpl,&
       fracibstot,fracibsx,ampbtf,phibtf,harmbtf,&
       alint,dqmin,phipk,sigI,betaS(10),theta(10),&
       tradlong,tradperp,siglong,sigperp,timeRatio,IPmult(10),&
       betx(nE),bety(nE),dispx(nE),alfx(nE),alfy(nE),dispxP(nE),&
       leng(nE),dispy(nE),dispyP(nE),bendAng(nE),K1L(nE),K1S(nE),xcut,ycut,&
       xmin,xmax,ymin,ymax,zmin,zmax,nSigCutBeta,nSigCutMom,dispxMom,&
       Ap(nIbs,nIbs,nIbs),Ax(nIbs,nIbs,nIbs),Ay(nIbs,nIbs,nIbs),&
       pnumInterp,eqTime,rho0,bunchLenPrecis,refEmxy,refMomSpread,betaxMom,&
       snrinv,coffsetx,coffsety,fracx,fracy,fracs, dispu, vmaxreal, &
       betaSMin(10),lumimax(10)

      real :: gTab(1502,2)

      common/inputi/nturns,nMacro1,nMacro2,nwrite,nresamp,nbunches,eqTime,rho0,bunchLenPrecis,dispxMom,&
        nSigCutMom,betaxMom,nharm,nharm2,nBins,nextra,betaS,theta,sigI,angleSwitch,cutoffAmpl,&
        nIPs,nAperCuts,longIntBins,xcut,ycut,betax,betay,filebeam1,filebeam2,longCoordMethod,&
        timeRatio,coupleIBS,RFswitch,betatronSwitch,raddampSwitch,collimAvgSwitch,coulombLog,&
        ibsSwitch,collimationSwitch,collisionSwitch,writeMountSwitch,IPmult,collRoutine,&
        stochCoolSwitch,&
        betaSMin,lumimax,levellingSwitch,levIPSwitch


      common/input2/gTab,nElem,betx,bety,dispx,alfx,&
           alfy,dispxP,leng,dispy,nSigCutBeta,writeAllCoordSwitch,&
           dispyP,bendAng,K1L,K1S,ibsMethod,radMethod,emitMethod

      common/input3/xmin,xmax,ymin,ymax,zmin,zmax,Ap,Ax,Ay,refEmxy,refMomSpread,&
        xbins,ybins,zbins,pnumInterp

      common/rfstuff/gammat,circ,gamma0,vrf,trf,omegarf,tcoeff,&
       pcoeff,frf,tunex,dqmin
      common/rad/tradlong,tradperp,siglong,sigperp
      common/inject/xinject,vrf2
      common/dynamics/beta,vrev,trev
      common/tran_in/tuney,chromx,chromy,taupart,alint
      common/randomstuff/iseed
      common/species/aatom,qatom,power
!      common/tunebeta/xtune(nb)
      common/resden/tlob1,tlob2,thib
!      common/diag/avgchrom,avgtune,csfull,dptrf(nb)
!      common/multi/zcbm,xavg0,pavg0
!      common/imped/wakey(nb),wakefield(nb),wakes(nb),wakefields(nb),
!     : waked(nb),wakefieldd(nb)
!      common/radstuff/tradperp,sigperp,tradlong,siglong
      common/tstuff/kturns,nturnon,nperturn,nwrite0
!      common/diag/coherecs,coherecss,
!     : avgline(nb)
      common/feedbackx/gainx0,gainx1,fimpedx0,fimpedx1,fracmixbad
      common/feedbacky/gainy0,gainy1,fimpedy0,fimpedy1
      common/feedbacks/gains0,gains1,fimpeds0,fimpeds1,phipks,vmaxreal
      common/ibstuff/fracibstot,fracibsx,ampbtf,phibtf,harmbtf
!      common/lumistuff/peakilumi,fwhmlumi
      common/coolingstuff/ntdelay,snrinv,writeStochCoolSwitch,coffsetx,coffsety,fracx,fracy,fracs,dispu
   
