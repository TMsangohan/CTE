!     version 1.1

      implicit none
      integer, parameter :: nb=6000000
      integer, parameter :: nE=13500
      integer, parameter :: nIbs=120
      real, parameter :: pi=3.1415926536
      real, parameter ::  clight = 2.998e8

      integer :: nturns,nMacro,nwrite,nresamp,nperturn,kturns,writeMountSwitch, &
        nturnon,nharm,nharm2,nBins,nextra,iseed,nwrite0,writeAllCoordSwitch,&
        blowupSwitch,ntdelay,nbunches,nIPs,longIntBins,longCoordMethod,transvCoordMethod,coupleIBS,&
        RFswitch,betatronSwitch,raddampSwitch,ibsSwitch,collimationSwitch, &
        collisionSwitch,angleSwitch,collimAvgSwitch,nElem,xbins,ybins,zbins, nAperCuts

      character*200 :: filebeam1,filebeam2
      character*2 :: collRoutine
      character*10 :: ibsMethod
      character*6 :: radMethod
      character*5 :: emitMethod,blowupMethod

      real :: gammat,circ,gamma0,vrf,trf,omegarf,tcoeff,betax,betay, &
       pcoeff,frf,tunex,tuney,xinject,vrf2,beta,coulombLog,k2L,k2Lskew,&
       vrev,trev,chromx,chromy,taupart,aatom,qatom,power,&
       gain0,gain1,fimped0,fimped1,fracmixbad,tlob,thib,&
       gains0,gains1,fimpeds0,fimpeds1,phipks,cutoffAmpl,&
       fracibstot,fracibsx,ampbtf,phibtf,harmbtf,&
       alint,dqmin,phipk,sigI,betaS(10),theta(10),pxKickFac,pyKickFac,&
       tradlong,tradperp,siglong,sigperp,timeRatio,IPmult(10),&
       betx(nE),bety(nE),dispx(nE),alfx(nE),alfy(nE),dispxP(nE),&
       leng(nE),dispy(nE),dispyP(nE),bendAng(nE),K1L(nE),xcut,ycut,&
       xmin,xmax,ymin,ymax,zmin,zmax,nSigCutBeta,nSigCutMom,dispxMom,&
       Ap(nIbs,nIbs,nIbs),Ax(nIbs,nIbs,nIbs),Ay(nIbs,nIbs,nIbs),&
       pnumInterp,eqTime,rho0,bunchLenPrecis,refEmxy,refMomSpread,betaxMom

      real :: gTab(1502,2)

      common/inputi/nturns,nMacro,nwrite,nresamp,nbunches,eqTime,rho0,bunchLenPrecis,dispxMom,&
        nSigCutMom,betaxMom,nharm,nharm2,nBins,nextra,betaS,theta,sigI,angleSwitch,cutoffAmpl,&
        nIPs,nAperCuts,longIntBins,xcut,ycut,betax,betay,filebeam1,filebeam2,longCoordMethod,transvCoordMethod,&
        timeRatio,coupleIBS,RFswitch,betatronSwitch,raddampSwitch,collimAvgSwitch,coulombLog,&
        pxKickFac,pyKickFac,&
        ibsSwitch,blowupSwitch,collimationSwitch,collisionSwitch,writeMountSwitch,IPmult,collRoutine

      common/input2/gTab,nElem,betx,bety,dispx,alfx,&
           alfy,dispxP,leng,dispy,nSigCutBeta,writeAllCoordSwitch,&
           dispyP,bendAng,K1L,ibsMethod,radMethod,emitMethod,blowupMethod

      common/input3/xmin,xmax,ymin,ymax,zmin,zmax,Ap,Ax,Ay,refEmxy,refMomSpread,&
        xbins,ybins,zbins,pnumInterp

      common/rfstuff/gammat,circ,gamma0,vrf,trf,omegarf,tcoeff,&
       pcoeff,frf,tunex,dqmin,k2L,k2Lskew
      common/rad/tradlong,tradperp,siglong,sigperp
      common/inject/xinject,vrf2
      common/dynamics/beta,vrev,trev
      common/tran_in/tuney,chromx,chromy,taupart,alint
      common/randomstuff/iseed
      common/species/aatom,qatom,power
!      common/tunebeta/xtune(nb)
      common/resden/tlob,thib
!      common/diag/avgchrom,avgtune,csfull,dptrf(nb)
!      common/multi/zcbm,xavg0,pavg0
!      common/imped/wakey(nb),wakefield(nb),wakes(nb),wakefields(nb),
!     : waked(nb),wakefieldd(nb)
!      common/radstuff/tradperp,sigperp,tradlong,siglong
      common/tstuff/kturns,nturnon,nperturn,nwrite0
!      common/diag/coherecs,coherecss,
!     : avgline(nb)
      common/feedback/gain0,gain1,fimped0,fimped1,fracmixbad
      common/feedbacks/gains0,gains1,fimpeds0,fimpeds1,phipks
      common/ibstuff/fracibstot,fracibsx,ampbtf,phibtf,harmbtf
!      common/lumistuff/peakilumi,fwhmlumi


