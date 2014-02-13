C  UNIX VERSION OF NBS CODES
*QFS
C  CALCULATE QUASIELASTIC, TWO-NUCLEON, RESONANCE, AND X-SCALING
C  ELECTRON SCATTERING CROSS SECTION FOR ARBITRARY NUCLEI
C
C  WRITTEN BY J.W. LIGHTBODY JR. AND J.S. O'CONNELL
C  NATIONAL BUREAU OF STANDARDS, GAITHERSBURG, MD 20899
C  OCTOBER 1987
! all the variables finishing by al were added to the original code
! and all the comment with ! 
! and all the non-capital text. 
! CEBAF A. Deur 07/31/98
! the cebaf spectrometer acceptance is taken in account: horizontal:+/-28 mr
                                                       ! vertical: +/-60 mr
C Modified by C. Gu for a library used by Geant4 similation

      subroutine qfs(tgt_Z,tgt_A,Ei,Ep,ang,xs,EPS1,EPSD1,PF1,Tb1,Ta1)

      IMPLICIT double precision (A-H,O-Z)

      integer tgt_Z,tgt_A 
      integer err,erre
      COMMON /PAR/E,TH,W,Z,A,EPS,EPSD,PF,SPENCE
      COMMON /ADDONS/SCALE,Tb,Ta
      
      SCALE=1.D-26 ! cm^2/MeV-sr
      PM=939.D0
      DM=1232.D0
      ALPH=1./137.03604D0
      HBARC=197.32858D0
      PI=ACOS(-1.)
      SIGAVE=0.0
     
      W=Ei-Ep
      E=Ei
      Z=dble(tgt_Z)
      A=dble(tgt_A)
      EPS=EPS1
      EPSD=EPSD1
      PF=PF1
      TH=ang
      Tb=Tb1
      Ta=Ta1

      SIGQFZA=SIGQFS(E,TH,W,Z,A,EPS,PF)*SCALE
      SIGDA=SIGDEL(E,TH,W,A,EPSD,PF)*SCALE
      SIGXA=SIGX(E,TH,W,A)*SCALE
      SIGR1A=SIGR1(E,TH,W,A,PF)*SCALE
      SIGR2A=SIGR2(E,TH,W,A,PF)*SCALE
      SIG2NA=SIG2N(E,TH,W,Z,A,PF)*SCALE
      SIG=SIGQFZA+SIGDA+SIGXA+SIGR1A+SIGR2A+SIG2NA
      xs=SIG*1.D+30 ! ub/MeV-sr
      if ((TB.gt.0.0D0).and.(Ta.gt.0.0D0)) then
      CALL RADIATE(E,TH,W,SIG,SIGRAD)
      xs=SIGRAD*1.D+30
      endif
      err=erre
      END
*FD
      double precision FUNCTION FD(QMS,A)
      IMPLICIT double precision (A-H,O-Z)
      FD=1./(1.+QMS/A**2)**2
      RETURN
      END
*FM
      double precision FUNCTION FM(QMS,A)
      IMPLICIT double precision (A-H,O-Z)
      FM=1./(1.+QMS/A**2)
      RETURN
      END
*FPHENOM
      double precision FUNCTION FPHENOM(QMS)
      IMPLICIT double precision (A-H,O-Z)
      A1=.55
      A2=20./1.E6
      B1=.45
      B2=.45/1.E6
      C1=0.03
      C2=0.2/1.E12
      FPHENOM=A1*EXP(-A2*QMS)+B1*EXP(-B2*QMS)
      FPHENOM=FPHENOM+C1*EXP(-C2*(QMS-4.5E6)**2)
      FPHENOM=SQRT(FPHENOM)
      RETURN
      END
*FYUKAWA
      double precision FUNCTION FYUKAWA(QMS,A)
      IMPLICIT double precision (A-H,O-Z)
      IF(QMS.LT.1.E-5.OR.A.LT.1.E-5)THEN
      FYUKAWA=0.
      ELSE
      ARG=SQRT(QMS/2.)/A
      FYUKAWA=ATAN(ARG)/ARG
      ENDIF
      RETURN
      END
*SIGMOT
      double precision FUNCTION SIGMOT(E,THR)
      IMPLICIT double precision (A-H,O-Z)
      ALPH=1./137.03604
      HBARC=197.3286
      SIGMOT=(ALPH*HBARC*COS(THR/2.)/2./E/SIN(THR/2.)**2)**2
C  FM**2/SR
      RETURN
      END
*RECOIL
      double precision FUNCTION RECOIL(E,THR,TM)
      IMPLICIT double precision (A-H,O-Z)
      RECOIL=1./(1.+2.*E*SIN(THR/2.)**2/TM)
      RETURN
      END
*SIGX
      double precision FUNCTION SIGX(E,TH,W,A)
      IMPLICIT double precision (A-H,O-Z)
      ALPH=1./137.03604
      PI=ACOS(-1.)
C     SIG0=111.*1.E-30
      SIG0=100.D-4
C     SIG1=60.*1.E-27
      SIG1=54.*1.D-1
      PIMASS=140.
      PM=939.
C     GAM0=550.
      GAM0=650.
C     R=0.10
      AQ=250.
      THR=TH*PI/180.
      IF(W.LT.1.E-5)GO TO 4
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      ARG0=W-QMS/2./PM-PIMASS-PIMASS**2/2./PM
      ARG1=ARG0/GAM0
      ARG=ARG1**2/2.
      IF(ARG1.GT.8.)THEN
      SHAPE=1.+SIG1/SIG0/ARG0
      ELSEIF(ARG1.LT.1.E-5)THEN
      SHAPE=0.
      ELSEIF(ARG1.LT.0.1)THEN
      SHAPE=SIG1*ARG0/2./GAM0**2/SIG0
      ELSE
      SHAPE=(1.-EXP(-ARG))*(1.+SIG1/SIG0/ARG0)
      ENDIF
      EKAPPA=W-QMS/2./PM
      SIGGAM=SIG0*SHAPE
      QS=QMS+W**2
      EPS=1./(1.+2.*QS*TAN(THR/2.)**2/QMS)
      FLUX=ALPH*EKAPPA*(E-W)/2./PI**2/QMS/E/(1.-EPS)
      IF(FLUX.LT.1.E-20)FLUX=0.
      SIGEE=FLUX*SIGGAM*FPHENOM(QMS)**2
C     SIGEE=FLUX*SIGGAM
      R=0.56*1.E6/(QMS+PM**2)
      FACTOR1=1.+EPS*R
      SIGEE=SIGEE*FACTOR1
 4    SIGX=A*SIGEE
      RETURN
      END
*SIGR1
      double precision FUNCTION SIGR1(E,TH,W,A,PF)
      IMPLICIT double precision (A-H,O-Z)
      PI=ACOS(-1.)
      PM=939.
      PIMASS=140.
      THR=TH*PI/180.
      PFR=230.
      RM=1500.
      EPSR=0.
      AR0=1000.
      AR1=1000.
      GAMQFR=120.
      GAMSPRD=140.
      GAMR=110.
      GAMPI=5.
      QFRP=1.20D-7
      QMSQFR=4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSQFR=QMSQFR+115.**2
      QMSRR=4.*10000.*(10000.-1240.)*SIN(6.*PI/180./2.)**2
      QVSRR=QMSRR+1240.**2
      SIGREF=FD(QMSRR,AR0)**2*QVSRR
      SIGREF=SIGREF*(QMSRR/2./QVSRR+TAN(6.*PI/180./2.)**2)
      SIGREF=SIGREF*SIGMOT(10000.D0,6.*PI/180.)
      NA=INT(A)
      IF(NA.EQ.1)THEN
      QFR=QFRP
      GSPRDA=0.
      AR=AR0
      ELSEIF(NA.LT.4)THEN
      QFR=QFRP
      GSPRDA=(A-1.)*GAMSPRD/3.
      AR=AR0+(A-1.)*(AR1-AR0)/3.
      ELSE
      AR=AR1
      GSPRDA=GAMSPRD
      QFR=QFRP
      ENDIF
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      IF(NA.GT.1)THEN
      GAMQ=GAMQFR*PF*SQRT(QVS)/PFR/SQRT(QVSQFR)
      ELSE
      GAMQ=0.
      ENDIF
      CMTOT2=PM**2+2.*PM*W-QMS
      WTHRESH=4.*E**2*SIN(THR/2.)**2+PIMASS**2+2.*PIMASS*PM
      WTHRESH=WTHRESH/2./PM
      THRESHD=1.+PF/PM+PF**2/2./PM**2+2.*E*SIN(THR/2.)**2/PM
      WTHRESH=WTHRESH/THRESHD
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
      THRESH=0.
      ENDIF
      EPR=E-(RM-PM)*(RM+PM)/2./PM
      EPR=EPR/(1.+2.*E*SIN(THR/2.)**2/PM)
      EPR=EPR-EPSR
      WR=E-EPR
      GAM=SQRT(GAMR**2+GAMQ**2+GSPRDA**2)
      SIGR=QFR*(GAMR/GAM)/SIGREF
      SIGR=SIGR*CMTOT2*GAM**2
      SIGR=SIGR/((CMTOT2-(RM+EPSR)**2)**2+CMTOT2*GAM**2)
      SIGR=SIGR*QVS*FD(QMS,AR)**2
      SIGR=SIGR*(QMS/2./QVS+TAN(THR/2.)**2)
      SIGR=SIGR*SIGMOT(E,THR)
      SIGR1=A*THRESH*SIGR
      RETURN
      END
*SIGR2
      double precision FUNCTION SIGR2(E,TH,W,A,PF)
      IMPLICIT double precision (A-H,O-Z)
      PI=ACOS(-1.)
      PM=939.
      PIMASS=140.
      THR=TH*PI/180.
      PFR=230.
      RM=1700.
      EPSR=0.
      AR0=1200.
      AR1=1200.
      GAMQFR=120.
      GAMSPRD=140.
      GAMR=110.
      GAMPI=5.
      QFRP=0.68D-7
      QMSQFR=4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSQFR=QMSQFR+115.**2
      QMSRR=4.*10000.*(10000.-1520.)*SIN(6.*PI/180./2.)**2
      QVSRR=QMSRR+1520.**2
      SIGREF=FD(QMSRR,AR0)**2*QVSRR
      SIGREF=SIGREF*(QMSRR/2./QVSRR+TAN(6.*PI/180./2.)**2)
      SIGREF=SIGREF*SIGMOT(10000.D0,6.*PI/180.)
      NA=INT(A)
      IF(NA.EQ.1)THEN
      QFR=QFRP
      GSPRDA=0.
      AR=AR0
      ELSEIF(NA.LT.4)THEN
      QFR=QFRP
      GSPRDA=(A-1.)*GAMSPRD/3.
      AR=AR0+(A-1.)*(AR1-AR0)/3.
      ELSE
      AR=AR1
      GSPRDA=GAMSPRD
      QFR=QFRP
      ENDIF
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      IF(NA.GT.1)THEN
      GAMQ=GAMQFR*PF*SQRT(QVS)/PFR/SQRT(QVSQFR)
      ELSE
      GAMQ=0.
      ENDIF
      CMTOT2=PM**2+2.*PM*W-QMS
      WTHRESH=4.*E**2*SIN(THR/2.)**2+PIMASS**2+2.*PIMASS*PM
      WTHRESH=WTHRESH/2./PM
      THRESHD=1.+PF/PM+PF**2/2./PM**2+2.*E*SIN(THR/2.)**2/PM
      WTHRESH=WTHRESH/THRESHD
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
      THRESH=0.
      ENDIF
      EPR=E-(RM-PM)*(RM+PM)/2./PM
      EPR=EPR/(1.+2.*E*SIN(THR/2.)**2/PM)
      EPR=EPR-EPSR
      WR=E-EPR
      GAM=SQRT(GAMR**2+GAMQ**2+GSPRDA**2)
      SIGR=QFR*(GAMR/GAM)/SIGREF
      SIGR=SIGR*CMTOT2*GAM**2
      SIGR=SIGR/((CMTOT2-(RM+EPSR)**2)**2+CMTOT2*GAM**2)
      SIGR=SIGR*QVS*FD(QMS,AR)**2
      SIGR=SIGR*(QMS/2./QVS+TAN(THR/2.)**2)
      SIGR=SIGR*SIGMOT(E,THR)
      SIGR2=A*THRESH*SIGR
      RETURN
      END
*SIG2N
      double precision FUNCTION SIG2N(E,TH,W,Z,A,PF)
      IMPLICIT double precision (A-H,O-Z)
      PI=ACOS(-1.)
      THR=TH*PI/180.
      DM=1232.
      PIMASS=140.
      PM=940.
      A2=550.
      PFR=60.
      GAM2N=20.
      GAMQFR=40.
      GAMREF=300.
      GAMR=GAMREF
      SIGREF=0.20D-7
      QMSR=4.*596.8*(596.8-380.)*SIN(60.*PI/180./2.)**2
      QVSR=QMSR+380.**2
      SIGKIN=0.5*SIGMOT(596.8D0,60.*PI/180.)
      SIGKIN=SIGKIN*(QMSR/2./QVSR+TAN(60.*PI/180./2.)**2)
      SIGKIN=SIGKIN*QVSR*FD(QMSR,A2)**2
      SIGKIN=SIGKIN*GAMR/GAMREF
      SIGCON=SIGREF/SIGKIN
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      GAMQF=GAMQFR*(PF/PFR)*(SQRT(QVS)/SQRT(QVSR))
      EFFMASS=(PM+DM)/2.
      SIG=(Z*(A-Z)/A)*SIGMOT(E,THR)
      SIG=SIG*(QMS/2./QVS+TAN(THR/2.)**2)
      SIG=SIG*QVS*FD(QMS,A2)**2
      EKAPPA=W-QMS/2./PM
      CMTOT2=PM**2+2.*PM*EKAPPA
C     GAM=SQRT(GAMR**2+GAMQF**2)
      GAM=GAMR
      SIG=SIG*CMTOT2*GAM**2
      SIG=SIG/((CMTOT2-EFFMASS**2)**2+CMTOT2*GAM**2)
      SIG=SIG*(GAMR/GAM)*SIGCON
      SIG2N=SIG
      WTHRESH=QMS/4./PM
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAM2N)
      ELSE
      THRESH=0.
      ENDIF
      SIG2N=SIG2N*THRESH
      RETURN
      END
*SIGDEL
      double precision FUNCTION SIGDEL(E,TH,W,A,EPSD,PF)
      IMPLICIT double precision (A-H,O-Z)
      PM=939.
      PIMASS=140.
      DM=1219.
      AD1=700.
      AD0=774.
      PI=ACOS(-1.)
      ALPH=1./137.03604
      HBARC=197.32858
      GAMDP=110.
      GAMSPRD=140.
      GAMR=120.
      GAMPI=5.
      QFDP=1.02D-7
      PFR=230.
      QMSR=4.*730.*(730.-390.)*SIN(37.1*PI/180./2.)**2
      QVSR=QMSR+390.**2
      QMSRQ=4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSRQ=QMSRQ+115.**2
      NA=INT(A)
      IF(NA.EQ.1)THEN
      QFD=QFDP
      GSPRDA=0.
      AD=AD0
      ELSEIF(NA.LT.4)THEN
      QFD=QFDP
      GSPRDA=(A-1.)*GAMSPRD/3.
      AD=AD0+(A-1.)*(AD1-AD0)/3.
      ELSE
      AD=AD1
      GSPRDA=GAMSPRD
      QFD=QFDP
      ENDIF
      THR=TH*PI/180.
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      EKAPPA=W-QMS/2./PM
      CMTOT2=PM**2+2.*PM*EKAPPA
C  BEGIN DELTA CALCULATION
      IF(NA.GT.1)THEN
      GAMQ=GAMR*PF*SQRT(QVS)/PFR/SQRT(QVSRQ)
      ELSE
      GAMQ=0.
      ENDIF
      EPD=E-(DM-PM)*(DM+PM)/2./PM
      EPD=EPD/(1.+2.*E*SIN(THR/2.)**2/PM)
      EPD=EPD-EPSD
      WD=E-EPD
      QMSPK=4.*E*EPD*SIN(THR/2.)**2
      QVSPK=QMSPK+WD**2
C
C NOTE WIDTH INCLUDES E-DEPENDENCE,FERMI BROADENING,& SPREADING
C
      WTHRESH=4.*E**2*SIN(THR/2.)**2+PIMASS**2+2.*PIMASS*PM
      WTHRESH=WTHRESH/2./PM
      THRESHD=1.+PF/PM+PF**2/2./PM**2+2.*E*SIN(THR/2.)**2/PM
      WTHRESH=WTHRESH/THRESHD
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAMPI)
      ELSE
      THRESH=0.
      ENDIF
      GAMD=GAMDP
      GAM=SQRT(GAMD**2+GAMQ**2+GSPRDA**2)
      SIGD=QFDP*(GAMDP/GAM)
      SIGD=SIGD*CMTOT2*GAM**2
      SIGD=SIGD/((CMTOT2-(DM+EPSD)**2)**2+CMTOT2*GAM**2)
      SIGD=SIGD*FD(QMS,AD)**2/FD(QMSR,AD)**2
      TEST=QVS/QVSR
      SIGD=SIGD*TEST
      SIGD=SIGD*(QMS/2./QVS+TAN(THR/2.)**2)
      SIGD=SIGD/(QMSR/2./QVSR+TAN(37.1*PI/180./2.)**2)
      SIGD=SIGD*SIGMOT(E,THR)/SIGMOT(730.D0,37.1*PI/180.)
      SIGD=SIGD*A
      SIGD=SIGD*THRESH
      SIGDEL=SIGD
      RETURN
      END
*SIGQFS
      double precision FUNCTION SIGQFS(E,TH,W,Z,A,EPS,PF)
      IMPLICIT double precision (A-H,O-Z)
      PM=939.
      UP=2.7928456
      UN=-1.91304184
      AP0=840.
      AP1=750.
      ALPH=1./137.03604
      HBARC=197.32858
      PI=ACOS(-1.)
      GAMR=120.
      PFR=230.
      QMSRQ=4.*730.*(730.-115.)*SIN(37.1*PI/180./2.)**2
      QVSRQ=QMSRQ+115.**2
      NA=INT(A)
      IF(NA.EQ.1)THEN
      AP=AP0
      ELSEIF(NA.LT.4)THEN
      AP=AP0+(A-1.)*(AP1-AP0)/3.
      ELSE
      AP=AP1
      ENDIF
C     PRINT 200
C  200 FORMAT(' ENTER DE-E[MEV],DOMEGA-E[SR],B-LUMINOSITY[CM-2*S-1]')
C     READ *,DEE,DWE,BLUM
      THR=TH*PI/180.
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      EKAPPA=W-QMS/2./PM
      IF(EKAPPA.GT.-PM/2.)THEN
      CMTOT=SQRT(PM**2+2.*PM*EKAPPA)
      ELSE
      CMTOT=PM
      ENDIF
C  START QFS SECTION
      SIGNS=SIGMOT(E,THR)*RECOIL(E,THR,PM)
      FORMP=1.+QMS*UP**2/4./PM**2
      FORMP=FORMP/(1.+QMS/4./PM**2)
      FORMP=FORMP+QMS*UP**2*TAN(THR/2.)**2/2./PM**2
      FORMP=FORMP*FD(QMS,AP)**2
      SIGEP=SIGNS*FORMP
      FALLOFF=(1.+5.6*QMS/4./PM**2)**2
      FORMN=(QMS*UN/4./PM**2)**2/FALLOFF+QMS*UN**2/4./PM**2
      FORMN=FORMN/(1.+QMS/4./PM**2)
      FORMN=FORMN+QMS*UN**2*TAN(THR/2.)**2/2./PM**2
      FORMN=FORMN*FD(QMS,AP)**2
      SIGEN=SIGNS*FORMN
      EPQ=4.*E**2*SIN(THR/2.)**2/2./PM
      EPQ=EPQ/(1.+2.*E*SIN(THR/2.)**2/PM)+EPS
      EPQ=E-EPQ
      IF(INT(A).EQ.1)THEN
      ARG=(E-W-EPQ)/SQRT(2.)/1.
      DEN=2.51
      ELSE
      GAMQ=GAMR*PF*SQRT(QVS)/PFR/SQRT(QVSRQ)
      ARG=(E-W-EPQ)/1.20/(GAMQ/2.)
      DEN=2.13*(GAMQ/2.)
      ENDIF
      NQ=INT(ARG)
      IF(ABS(NQ).GT.10)THEN
      SIGQ=0.
      ELSE
      SIGQ=(Z*SIGEP+(A-Z)*SIGEN)*EXP(-ARG**2)/DEN
      ENDIF
      SIGQFS=SIGQ
      RETURN
      END
* RADIATE
      SUBROUTINE RADIATE(E,TH,W,SIGNR,SIGRAD)
C     DOES NOT INCLUDE CONTRIBUTION FROM ELASTIC SCATTERING
C
C-----K. Slifer. 09/16/02
C
C     Rewrote subroutine to include external bremsstrahlung using
C     formalism from S. Stein et al. Phys. Rev. D 12 7. Equation (A82)
C     Where possible the equation number is given with each expression.
C----------------------------------------------------------------------
      IMPLICIT double precision (A-H,O-Z)
      COMMON/PAR/E0,TH0,W0,Z,A,EPS,EPSD,PF,SPENCE
      COMMON/ADDONS/SCALE,Tb,Ta
      DEL   = 10
      PREC  = .001  ! 0.001 gets rid of glitches

      xb    = 4./3.  ! A45
      XM    = 931.49 ! Mass of the nucleon
      XMT   = A*XM   ! Mass of the target
      ALPH  = 1./137.03604
      EMASS = 0.511
      PI    = ACOS(-1.)
      THR   = TH*PI/180.
      ARG   = COS(THR/2.)**2

      SPENCE= PI**2/6.-LOG(ARG)*LOG(1.-ARG)
      DO 10 NSP=1,50
 10     SPENCE = SPENCE-ARG**NSP/FLOAT(NSP)**2

      QMS   = 4.*E*(E-W)*SIN(THR/2.)**2

      D1=(2.*ALPH/PI)*(LOG(QMS/EMASS**2)-1.)      ! =2b*tr (A57)
      tr=D1/2./xb

      D2 = 13.*(LOG(QMS/EMASS**2)-1.)/12.-17./36. ! this term dominates D2
      D2 = D2 +0.5*(PI**2/6.-SPENCE)
      D2 = D2 -1./4.*( LOG( E/(E-W) ) )**2        ! Correct. to peak. appr.
      D2 = D2*(2.*ALPH/PI)                 
      D2 = D2+0.5772*xb*(Tb+Ta)                   ! Here D2= F-1
      xF = (1.+D2)                                ! (A44)

      Tpb = tr + Tb
      Tpa = tr + Ta  
   
      R   = ( XMT+E*(1-COS(THR)) )/( XMT-(E-W)*(1-COS(THR)) ) ! (A83)
      eta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )    ! (A46)
      xi  = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi  = xi/( (Z+eta)*LOG(183.*Z**(-1./3.)) )              ! (A52)

      SIGRAD = SIGNR * xF
      SIGRAD = SIGRAD*( (R*DEL/E  )**(xb*Tpb) )
      SIGRAD = SIGRAD*( (DEL/(E-W))**(xb*Tpa) )
      SIGRAD = SIGRAD*(1. - xi/DEL/( 1.-xb*(Tpb+Tpa)) )

      TERM1=(R*DEL/E  )**(xb*Tpb)
      TERM2=(DEL/(E-W))**(xb*Tpa)
      TERM3=(1. - xi/DEL/( 1.-xb*(Tpb+Tpa)) )
      TERM4=xF
C
C-----Stein's 1st integral wrt dEs' (A82)
C
C     limits of 0 to W-DEL give almost same results
C
      X1LO   = (E-W)*( XMT/( XMT-2.*(E-W)*(SIN(THR/2.))**2) -1.0 )
      X1HI   = W-R*DEL
      ANS_Es = 0.
      IF (X1HI.GT.X1LO) THEN
        CALL ROM(X1LO,X1HI,PREC,ANS_Es,KF,1)
      ENDIF
      ANS_Es = ANS_Es * SCALE
C
C-----Stein's 2nd integral wrt dEp' (A82)
C
C     limits of 0 to W-DEL give almost same results
C
      X2LO   = E*( 1.0-XMT/( XMT+2.*E*(SIN(THR/2.))**2) )
      X2HI   = W-DEL 
      ANS_Ep = 0.
      IF (X2HI.GT.X2LO) THEN
         CALL ROM(X2LO,X2HI,PREC,ANS_Ep,KF,2)
      ENDIF
      ANS_Ep = ANS_Ep * SCALE
      
      SIGRAD = SIGRAD + ANS_Es + ANS_Ep

      RETURN
      END

*VALY
      SUBROUTINE VALY(X,F,IFUNC)
      IMPLICIT double precision (A-H,O-Z)
C     REAL*8 FUNCTION F(X)
      COMMON/PAR/E,TH,W,Z,A,EPS,EPSD,PF,SPENCE
      COMMON/ADDONS/SCALE,Tb,Ta

      ALPH  = 1./137.03604
      EMASS = 0.511
      PI    = ACOS(-1.)
      THR   = TH*PI/180.
      xb    = 4./3.
      XM    = 931.49 ! Mass of the nucleon
      XMT   = A*XM   ! Mass of the target

      eta = LOG(1440.*Z**(-2./3.) )/LOG(183.*Z**(-1./3.) )
      xi  = (PI*EMASS/2./ALPH)*(Ta+Tb)
      xi  = xi/( (Z+eta)*LOG(183.*Z**(-1./3.)) )
      R   = ( XMT+E*(1-COS(THR)) )/( XMT-(E-W)*(1-COS(THR)) )
C
C-----Stein's 2nd integral dEp'
C
      QMS2  = 4.*E*(E-X)*SIN(THR/2.)**2  ! 1/15/03
      tr2   = 1./xb*(ALPH/PI)*(LOG(QMS2/EMASS**2)-1.)
      Tpb   = tr2 + Tb
      Tpa   = tr2 + Ta

      D2    = 13.*(LOG(QMS2/EMASS**2)-1.)/12.-17./36.
      D2    = D2 - 1./4.*( LOG( E/(E-W) ) )**2 !KS. Correction to peak. approx.
      D2    = D2 + 0.5*(PI**2/6.-SPENCE)
      D2    = D2 * (2.*ALPH/PI)
      D2    = D2 + 0.5772*xb*(Tb+Ta)

      SIG2  = SIGQFS(E,TH,X,Z,A,EPS,PF)
      SIG2  = SIG2 + SIGDEL(E,TH,X,A,EPSD,PF)
      SIG2  = SIG2 + SIGX(E,TH,X,A)
      SIG2  = SIG2 + SIGR1(E,TH,X,A,PF)
      SIG2  = SIG2 + SIGR2(E,TH,X,A,PF)
      SIG2  = SIG2 + SIG2N(E,TH,X,Z,A,PF)
      F2    = ( xb*Tpa/(W-X) ) *phi((W-X)/(E-X))
      F2    = F2 + xi/(2.*(W-X)**2)
      F2    = F2 * SIG2*(1.+D2)
      F2    = F2 * ( (W-X)/(E-X) )**(xb*Tpa)
      F2    = F2 * ( (W-X)*R/(E) )**(xb*Tpb)
C
C-----Stein's 1st integral dEs'
C
      QMS1  = 4.*(E-W+X)*(E-W)*SIN(THR/2.)**2    !    1/15/03
      tr1   = 1./xb*(ALPH/PI)*(LOG(QMS1/EMASS**2)-1.)
      Tpb   = tr1 + Tb
      Tpa   = tr1 + Ta

      D2    = 13.*(LOG(QMS1/EMASS**2)-1.)/12.-17./36.
      D2    = D2 - 1./4.*( LOG( E/(E-W) ) )**2 !Corr. to peak. approx.
      D2    = D2 + 0.5*(PI**2/6.-SPENCE)
      D2    = D2 * (2.*ALPH/PI)
      D2    = D2 + 0.5772*xb*(Tb+Ta) ! 1/14/02

      SIG1  = SIGQFS(E-W+X,TH,X,Z,A,EPS,PF)
      SIG1  = SIG1 + SIGDEL(E-W+X,TH,X,A,EPSD,PF)
      SIG1  = SIG1 +   SIGX(E-W+X,TH,X,A)
      SIG1  = SIG1 +  SIGR1(E-W+X,TH,X,A,PF)
      SIG1  = SIG1 +  SIGR2(E-W+X,TH,X,A,PF)
      SIG1  = SIG1 +  SIG2N(E-W+X,TH,X,Z,A,PF)

      F1    = ( xb*Tpb/(W-X) ) *phi((W-X)/(E))   ! 
      F1    = F1 + xi/(2.*(W-X)**2)
      F1    = F1 * SIG1*(1.+D2)
      F1    = F1 * ( (W-X)/((E-W)*R) )**(xb*Tpa)
      F1    = F1 * ( (W-X)/ (E)      )**(xb*Tpb) ! 

      IF(IFUNC.EQ.2) THEN      ! dEp'
        F=F2
      ELSEIF (IFUNC.EQ.1) THEN ! dEs'
        F=F1
      ENDIF
      RETURN
      END

      SUBROUTINE ROM(A,B,EPS,ANS,K,IFUNC)
      IMPLICIT double precision (A-H,O-Z)
      COMMON/ADDONS/SCALE,Tb,Ta
C  ROMBERG METHOD OF INTEGRATION
      DIMENSION W(2,50)
      H=B-A
      K=0
      CALL VALY(A,FA,IFUNC)
      CALL VALY(B,FB,IFUNC)
      W(2,1)=(FA+FB)*H/2.
    4 K=K+1
      IF(K.GE.49)GO TO 5
      DO 100 I=1,K
  100 W(1,I)=W(2,I)
      H=H/2.
      SIG=0.
      M=2**(K-1)
      DO 1 J=1,M
      J1=2*J-1
      X=A+FLOAT(J1)*H
      CALL VALY(X,F,IFUNC)
    1 SIG=SIG+F
      W(2,1)=W(1,1)/2.+H*SIG
      DO 2 L=1,K
    2 W(2,L+1)=(4.**(L)*W(2,L)-W(1,L))/(4.**(L)-1.)
      E=(W(2,K+1)-W(1,K))/W(2,K+1)
      IF(ABS(E)-EPS) 3,3,4
    3 ANS=W(2,K+1)
      RETURN
    5 ANS=W(2,K+1)
      END
      
      double precision function phi(x)
      IMPLICIT double precision (A-H,O-Z)
      phi=1.0-x+3./4.*x**2
      RETURN
      END
      
