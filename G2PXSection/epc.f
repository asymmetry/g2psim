C  ELECTROPRODUCTION YIELDS OF NUCLEONS AND PIONS 
C  WRITTEN BY J.S. O'CONNELL AND J.W. LIGHTBODY, JR.
C  NATIONAL BUREAU OF STANDARDS 
C  SEPTEMBER 1987 
! all the variables finishing by al were added to the original code
! and all the comment with ! 
! and all the non-capital text. 
! CEBAF A. Deur 07/30/98
C
C  APRIL 1988
C  TRANSVERSE SCALING REGION ADDED
!z
!z modified by X. Zheng, 12/13/2000
!z based on 1988 version of EPC and modified the input method 
!z for e99117
C
C modified by C. Gu 04/01/2011
C for a library used by Geant4 similation

      subroutine epc(PART,tgt_Z,tgt_N,E1,PTP1,THP,xs)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer err,tgt_Z,tgt_N
      double precision N
      INTEGER PART
      logical debug

      COMMON/QD/QDF 
      COMMON/DEL/IP 
      COMMON/M/Z,N
      COMMON/SP/IA

      DATA PI/3.1415926536/,AM/938.272046/
      DATA AMD/1876.123846/,AMP/139.57018/
      
      Z=dble(tgt_Z)
      N=dble(tgt_N)
      PTP=PTP1
      IA=Z+N
C  'AN' IS EFFECTIVE NUMBER OF NUCLEONS FOR PION PRODUCTION 
      IF(PART.EQ.2212)THEN 
         AN=N/3.+2.*Z/3.
         IP=1 
      ELSEIF(PART.EQ.2112)THEN 
         AN=Z/3.+2.*N/3.
         IP=-1
      ELSEIF(PART.EQ.211)THEN 
         AN=Z/3.
         IP=2 
      ELSEIF(PART.EQ.-211)THEN 
         AN=N/3.
         IP=2 
      ELSEIF(PART.EQ.111)THEN 
         AN=2.*(N+Z)/3. 
         IP=2 
      ELSE
      ENDIF
      IF(ABS(IP).EQ.1)THEN
         IF(IA.GT.1.AND.IA.LT.5)THEN
            DLF=IA
         ELSE 
            DLF=7.
         ENDIF
         AL=DLF
         QDF=AL*N*Z/FLOAT(IA) 
      ELSE
      ENDIF 
      TH=THP*PI/180.
      IF(ABS(IP).EQ.1)THEN
         E=SQRT(PTP**2+AM**2) 
         TP=PTP**2/(E+AM) 
         AJ=PTP/E 
         IF(IA.EQ.1)THEN
            D2QD=0. 
            D2QF=0. 
         ELSEIF(IA.GT.1)THEN
            CALL DEP(E1,TP,TH,IP,D2QD)
            D2QD=D2QD*AJ
            CALL EP(E1,TP,TH,D2QF)
            D2QF=D2QF*AJ
         ENDIF
      ELSEIF(ABS(IP).EQ.2.OR.IP.EQ.0)THEN
         E=SQRT(PTP**2+AMP**2)
         TP=PTP**2/(E+AMP)
         AJ=PTP/E 
         D2QD=0.
         D2QF=0.
      ELSE
      ENDIF

      CALL DELTA(E1,TP,TH,D2DEL)
      D2DEL=AJ*AN*D2DEL
      
      TOTAL=D2QD+D2QF+D2DEL
      xs=TOTAL
      
C      print *,D2QD,D2QF,D2DEL

      END 
*----------------------------------------------------------------
*VTP
      SUBROUTINE VTP(AMT,AM1,EI,W0,TP,TH,GN)
C  TIATOR-WRIGHT VIRTUAL PHOTON SPECTRUM
C  PHYS. REV. C26,2349(1982) AND NUC. PHYS. A379,407(1982)
      IMPLICIT double precision(A-H,O-Z) 
      DATA AME/.510998928/,PI/3.1415926536/ 
      EF0=EI-W0
C      print *,EF0
      AKI=SQRT(EI**2-AME**2)
      if(EF0.LT.AME)GO TO 1
      AKF0=SQRT(EF0**2-AME**2)
      AKP=SQRT(TP**2+2.*AM1*TP) 
      EP=TP+AM1 
      AR=EI+AMT-EP
      BR=EF0*(AKP*COS(TH)-AKI)/AKF0 
      BRP=(AKF0/EF0)**2*BR
      A=AME**2-EI*EF0 
      B=AKI*AKF0
      D=-AME**2*BR*(EI/EF0-1.)/AR 
      AP=A-D
      BP=B+D
      AN1=1./137./2./PI*W0**2/AKI**2
      APB=-AME**2*(AKI-AKF0)**2/(AME**2+EI*EF0+AKI*AKF0)
      AN1=AN1*B/BP*(AR+BR)/(AR-AP/BP*BR)
      AN2=1.-2.*A/W0**2 
      AN4=((AP-BP)*(AR+BR)/(AP+BP)/(AR-BR)) 
C      print *,AN1,AN2,AN3,AN4
      IF(AN4.LE.0.)GO TO 1
      AN2=AN2*LOG(AN4)
      AN3=-4.*B/W0**2 
      ANE=AN1*(AN2+AN3) 
      D0=AR+BR
      R=(AMT+W0-EP/AKP*W0*COS(TH))/D0 
      GN=ANE*R/W0
      IF(GN.LT.0.)GN=0. 
      RETURN
    1 GN=0. 
      RETURN
      END 
*DEP
      SUBROUTINE DEP(E1,TP,TH,IP,D2N) 
C  QUASI-DEUTERON CROSS SECTION 
      IMPLICIT double precision (A-H,O-Z) 
      COMMON/QD/QDF 
C      DIMENSION C(0:4,8),A(0:4) 
      DATA PI/3.1415926536/,AM/938.272046/,AMD/1876.123846/ 
      PN=SQRT(TP**2+2.*AM*TP) 
      EP=TP+AM
      CALL KINE(AMD,AM,AM,PN,TH,W0,THC) 
      IF(W0.GE.E1)GO TO 1 
      IF(W0.LE.0.)GO TO 1 
      W0G=W0/1000.
      CALL SIGD(W0G,THC,IP,DSQD)
      CALL PART(AMD,AM,AM,PN,TH,AJT,AJW)
      DSQD=AJT*DSQD 
C  CROSS SECTION IN UB/MEV-SR 
      CALL VTP(AMD,AM,E1,W0,TP,TH,PHI)
      D2N=QDF*PHI*DSQD
      RETURN
    1 D2N=0.
      RETURN
      END 
*SIGD 
      SUBROUTINE SIGD(E,TH,IP,DSQD) 
C  DEUTERON CROSS SECTION 
C  BASED ON FIT OF THORLACIUS & FEARING 
C  PHYS. REV. C33,1830(1986)
C  ENERGY RANGE 10 - 625 MEV
C 
C  E[GEV] IN LAB SYSTEM 
C  TH[RAD] & DSQD[UB/SR] IN CENTER-OF-MOMENTUM SYSTEM 
C 
      IMPLICIT double precision (A-H,O-Z) 
      DIMENSION C0(8),C1(4),C2(4),C3(4),C4(4) 
      DIMENSION A(0:4),B(4,4)
      DATA C0/2.61E2,-1.10E2,2.46E1,-1.71E1,5.76E0,-2.05E0,2.67E-1, 
     1    1.13E2/ 
      DATA C1/1.68E1,-4.66E1,2.56E0,-4.72E0/
      DATA C2/-2.03E2,-8.12E1,-4.05E0,-5.99E0/
      DATA C3/-1.77E1,-3.74E1,-5.07E-1,-5.40E0/ 
      DATA C4/-2.05E0,-7.05E0,9.40E-1,-2.05E0/
      X=COS(TH)
      IF(E.LE.0.625)THEN 
C  TEST FOR NEUTRON 
      X=IP*X
C  COEFICIENTS
      A(0)=C0(1)*EXP(C0(2)*E)+ C0(3)*EXP(C0(4)*E) 
      A(0)=A(0)+(C0(5)+C0(6)*E)/(1.+C0(8)*(E-C0(7))**2) 
      DSQD=A(0)*P(0,X)
      DO 2 L=1,4
      B(1,L)=C1(L)
      B(2,L)=C2(L)
      B(3,L)=C3(L)
    2 B(4,L)=C4(L)
      DO 1 L=1,4
      A(L)=B(L,1)*EXP(B(L,2)*E)+ B(L,3)*EXP(B(L,4)*E) 
    1 DSQD=DSQD+A(L)*P(L,X) 
      ELSEIF(E.LT.0.700)THEN 
         DSQD=.3
      ELSEIF(E.LT.0.800)THEN 
         DSQD=.15 
      ELSEIF(E.LT.0.900)THEN 
         DSQD=.1
      ELSE
         DSQD=55./(E*1000.0-350.)
      ENDIF 
      RETURN
      END 
*LEG
      double precision FUNCTION P(L,X)
C  LEGENDRE POLYNOMIALS 
      IMPLICIT double precision (A-H,O-Z) 
      IF(L.EQ.0)THEN
         P=1. 
      ELSEIF(L.EQ.1)THEN
         P=X
      ELSEIF(L.EQ.2)THEN
         P=.5*(3.*X**2-1.)
      ELSEIF(L.EQ.3)THEN
         P=.5*(5.*X**3-3.*X)
      ELSEIF(L.EQ.4)THEN
         P=1./8.*(35.*X**4-30.*X**2+3.) 
      ELSE
         P=0. 
      ENDIF 
      RETURN
      END 
*DELTA
      SUBROUTINE DELTA(E1,TP,TH,D2DEL)
C  PHOTOPRODUCTION OF NUCLEONS AND PIONS VIA DELTA
      IMPLICIT double precision (A-H,O-Z) 
      COMMON/DEL/IP 
      DATA PI/3.1415926536/,AM/938.272046/,AMP/139.57018/
      IF(ABS(IP).EQ.1)THEN
         AM1=AM 
         AM2=AMP
      ELSE
         AM1=AMP
         AM2=AM 
      ENDIF 
      EP=TP+AM1 
      PN=SQRT(EP**2-AM1**2) 
      CALL KINE(AM,AM1,AM2,PN,TH,W,TC)
C      PRINT *, W
      IF(W.LE.0.)GO TO 1
      IF(W.GE.E1)GO TO 1
      CALL PART(AM,AM1,AM2,PN,TH,AJT,AJW) 
C      print *,AJT
      CALL SIGMA(W,TC,DSIGG)
C      print *,W,DSIGG
      CALL VTP(AM,AM1,E1,W,TP,TH,PHI)
C      print *,PHI
      D2DEL=PHI*DSIGG*AJT 
C CROSS SECTION IN UB/MEV-SR
      RETURN
    1 D2DEL=0.
      RETURN
      END 
*PART 
      SUBROUTINE PART(AMT,AM1,AM2,PN,TN,AJT,AJW)
      IMPLICIT double precision (A-H,O-Z) 
C  PARTIAL DERIVATIVES
      DATA PI/3.1415926536/ 
      DT=PI/50. 
      DP=10.
C  ANGLE
      TNP=TN+DT 
      TNM=TN-DT 
      CALL KINE(AMT,AM1,AM2,PN,TNP,W,TCP) 
      CALL KINE(AMT,AM1,AM2,PN,TNM,W,TCM) 
      AJT=(COS(TCP)-COS(TCM))/(COS(TNP)-COS(TNM)) 
      AJT=ABS(AJT)
C  ENERGY 
      PNP=PN+DP 
      PNM=PN-DP 
      CALL KINE(AMT,AM1,AM2,PNP,TN,WP,TC) 
      CALL KINE(AMT,AM1,AM2,PNM,TN,WM,TC) 
      AJW=(WP-WM)/(PNP-PNM) 
      AJW=ABS(AJW)
      RETURN
      END 
*KINE 
      SUBROUTINE KINE(AMT,AM1,AM2,PN,TH,W,TC) 
C  COMPUTES CM VARIABLES FROM LAB VARIABLES 
      IMPLICIT double precision (A-H,O-Z) 
      EP=SQRT(PN**2+AM1**2) 
      PNT=PN*SIN(TH)
      PNL=PN*COS(TH)
      ANUM=PN**2+AM2**2-(AMT-EP)**2 
      DEN=2.*(PNL+AMT-EP) 
      W=ANUM/DEN
C     PRINT *,W,ANUM, PN, AMT, EP, DEN
      IF(W.LE.0.)W=0. 
C  INVARIANT MASS 
      WW=SQRT(AMT**2+2.*W*AMT)
C  CM VARIABLES 
      PCT=PNT 
      B=W/(AMT+W) 
      G=(W+AMT)/WW
      PCL=G*(PNL-B*EP)
      PCS=PCL**2+PCT**2 
      PC=SQRT(PCS)
      CTHC=PCL/PC 
      TC=ACOS(CTHC) 
      RETURN
      END 
*SIGMA
      SUBROUTINE SIGMA(E,THRCM,SIGCM) 
      IMPLICIT double precision (A-H,O-Z) 
      GAM=100.
      PI=3.1415926536
      IF(E.GT.420.)THEN
        SIGCM=90./4./PI/(1.+420./E)
      ELSE
        SIGCM=360.*(5.-3.*COS(THRCM)**2)/16./PI/(1.+(E-320.)**2/GAM**2) 
      ENDIF 
      RETURN
      END 
C  BEGIN QUASI-FREE SCATTERING CROSS SECTION
*EP 
      SUBROUTINE EP(E1,TP,THP,DSEP) 
C  ELECTRO PROTON PRODUCTION CROSS SECTIONS 
      IMPLICIT double precision (A-H,O-Z) 
      COMMON/Pcm/PH(10),WPH(10) 
      DATA AML/.510998928/,PI/3.1415926536/ 
      CALL GAUSAB(10,PH,WPH,0.D0,2.*PI,PI)
      AK=SQRT(E1**2-AML**2) 
      CALL SEP(AK,TP,THP,DSEP)
      DSEP=DSEP*1.E4
C  CROSS SECTION IN UB/MEV-SR 
      END 
*DOT
      double precision FUNCTION  DOT(V,U) 
      IMPLICIT double precision (A-H,O-Z) 
      DIMENSION V(3),U(3) 
      DOT=0.
      DO 1 I=1,3
    1 DOT=DOT+V(I)*U(I) 
      RETURN
      END 
*CROSS
      SUBROUTINE CROSS(V,U,W) 
      IMPLICIT double precision (A-H,O-Z) 
      DIMENSION V(3),U(3),W(3)
      W(1)=V(2)*U(3)-V(3)*U(2)
      W(2)=V(3)*U(1)-V(1)*U(3)
      W(3)=V(1)*U(2)-V(2)*U(1)
      RETURN
      END 
*GAUSAB 
C   SUBROUTINE GAUSAB 
C 
      SUBROUTINE GAUSAB(N,E,W,A,B,C)
      IMPLICIT double precision (A-H,O-Z) 
      DIMENSION E(24),W(24) 
      DATA PI/3.141592653589793238462643D0/,EPS/1.D-16/ 
      AL=(C*(A+B)-2*A*B)/(B-A)
      BE=(A+B-2*C)/(B-A)
      M=(N+1)/2 
      DN=N
      DO 5 I=1,M
       DI=I 
       X=PI*(4.D0*(DN-DI)+3.D0)/(4.D0*DN+2.D0)
       XN=(1.D0-(DN-1.D0)/(8.D0*DN*DN*DN))*COS(X) 
       IF(I.GT.N/2) XN=0
       DO 3 ITER=1,10 
        X=XN
        Y1=1.D0 
        Y=X 
        IF(N.LT.2) GO TO 2
        DO 1 J=2,N
         DJ=J 
         Y2=Y1
         Y1=Y 
    1    Y=((2.D0*DJ-1.D0)*X*Y1-(DJ-1.D0)*Y2)/DJ
    2   CONTINUE
        YS=DN*(X*Y-Y1)/(X*X-1.D0) 
        H=-Y/YS 
        XN=X+H
        IF(ABS(H).LT.EPS) GO TO 4 
    3   CONTINUE
    4  E(I)=(C+AL*X)/(1.D0-BE*X)
       E(N-I+1)=(C-AL*X)/(1.D0+BE*X)
       GEW=2.D0/((1.D0-X*X)*YS*YS)
       W(I)=GEW*(AL+BE*C)/(1.D0-BE*X)**2
       W(N-I+1)=GEW*(AL+BE*C)/(1.D0+BE*X)**2
    5  CONTINUE 
      RETURN
      END 
*VECT 
      SUBROUTINE VECT(THP,THE,PHI,P,AK1,AK2)
C  CARTESIAN COMPONENTS OF ELECTRON AND PROTON VECTORS
      IMPLICIT double precision (A-H,O-Z) 
      COMMON/V/AK1V(3),AK2V(3),QV(3),PV(3),PP(3)
      PV(1)=P*SIN(THP)
      PV(2)=0.
      PV(3)=P*COS(THP)
      AK1V(1)=0.
      AK1V(2)=0.
      AK1V(3)=AK1 
      AK2V(1)=AK2*SIN(THE)*COS(PHI) 
      AK2V(2)=AK2*SIN(THE)*SIN(PHI) 
      AK2V(3)=AK2*COS(THE)
      QV(1)=AK1V(1)-AK2V(1) 
      QV(2)=AK1V(2)-AK2V(2) 
      QV(3)=AK1V(3)-AK2V(3) 
      PP(1)=PV(1)-QV(1) 
      PP(2)=PV(2)-QV(2) 
      PP(3)=PV(3)-QV(3) 
      RETURN
      END 
*AMAG 
      double precision FUNCTION AMAG(V) 
      IMPLICIT double precision (A-H,O-Z) 
      DIMENSION V(3)
      AMAG=0. 
      DO 1 I=1,3
    1 AMAG=AMAG+V(I)**2 
      AMAG=SQRT(AMAG) 
      RETURN
      END 
*LEPT 
      SUBROUTINE LEPT(E1,E2,AK1,AK2,AML,QS,QUS,THE,V) 
C  LEPTON FACTORS FOR COINCIDENCE CROSS SECTION 
      IMPLICIT double precision (A-H,O-Z) 
      DIMENSION V(5)
      V(1)=(QUS/QS)**2*(E1*E2+AK1*AK2*COS(THE)+AML**2)
      X=AK1*AK2*SIN(THE)
      V(2)=X**2/QS+QUS/2. 
      V(3)=QUS/QS*X/SQRT(QS)*(E1+E2)
      V(4)=X**2/QS
      V(5)=0. 
      RETURN
      END 
*D4S
      SUBROUTINE D4S(AK1,AK2,THE,P,PP,THQP,CPHIP,DSIG)
C  FULLY DIFFERENTIAL CROSS SECTION 
      IMPLICIT double precision (A-H,O-Z) 
      DIMENSION V(5),W(5) 
      DATA AM/938.272046/,AML/.510998928/,PI/3.1415926536/,A/855./
      QS=AK1**2+AK2**2-2.*AK1*AK2*COS(THE)
      E1=SQRT(AK1**2+AML**2)
      E2=SQRT(AK2**2+AML**2)
      QUS=2.*(E1*E2-AK1*AK2*COS(THE)-AML**2)
      SM=2.*(1.44)**2/QUS**2*AK2/AK1
      EP=SQRT(AM**2+P**2) 
      PS=EP*P 
      FNS=1./(1.+QUS/A**2)**4 
      CALL LEPT(E1,E2,AK1,AK2,AML,QS,QUS,THE,V) 
      CALL FORM(QS,QUS,P,THQP,CPHIP,W)
      SUM=0.
      DO 1 I=1,5
    1 SUM=SUM+V(I)*W(I) 
      DSIG=SM*PS*FNS*SUM*SGSL(PP) 
      RETURN
      END 
*STHE 
      SUBROUTINE STHE(D2S)
C  INTEGRAL OVER ELECTRON POLAR ANGLE 
      IMPLICIT double precision (A-H,O-Z) 
      COMMON/S/ AK1,AK2,THE,P,THP 
      COMMON/E/TH1(12),WT1(12),TH2(12),WT2(12)
      COMMON/E1/TH3(24),WT3(24) 
      D2S1=0. 
      DO 1 I=1,12 
      THE=TH1(I)
      CALL SPHI(D3S)
    1 D2S1=D2S1+D3S*WT1(I)*SIN(THE) 
      D2S2=0. 
      DO 2 I=1,12 
      THE=TH2(I)
      CALL SPHI(D3S)
    2 D2S2=D2S2+D3S*WT2(I)*SIN(THE) 
      D2S3=0. 
      DO 3 I=1,24 
      THE=TH3(I)
      CALL SPHI(D3S)
    3 D2S3=D2S3+D3S*WT3(I)*SIN(THE) 
      D2S=D2S1+D2S2+D2S3
      RETURN
      END 
*SPHI 
      SUBROUTINE SPHI(D3S)
C  INTEGRATE OVER ELECTRON AZIMUTHAL ANGLE
      IMPLICIT double precision (A-H,O-Z) 
      COMMON/S/ AK1,AK2,THE,P,THP 
      COMMON/V/AK1V(3),AK2V(3),QV(3),PV(3),PP(3)
      COMMON/Pcm/PH(10),WPH(10) 
      DIMENSION QXP(3),AK1X2(3) 
      D3S=0.
      DO 1 I=1,10 
      PHI=PH(I) 
      CALL VECT(THP,THE,PHI,P,AK1,AK2)
      CALL CROSS(QV,PV,QXP) 
      CALL CROSS(AK1V,AK2V,AK1X2) 
C  PROTON THETA 
      CTHEP=DOT(PV,QV)/AMAG(PV)/AMAG(QV)
      THQP=ACOS(CTHEP)
C  PROTON PHI 
      CPHIP=DOT(QXP,AK1X2)/AMAG(QXP)/AMAG(AK1X2)
      PPM=AMAG(PP)
      CALL D4S(AK1,AK2,THE,P,PPM,THQP,CPHIP,DSIG) 
    1 D3S=D3S+DSIG*WPH(I) 
      RETURN
      END 
*FORM 
      SUBROUTINE FORM(QS,QUS,P,THQP,CPHIP,W)
C  NUCLEAR FORM FACTORS 
      IMPLICIT double precision (A-H,O-Z)
      double precision NN
      COMMON/M/ZZ,NN
      COMMON/DEL/IP 
      DIMENSION W(5)
      DATA AM/938.272046/,UP/2.79/,UN/-1.91/
      IF(IP.EQ.1)THEN 
         Z=ZZ 
         N=0. 
      ELSEIF(IP.EQ.-1)THEN
         Z=0. 
         N=NN 
      ELSE
         Z=0. 
         N=0. 
      ENDIF 
      Y=P/AM*SIN(THQP)
      W(1)=Z
      W(2)=Z*Y**2 
      W(2)=W(2)+(Z*UP**2+N*UN**2)*QS/2./AM**2 
      W(3)=-2.*Z*Y*CPHIP
      W(4)=Z*Y**2*(2.*CPHIP**2-1.)
      W(5)=0. 
      RETURN
      END 
*SEP
      SUBROUTINE SEP(AK,TP,THPP,D2S)
      IMPLICIT double precision (A-H,O-Z) 
      COMMON/S/AK1,AK2,THE,P,THP
      COMMON/E/TH1(12),WT1(12),TH2(12),WT2(12)
      COMMON/E1/TH3(24),WT3(24) 
      DATA AM/938.272046/,AML/.510998928/,BE/16./
      DATA PI/3.1415926536/
      THP=THPP
      AK1=AK
      AK2=AK1-TP-BE 
C  GAUSSIAN POINTS FOR THE
      THEMAX=AML*(AK1-AK2)/AK1/AK2
      CALL GAUSAB(12,TH1,WT1,0.D0,2.*THEMAX,THEMAX) 
      CALL GAUSAB(12,TH2,WT2,2.*THEMAX,100.*THEMAX,10.*THEMAX)
      A3=100.*THEMAX
      C3=A3+(PI-A3)/10. 
      CALL GAUSAB(24,TH3,WT3,A3,PI,A3+(PI-A3)/10.)
      P=SQRT(TP**2+2.*AM*TP)
      CALL STHE(D2S)
      IF(AK2.LE.0.)D2S=0. 
      RETURN
      END 
*SGSL 
      double precision FUNCTION SGSL(P) 
C  P INTEGRAL OVER SGSL NORMALIZED TO 1/4PI 
      IMPLICIT double precision (A-H,O-Z) 
      COMMON/SP/IA
      IF(IA.EQ.2)THEN 
C  BEGIN 2-H
         PP=P/197.3 
         SGS=3.697-7.428*PP-2.257*PP**2 
         SGS=SGS+3.618*PP**3-1.377*PP**4+.221*PP**5-.013*PP**6
         IF(SGS.LT.-293.)GO TO 1
         SGS=EXP(SGS) 
         SGS=SGS/.18825/4./3.1416/(197.3)**3
         SGSL=SGS/1.
      ELSEIF(IA.EQ.3)THEN 
C  BEGIN 3-HE 
         IF(-(P/33)**2.LT.-293.)GO TO 1 
         SGS=2.4101E-6*EXP(-P/33) 
         SGS=SGS-1.4461E-6*EXP(-(P/33)**2)
         SGS=SGS+1.6871E-10*EXP(-(P/493)**2)
         SGSL=SGS/2.
      ELSEIF(IA.EQ.4)THEN 
C   BEGIN 4-HE
         IF(-(P/113.24)**2.LT.-293.)GO TO 1 
         SGS=1.39066E-6*EXP(-(P/113.24)**2) 
         SGS=SGS+3.96476E-9*EXP(-(P/390.75)**2) 
         SGSL=SGS/2.
         SGSL=SGSL/2./3.1416
      ELSEIF(IA.EQ.12)THEN
C  BEGIN 12-C 
         IF(-(P/127)**2.LT.-293.)GO TO 1
         SGS=1.7052E-7*(1.+(P/127)**2)*EXP(-(P/127)**2) 
         SGS=SGS+1.7052E-9*EXP(-(P/493)**2) 
         SGSL=SGS/6.
      ELSE
C  BEGIN 16-O 
         IF(-(P/120)**2.LT.-293.)GO TO 1
         SGS=3.0124E-7*(1.+(P/120)**2)*EXP(-(P/120)**2) 
         SGS=SGS+1.1296E-9*EXP(-(P/493)**2) 
         SGSL=SGS/8.
      ENDIF 
      RETURN
    1 SGSL=0. 
      RETURN
      END 
