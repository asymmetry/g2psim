      SUBROUTINE RAVPRO(GAMP,ELAB,PCM,ATWT,CONV,X0,STEP,OTESTD,
     1 THELAB,THECM,SIGLAB,SIGCM,NANG,NF,NFL,DELTH,DELPH,XMS,NCM,RSR)
C
C    - * -  R A V E N H A L L  P R O G R A M  - * -
C
C    This version which runs under Unix was generated from the VAX
C    VMS version which came from the  Perkin-Elmer version used
C    at Illinois.  That version came from the NORD-500 version used
C    at Saclay, which in turn came from a conversion of the Cyber
C    175 version used at Illinois.
C
C    In version 1B, the standard DEQ calls have been used.  The DEQ
C    routine has been removed (its in the MATH library).
C
C    In version 1C, all stop statements have been modified to print
C    more informative messages.  There have been some minor cosmetic
C    mods to the comments.  Finally, the call to the GAMMA function
C    routine has been modified in anticipation of the inclusion of
C    a scale factor to avoid exponent overflows.  The GAMMA function
C    routine has not been changed in this version however.  This
C    version used to be called RAVPRO-TEST.
C
C    In version 1D, the size of the array POT is now set by the parameter
C    NPOTMX.  The sizes of the arrays: TND, FRN, FIN, GRN, and GIN are
C    set by the parameter NFMX.  These parameters are defined in the
C    included file ravpro.inc.
C
C    In version 1E, the GAMMA function has been renamed to CDGAM to
C    avoid a conflict with the intrinsic GAMMA function included in
C    some compiler's versions of Fortran.
C
C    Input specifications:
C
C    Input data contained in the subroutine call:
C          GAMP    Atomic number times the fine structure constant
C          ELAB    Total lab energy of electron (MeV). Used if NCM
C            is 0.
C          PCM     C/M momentum of electron (MeV/c). Used if NCM is 1.
C            Calculated if NCM is 0.
C          ATWT    Atomic weight of target
C          CONV    Absolute magnitude of phase shifts used to halt
C            calculation of subsequent phase shifts.
C          X0      Fitting-on radius
C          STEP    Step size for integration of Dirac equations
C          OTEST   Maximum angle for which small angle calculation
C            will be done (in degrees).
C          THELAB  Lab angles (in degrees). Used if NCM is 0.
C          THECM   C/M angles (in degrees). Used if NCM is 1,
C            calculated if NCM is 0.
C          NANG    Number of angles for which cross sections should
C            be calculated.
C          NF      Maximum number of phase shifts
C          NFL     Number of terms in Legendre series.
C          DELTH   The half-angular resolution in the scattering
C            plane.  Used for folding when NCM is 0.
C          DELPH   The half-angular resolution perpendicular to the
C            scattering plane.  Used for folding when NCM is 0.
C          XMS     If positive, XMS is the RMS multiple scattering
C            angle in degrees.  If negative, XMS is the target
C            thickness in grams per square cm and the multiple
C            scattering angle is calculated using a formulae from
C            Larry Cardman's thesis (appendix B).  If XMS is
C            0, no multiple scattering correction is made.
C            Multiple scattering corrections are made only if
C            NCM is 0.
C          NCM     If NCM is 0, make all necessary calculations for the
C            transformation from the C/M to the lab frame.
C            Cross sections are calculated for nice lab angles.
C            If NCM is non-zero, no transformations are made.
C            Cross sections are calculated at nice C/M angles.
C
C    Input data in COMMON /COUL/ (set by COULS, which must be
C          called before RAVPRO using the same X0, NF, and GAMMA)
C          FRN     F-regular
C          GRN     G-regular
C          FIN     F-irregular
C          GIN     G-irregular
C
C    Input data in COMMON /IO/
C          NPRINT causes the following data to be printed when set:
C            1.  COULS coulomb function subroutine input
C            2.  RAVPRO input
C            3.  Charge density information (RHO, V, NORM, ...)
C            4.  Dirac equation solutions
C            5.  Phase shifts
C            6.  Energy, angles, cross section in C/M,
C            form factors.
C            7.  COULS integration results
C            8.  Intermediate results of integration and angle
C            dependant calculation
C            9.  Intermediate integration results
C            10. Folding integration results.
C          NPUNCH causes the following data to be punched when set:
C            1.  Coulomb phase shifts
C            2.  Amplitudes
C            3.  Cross sections and angles
C            4.  Unused
C            5.  Unused
C
C    Input data in COMMON /CONST/
C          XMN     One atomic mass unit (AMU) in MeV/C**2
C          HBARC   H-bar times c in MeV-Fermies.
C          Alpha   The fine structure constant
C          XME     The mass of the electron in MeV/C**2
C          PI      The mathematical constant PI.
C          DTR     PI/180, conversion factor for degrees to radians
C          RAD     180/PI, conversion factor for radians to degrees
C
C    Input data in COMMON /CHARGE/
C          MODELN  The charge density model number
C          PARAM   The parameters of the charge density.
C            Distance units should be in Fermies.
C
C    Data returned:
C
C          SIGLAB  Lab cross sections at specified angles (corrected
C            for finite mass of the electron)
C          SIGCM   C/M cross sections at specified angles.
C          RSR     RMS radius calculated from charge density
C          THCM    C/M angles calculated if NCM is 0.
C          PCM     C/M momentum calculated if NCM is 0
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Include parameter file to define array sizes
      INCLUDE 'ravpro.inc'
      EXTERNAL RAVDEQ
      CHARACTER RHONA*72
      PARAMETER (NPMX=29)
      DOUBLE PRECISION ZV10(4), ZV11(4), ZV12(3), ZV13(4)
      DOUBLE PRECISION ZV14(11), ZV15(10)
      INTEGER NZ001(10), NZ002(11)
      DOUBLE PRECISION TND(NFMX),FRN(NFMX),FIN(NFMX),GRN(NFMX),GIN(NFMX)
      DOUBLE PRECISION P(101)
      DOUBLE PRECISION S(3)
      DOUBLE PRECISION THELAB(NANG),THECM(NANG),SIGLAB(NANG)
      DOUBLE PRECISION SIGCM(NANG),DELTH(NANG),DELPH(NANG),PARAM(NPMX)
      DOUBLE PRECISION XMS(NANG)
      INTEGER NPRINT(10), NPUNCH(5)
      DOUBLE COMPLEX AA,BB,TT,TTT,BETA,BETA0,ETA,ETN,PS,PSN,PSNM
      DOUBLE COMPLEX SUM,FO,CFO,GAMC,GAMM,A(101),B(101)
      DOUBLE COMPLEX CDGAM,CDEXPI,CDIVIDE
      EQUIVALENCE (X01,ZV10(1)),(X02,ZV10(2)),(X03,ZV10(3))
      EQUIVALENCE (X04,ZV10(4)),(DX01,ZV11(1)),(DX02,ZV11(2))
      EQUIVALENCE (DX03,ZV11(3)),(DX04,ZV11(4)),(XX,ZV12(1))
      COMMON /DEQ/ ZV10,ZV11,ZV12,ZV13,ZV14,ZV15,
     1 GAM,FN,FFI,FFJ,FUN,CD,XDEQ,MDEQ,N,NST,FK
      COMMON /CHARGE/ MODELN,PARAM
      COMMON /IO/ NPRINT,NPUNCH
      COMMON /COUL/ FRN,GRN,FIN,GIN
      COMMON /CONST/ XMN,HBARC,ALPHA,XME,PI,DTR,RAD
      COMMON /POT/ POT(NPOTMX)
C
      NST=1
      NS=0
      AMEV=ATWT*XMN
      GAM=GAMP
      GAMC=DCMPLX(0.0D0,GAM)
      ZV15(1)=STEP
      IF (NCM.NE.0) GO TO 7701
      PCM=DSQRT(AMEV/(AMEV+2.0D0*ELAB))
      ECMN=PCM*(AMEV+ELAB)
      PCM=DSQRT((ELAB*PCM)**2-XME**2)
7701      FK=PCM/HBARC
C
C    Set up initial default conditions
C    Integration parameters
      ZV15(2)=0.00D0
      NZ002(1)=49
      NZ002(2)=0
      ZV14(1)=0.0D0
      ZV14(2)=X0
      IZ2=2
C
C    Read in integration variables
C          IZ2     Number of regions
C          ZV14(i) Starting value for region i
C          ZV15(i) Step size for region i
C          NZ002(i)Print frequency for region i
C    In order to read in the integration parameters, remove the C
C    in column 1 of the following lines:
C    READ (7,9050) IZ2,(ZV14(IZ1), ZV15(IZ1),NZ002(IZ1),IZ1 = 1,IZ2)
C9050     FORMAT (1I7/(2D16.8,I7))
      IZF=IZ2-1
      DO 4449 I=1,IZF
      IF (DABS(ZV15(I))-1.0D-08) 5551,5551,4449
5551      WRITE (6,5553)
5553     FORMAT (' Zero step error')
      STOP 'RAVPRO aborted -- zero step error'
4449      CONTINUE
C
C    Write out the input
      IF (NPRINT(2).NE.1) GO TO 7709
      WRITE (6,8012)
8012     FORMAT ('1',T36,'- * - R A V E N H A L L  P R O G R A M - * -')
      WRITE(*,*) NPRINT(2)
      CALL CHGDE3(NPAR)
      WRITE (6,7703) RHONA
7703     FORMAT (' Charge distribution is ',A)
      WRITE (6,7704) (I,PARAM(I),I=1,NPAR)
7704     FORMAT (' Charge distribution parameter',I3,5X,1PE23.16)
      IF (NCM) 7705,7705,7707
7705      WRITE (6,7706) ELAB,PCM,ECMN,ATWT,FK
7706     FORMAT ('0Transformation from lab angles to C/M angles',
     1 ' will be made.  The cross sections will be found'/
     2 ' in the C/M system and transformed to the lab system'/
     3 ' The energy of the electron (positron) in the lab system = ',
     4 1PE23.16,' MeV'/
     5 ' The momentum of the electron (positron) in the C/M ',
     6 'system = ',E23.16,' MeV/c'/
     7 ' The energy of the target in the C/M system =',E23.16,
     8 ' MeV'/
     9 ' The atomic weight of the target = ',0PF13.8/
     1 ' K = ',1PE23.16)
      GO TO 7709
7707      WRITE (6,7708) PCM,FK
7708     FORMAT ('0Transformations from lab to C/M will not be made'/
     1 ' The momentum of the electron (positron) in the C/M',
     2 'system = ',1PE23.16,' MeV/c'/
     3 ' K = ',E23.16)
7709      IF (NPRINT(2).NE.1) GO TO 8015
      WRITE (6,8014) GAM,ZV14(1),ZV14(2),NF,X0
8014     FORMAT (' GAMMA = ',1PE23.16/
     1 ' Integrate Dirac equations from ',E23.16,' to ',E23.16/
     2 1X,I3,' Phase shifts will be calculated at a fitting on ',
     3 'radius of ',E23.16)
8015      GAMM=DCONJG(GAMC)
      MDEQ=0
      IF (NPRINT(2).EQ.1) WRITE (6,5554) (ZV14(IZ1),ZV15(IZ1),IZ1=1,IZ2)
5554     FORMAT (' Integration variables:: X = ',1PE23.16,
     1 ', DELTAX =',E23.16)
C
C    Begin calculation of the nuclear phase shifts
      DO 3240 N=NST,NF
      XDEQ=ZV15(1)
      ZV14(1)=XDEQ
      IF (N.GT.NST) GO TO 20
      CALL CHGDE1 (FK,F1,F2,F3)
      FIG=F1/3.D0 + F2*XDEQ/4.D0 + F3*XDEQ**2/10.D0
      FJG=F1/2.D0 + F2*XDEQ/3.D0 + F3*XDEQ**2/8.D0
      V22=2.D0*(FIG-FJG)
      FIG=FIG*XDEQ**3
      FJG=FJG*XDEQ**2
      X01=FIG
      X02=FJG
      X03=XDEQ**5*(F1/5.0D0+(F2/6.0D0)*XDEQ+F3*XDEQ**2/14.0D0)
      NE=3
      GO TO 9030
C
C    Calculation of RHO
20      FN=N
      RHO=DSQRT(FN*FN-GAM*GAM)
      XDEQ=ZV15(1)
      NS=INT(FN*(10.0D0**(-11.0D0/FN)))
      IF (NS) 5015,5015,5012
C
C    Integration starting further away from the origin, XDEQ.
5012      STRT=NS
      ZV14(1)=STRT
      GO TO 2117
5015      ZV14(1)=XDEQ
2117      IF (NS) 5212,5212,5213
5212      V0=-GAM*FFJ/FFI
      V1=0.0D0
      V2=-GAM*V22/FFI
      B0=(1.0D0-V0)/(2.0D0*FN+1.0D0)
      A1=0.0D0
      B1=-V1/(2.0D0*FN+2.0D0)
      A2=(-(1.0D0-V0)/2.0D0)*B0
      B2=(1.0D0/(2.0D0*FN+3.0D0))*((1.0D0-V0)*A2-0.5D0*V2)
      A3=(-(1.0D0-V0)*B1+V1*B0)/3.0D0
      TF=XDEQ**MIN0(N,5)
      AF=(((A3*XDEQ+A2)*XDEQ+A1)*XDEQ+1.0D0)*TF
      BF=((B2*XDEQ+B1)*XDEQ+B0)*XDEQ*TF
      IF (N-NST) 9148,9148,9149
C
C    Initial values for the integration
9148      NE=4
      X01=FIG
      X02=FJG
      X03=AF
      X04=BF
      MDEQ=2
      GO TO 9030
C
C    Set up for integration with N greater than NST.
9149      NE=2
      X01=AF
      X02=BF
      MDEQ=2
      GO TO 9030
C
C    Starting at X not equal to XDEQ
5213      I=INT((STRT*2.0D0/XDEQ)+.51D0)
      IF (I.GT.NPOTMX) THEN
          WRITE (6,5217) I,NPOTMX
5217     FORMAT (' The array dimension for the POT array is ',
     1    'too small.'/
     2    ' An array dimension of at least',I8,' is required ',
     3    ' whereas POT is dimensioned at only',I8,'.')
          STOP 'RAVPRO aborted -- POT array index too big'
      ENDIF
      IF(NST-N)5215,5214,5214
5214      POT(I)=CD
      GO TO 5216
5215      CD=POT(I)
5216      V0=CD
C
C    Use first term only
      B0=(1.0D0-V0)/(2.0D0*FN+1.0D0)
      X01=1.0D-5
      X02=STRT*B0*1.0D-5
      MDEQ=2
      NE=2
C
9030      XX=ZV14(1)
      ZV12(2)=ZV15(1)
      ZV12(3)=ZV14(2)
      NZ001(9)=NZ002(1)
C
C    Call differential equations routine
C          ZV12    Independent variable array
C          ZV10    Array of solutions (see equivalence statements)
C          ZV11    Array of derivatives (see equivalence statements)
C          ZV13    Array used by DEQ
C          NE      The number of equations
C          X01     Solutions
C          X02     Solutions
C          X03     Solutions
C          X04     Solutions
C          DX01    Derivatives
C          DX02    Derivatives
C          DX03    Derivatives
C          DX04    Derivatives
C          XX      Independent variable, not suitabely named
      IF (ZV12(1).LT.ZV12(3)) GO TO 5100
      WRITE (6,5101) ZV12(1),ZV12(3)
5101     FORMAT (' Bad DEQ call in RAVPRO:',2(1X,1PE23.16))
      STOP 'RAVPRO aborted -- bad DEQ call'
5100      CALL INDEQ2(NZ001(9))
      NZ005=1
9011      CALL DEQ2(IRET,ZV12,ZV10,ZV11,ZV13,NE,RAVDEQ)
      IF (IRET.EQ.0) GO TO 9013
      IF (MDEQ.EQ.0) GO TO 6011
      IF (N-NST) 8117,8117,8118
8117      IF (XX.GT.XDEQ) GO TO 4117
      X=0.0D0
      CALL CHGDE2 (FK,X,FUN)
      CD=-GAM*FFJ/FFI
      IF(NPRINT(3).NE.0)WRITE(8,8119)X,CD,FUN
8119     FORMAT ('0Radius = ',1PE23.16,', potential = ',E23.16,
     1 ', charge distribution = ',E23.16)
      GO TO 6011
4117      IF (NPRINT(3).EQ.1) WRITE (6,8119) XX,CD,FUN
      NZ001(9)=50
      GO TO 6011
8118      IF(NPRINT(9).NE.0)WRITE(8,9119)XX,DX01,DX02,DX03,DX04,X01,X02,
     1 X03,X04,CD,FUN
9119     FORMAT ('0XX = ',1PE23.16,', DX01 = ',E23.16,', DX02= ',E23.16/
     1 ' DX03 = ',E23.16,', DX04 = ',E23.16,', X01 = ',E23.16/
     2 ' X02 = ',E23.16,', X03 = ',E23.16,', X04 = ',E23.16/
     3 ' POT= ',E23.16,', CD = ',E23.16)
6011      IF (N-NST) 5921,5921,5923
5921      IF (DABS(X03)-1.0D+10) 9011,8411,8411
C
C    Scaling of solutions
8411      X03=X03*1.0D-10
      X04=X04*1.0D-10
      GO TO 5924
5923      IF (DABS(X01)-1.0D+10) 9011,5925,5925
C
C    Scaling of solutions
5925      X01=X01*1.0D-10
      X02=X02*1.0D-10
5924      ZV14(1)=XX
      ZV15(1)=XDEQ
      NZ001(9)=NZ002(1)
      GO TO 9030
C
C    Re-initialization
9013      NZ005=NZ005+1
      ZV12(2)=ZV15(NZ005)
      ZV12(3)=ZV14(NZ005+1)
      NZ001(9)=NZ002(NZ005)
      IF (ZV15(NZ005)) 9011,9035,9011
9035      IF (MDEQ) 9134,9134,3030
9134      MDEQ=1
      FFI=X01
      FFJ=X02
      FFR=X03
      RSR=DSQRT(FFR/FFI)/FK
      IF (NPRINT(3).EQ.1) WRITE (6,3134) X01,X02,RSR
3134     FORMAT ('0Normalization: FFI = ',1PE23.16,5X,'FFJ = ',E23.16/
     1 ' Mean square radius = ',E23.16,' fm')
      GO TO 20
C
C    Exit from DEQ (results, G=X01 and F=X02)
C    The results of Coulomb function are in FRN, FIN, GRN, and GIN.
3030      CONTINUE
      IF (N-NST) 5031,5031,3032
C    Ratio of A to B
3032      AB=(FIN(N)/GIN(N)-X02/X01)*(GIN(N)/GRN(N))
      C1=X02/X01
      C2=FRN(N)/GRN(N)
      GO TO 3033
5031      AB=(FIN(N)/GIN(N)-X04/X03)*(GIN(N)/GRN(N))
      C1=X04/X03
      C2=FRN(N)/GRN(N)
3033      CONTINUE
C
C    Calculation of phase shifts and calculation of exponential of
C    sine and cosine argument.  Results in AA.
      A1=PI*(RHO-FN)
      A2=DTAN(A1)/DTANH(PI*GAM)
      AA=DCMPLX(1.0D0,A2)/DSQRT(1.0D0+A2*A2)*CDEXPI(A1)
      IF(NPRINT(8).NE.0)WRITE(8,3170)AA
3170     FORMAT ('0Exponential',4(1X,1PE23.16))
C
C    Calculation of tangent
      TND(N)=DIMAG(AA)*(C1-C2)/(DBLE(AA)*(C1-C2)+AB)
C    To prevent spurious convergence when Tandel changes sign
      IF (N.EQ.NST) GO TO 3190
      IF (DSQRT(TND(N)**2+TND(N-1)**2)-CONV) 3185,3185,3190
3185      NLST=N+1
      DO 3186 ITAN=NLST,NF
      TND(ITAN)=0.0D0
3186      CONTINUE
      GO TO 8888
3190      IF (NPRINT(4).EQ.0) GO TO 3240
      IF (N-NST) 3201,3201,3204
3201      WRITE (6,8007)
8007     FORMAT (1H1)
      WRITE (6,3202)
3202     FORMAT ('0',T31,'Solutions to the Dirac equations')
      WRITE (6,3203)
3203     FORMAT ('0',T11,'N',T21,'G solution',T47,'F solution')
3204      IF(NE-3)3210,3216,3216
3210      WRITE (6,3215) N,X01,X02
3215     FORMAT ('0',6X,I5,6X,1PE23.16,6X,E23.16)
      GO TO 3240
3216      WRITE (6,3215) N,X03,X04
3240      CONTINUE
8888      IF(NPRINT(5).EQ.0)GO TO 8245
      WRITE (6,8007)
      WRITE (6,8240)
8240     FORMAT ('0Phase shifts Ravenhall program'/'0',T13,'N',T32,
     1 'Phase shifts')
      WRITE(8,8242)(I,TND(I),I=1,NF)
8242     FORMAT ('0',7X,I5,7X,1PE23.16)
C
C    Calculation of the phase shifts completed.  Now calculate
C    exp(2IBETA-ZERO)
8245      BETA0=-CDGAM(GAMM,0,NTEST,0)/CDGAM(GAMC,0,NTEST,0)
      IF(NPRINT(8).NE.0)WRITE(8,570)BETA0
570     FORMAT ('0BETA0 =',2(1X,1PE23.16))
      FNS=FN
      NLAT=1
      NREC=1
      LEAF=0
      OTEST=OTESTD*RAD
      IF (NPRINT(6).NE.1) GO TO 7000
      WRITE (6,8007)
      WRITE (6,7011) FK
7011     FORMAT ('0The momentum is given in terms of K = ',1PE23.16,
     1 ' fm-1')
      IF(NCM.LE.0)WRITE(8,7719)
7719     FORMAT ('0The upper entry is the cross section in the lab ',
     1 'system corrected for effects due to the electron mass'/
     2 ' The lower entry is the cross section in the C/M system.')
      WRITE(8,7724)
7724     FORMAT ('0All cross sections in fm per steradian')
      WRITE (6,8243)
8243     FORMAT ('0Cross sections Ravenhall program'/
     1 '0',T9,'Angle in degrees',T40,'Cross section',T68,
     2 'Real amplitude',T92,'Imaginary amplitude')
7000      DO 7001 NAR=1,NANG
      NFLDL=2
      NFLDH=2
      DTH=0.D0
      IF (NCM.NE.0.OR.(DELTH(NAR).EQ.0.D0.AND.
     1 DELPH(NAR).EQ.0.D0.AND.XMS(NAR).EQ.0.D0)) GO TO 7020
      THMS=XMS(NAR)
      IF (THMS.GE.0.D0) GO TO 7019
C    Determine the multiple scattering angle from the target
C    thickness.
      Z=DABS(GAM/ALPHA)
      THMS=DABS(THMS)/DCOS(THELAB(NAR)*.5D0*DTR)
      THMS=RAD*DSQRT(.157D0*Z*(Z+1)*THMS/(ATWT*ELAB**2)*
     1 (-.047D0+2.583D0*DLOG10(7800.D0*(Z+1)*Z**.333333333D0*
     2 THMS/((1.D0+(XME/ELAB)**2)*ATWT*(1.D0+3.35D0*DABS(GAM)**2)))))
      IF (NPRINT(10).NE.0) WRITE (6,7702) THMS
7702     FORMAT (' RMS multiple scattering angle calculated to be ',
     1 1PE23.16,' degrees.')
C    Determine the step size for folding.  This is taken to be the smallest
C    of all of the non-0 folding angles.  If the step size is 0, then no
C    folding is done.
7019      DTH=1.D37
      IF (DELTH(NAR).GT.0.D0) DTH=DELTH(NAR)
      IF (DELPH(NAR).GT.0.D0) DTH=DMIN1(DTH,DELPH(NAR))
      IF (THMS.GT.0.D0) DTH=DMIN1(DTH,THMS)
      DTH=DTH*DTR
      NFLDL=1
      NFLDH=3
7020      DO 4365 NFLD=NFLDL,NFLDH
      IF (NCM.EQ.0) GO TO 7021
      O=THECM(NAR)*DTR
      GO TO 7717
7021      OL=THELAB(NAR)*DTR+DFLOAT(NFLD-2)*DTH
C
C    Change to the C/M angles
      O=DATAN2(DSIN(OL)*DSQRT(1.D0+2.D0*ELAB/AMEV),
     1 (1.D0+ELAB/AMEV)*DCOS(OL)-ELAB/AMEV)
C
C    Legendre function recurrence relation
7717      P(1)=1.0D0
      P(2)=DCOS(O)
      LF=NFL+1
      DO 660 L=2,NFL
      FL=L
660      P(L+1)=((2.0D0*(FL-1.0D0)+1.0D0)*P(2)*P(L)-(FL-1.0D0)*P(L-1))/FL
      IF(NPRINT(8).NE.0)WRITE(8,690)(P(L),L=1,LF)
690     FORMAT ('0Legendre function'/(4(1X,1PE23.16)))
C
C    Test on size of angle.  If less than OTEST, calculate constant
C    factor for F(O)
      IF (O-OTEST) 760,902,902
760      A1=(DSIN(O/2.0D0))**2
      A2=DLOG(A1)
      IF(NPRINT(8).NE.0)WRITE(8,5000)A2
5000     FORMAT (' LOG = ',1PE23.16)
      A3=GAM*A2
      A4=GAM/(2.0D0*FK*A1)
      TT=A4*BETA0*CDEXPI(A3)
      IF(NPRINT(8).NE.0)WRITE(8,860)BETA0,A1,A4,TT
860     FORMAT ('0Small O'/(4(1X,1PE23.16)))
      IF (NREC) 4071,4071,901
902      IF (NLAT) 4071,4071,1901
1901      NLAT=0
C
C    Now sum the Faxen-Holtzmark series.  First initialize for the
C    N loop.
901      NREC=0
      ETA=DCMPLX(0.0D0,0.0D0)
      ETN=DCMPLX(0.0D0,0.0D0)
      PSNM=DCMPLX(0.0D0,0.0D0)
      DO 3900 N=1,NFL
      BETA=DCMPLX(0.0D0,0.0D0)
      FN=N
      TAND=0.0D0
      IF(FN.LE.FNS)TAND=TND(N)
      RHO=DSQRT(FN*FN-GAM*GAM)
C
C    exp(2I ETA)
      AA=CDGAM(GAMC-(GAM*GAM)/(RHO+FN),N,NTEST,1)
      ETA=CDEXPI(PI*(FN-RHO))*(RHO+GAMM)*CDIVIDE(DCONJG(AA),FN*AA)
C
C    Calculate the coefficients in the series. Check the angle size.
      IF(O-OTEST)3480,3690,3690
C
C    Calculation of BETA(N).  The small angle calculation:
3480      BETA=CDIVIDE(CDGAM(GAMM,N,NTEST,1),CDGAM(GAMC,N,NTEST,1))
C
C    Calculation of A(L) (Legendre coefficients)
      ETN=FN*ETA
      TTT=BETA*(FN+FN-1.0D0)
C
C    Large angle coefficients for O larger than OTEST
C    Calculate the phase shift exponential.
3690      PS=ETA*DCMPLX(1.0D0-TAND**2,2.0D0*TAND)/(1.0D0+TAND**2)
C
C    Calculation of A(L) (Legendre coefficient)
      PSN=PS*FN
      BB=PSN+PSNM
      IF(O-OTEST)3852,3870,3870
3852      BB=BB-TTT
C
C    Coefficients for any angle have been found
3870      A(N+1)=BB
      IF(NPRINT(8).NE.0)WRITE(8,*)N,RHO,ETA,BETA,TTT,PS,PSN,BB
3900      PSNM=PSN
C
C    Perform the reduction on the coefficients
C    Recursion relations for reduction of the series for the
C    scattering amplitude.  (see Yennie, Ravenhall, and Wilson
C    Phys Rev 94)
      MF=NFL
C    NORD-500 Fortran has problems with some of the following:
      A(1)=DCMPLX(0.0D0,0.0D0)
      DO 3901 J=2,MF
      TEMP1=DFLOAT(J-1)/DFLOAT(J+J-1)
      TEMP2=DFLOAT(J-2)/DFLOAT(J+J-5)
3901      B(J)=A(J)-TEMP1*A(J+1)-TEMP2*A(J-1)
      IF(NPRINT(8).NE.0)WRITE(8,3960)(A(J),B(J),J=1,MF)
3960     FORMAT ('0(1-cos(O))F(O) = '/(4(1X,1PE23.16)))
      MF=NFL-1
      B(1)=DCMPLX(0.0D0,0.0D0)
      DO 3902 J=2,MF
      TEMP1=DFLOAT(J-1)/DFLOAT(J+J-1)
      TEMP2=DFLOAT(J-2)/DFLOAT(J+J-5)
3902      A(J)=B(J)-TEMP1*B(J+1)-TEMP2*B(J-1)
      IF(NPRINT(8).NE.0)WRITE(8,4010)(B(J),A(J),J=1,MF)
4010     FORMAT ('0(1-cos(O))2F(O) = '/(4(1X,1PE23.16)))
      MF=NFL-2
      A(1)=DCMPLX(0.0D0,0.0D0)
      DO 3903 J=2,MF
      TEMP1=DFLOAT(J-1)/DFLOAT(J+J-1)
      TEMP2=DFLOAT(J-2)/DFLOAT(J+J-5)
3903      B(J)=A(J)-TEMP1*A(J+1)-TEMP2*A(J-1)
      IF(NPRINT(8).NE.0)WRITE(8,4060)(A(J),B(J),J=1,MF)
4060     FORMAT ('0(1-cos(O))3F(O) = '/(4(1X,1PE23.16)))
C
C    Calculate the cross section
4071      SUM=DCMPLX(0.0D0,0.0D0)
C
C    Sum the series
      JF=MF-1
      DO 4160 J=1,JF
4160      SUM=SUM+B(J+1)*P(J)
C
C    (1-cos(O))**3F(O)
      A1=(1.0D0-P(2))**3
      CFO=SUM/DCMPLX(0.0D0,FK+FK)
      FO=CFO/A1
      IF(O-OTEST)4210,4240,4240
C
C    Take out the non-relativistic point amplitude for the small
C    angle calculation
4210      FO=FO+TT
4240      IF(NPRINT(8).NE.0)WRITE(8,4290)CFO,FO,A1,SUM
4290     FORMAT ('0Cross section'/4(1X,1PE23.16)/3(1X,E23.16))
C
C    Scattering amplitude found.
C          FO      Complex scattering amplitude
C          DCS     Differential cross section in C/M
      DCS=DCABS(FO)**2*(1.0D0+DTAN(O/2.0D0)**2)
      IF(NPRINT(8).NE.0)WRITE(8,4350)DCS
4350     FORMAT ('0Differential cross section = ',1PE23.16)
      ODEG=O*RAD
      IF (NCM) 7721,7721,7723
C
C    Change the cross section to the lab system
C          DCC     Cross section in C/M system
C          DCS     Cross section in lab system corrected for effects
C            due to electron mass
7721      DCC=DCS
      DCS=((DSQRT(PCM**2+XME**2)+ECMN)/(AMEV+ELAB*(1.0D0-DCOS(OL))))**2*
     1 (1.0D0+(XME/ELAB)**2)*(1.0D0+(XME* DTAN(OL/2.0D0)/ELAB)**2)*DCS
      IF(NPRINT(6).EQ.0)GO TO 7751
      IF (LEAF.LT.16) GO TO 1500
      WRITE (6,8007)
      WRITE (6,7011) FK
      WRITE (6,7719)
      WRITE (6,7724)
      WRITE (6,8243)
      LEAF=0
1500      OLDEG=OL*RAD
      WRITE (6,7722) OLDEG,DCS
7722     FORMAT ('0',4(1X,1PE23.16))
      LEAF=LEAF+1
      WRITE (6,4351) ODEG,DCC,FO
4351     FORMAT (' ',4(1X,1PE23.16))
      GO TO 7751
7723      IF(NPRINT(6).EQ.0)GO TO 7751
      IF (LEAF.LT.16) GO TO 1600
      WRITE (6,8007)
      WRITE (6,7011) FK
      WRITE (6,7724)
      WRITE (6,8243)
      LEAF=0
1600      WRITE (6,4351) ODEG,DCS,FO
      LEAF=LEAF+1
C
C    Fold the cross section
7751      IF (NCM.EQ.0) GO TO 7752
      SIGCM(NAR)=DCS
      IF (NPUNCH(3).NE.0) WRITE (9,4351) SIGCM(NAR),THECM(NAR)
      GO TO 7001
7752      IF (NFLD.NE.2) GO TO 4365
      OLS=OL
      SIGCM(NAR)=DCC
      THECM(NAR)=ODEG
      SIGLAB(NAR)=DCS
4365      S(NFLD)=DCS
      IF (DTH.EQ.0.D0) GO TO 7753
      CTHMS=(DTR*THMS)**2/4.D0
      DSDT=DLOG(S(3)/S(1))/(2.D0*DTH)
      D2SDT2=DLOG(S(3)*S(1)/SIGLAB(NAR)**2)/DTH**2+DSDT**2
      SIGLAB(NAR)=SIGLAB(NAR)*(1.D0 + DSDT*DCOS(OLS)*((DTR*DELPH(NAR))
     1 **2/6.D0 + CTHMS)/DSIN(OLS) + D2SDT2*((DTR*DELTH(NAR))
     2 **2/6.D0 + CTHMS))
      IF (NPRINT(10).NE.0) WRITE (6,7754) S(2),SIGLAB(NAR)
7754     FORMAT (' Unfolded cross section = ',1PE23.16/
     1 ' Folded cross section = ',E23.16)
7753      IF (NPUNCH(3).NE.0) WRITE (9,4351) SIGLAB(NAR),THELAB(NAR)
7001      CONTINUE
      RETURN
      END

      SUBROUTINE CFDEQF
C    Coulomb functions differential equations
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ZV10(4), ZV11(4), ZV12(3), ZV13(4)
      DIMENSION ZV14(11), ZV15(10)
      EQUIVALENCE (X01,ZV10(1)),(X02,ZV10(2)),(DX01,ZV11(1))
      EQUIVALENCE (DX02,ZV11(2)),(XX,ZV12(1))
      COMMON /DEQ/ ZV10,ZV11,ZV12,ZV13,ZV14,ZV15,
     1 GAM,FN,FFI,FFJ,FUN,CD,XDEQ,MDEQ,N,NST,FK
C
C    Differential equations
      VFUNCT=-GAM/XX
      DX01=(FN/XX)*X01-(1.0D0-VFUNCT)*X02
      DX02=-(FN/XX)*X02+(1.0D0-VFUNCT)*X01
      RETURN
      END

      SUBROUTINE RAVDEQ
C    Differential equations for the Ravenhall program
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Include parameter file to define array sizes
      INCLUDE 'ravpro.inc'
      DIMENSION ZV10(4), ZV11(4), ZV12(3), ZV13(4)
      DIMENSION ZV14(11), ZV15(10)
      EQUIVALENCE (X01,ZV10(1)),(X02,ZV10(2)),(X03,ZV10(3))
      EQUIVALENCE (X04,ZV10(4)),(DX01,ZV11(1)),(DX02,ZV11(2))
      EQUIVALENCE (DX03,ZV11(3)),(DX04,ZV11(4)),(XX,ZV12(1))
      COMMON /DEQ/ ZV10,ZV11,ZV12,ZV13,ZV14,ZV15,
     1 GAM,FN,FFI,FFJ,FUN,CD,XDEQ,MDEQ,N,NST,FK
      COMMON /POT/ POT(NPOTMX)
C
      IF (MDEQ) 20,10,20
10      CALL CHGDE2 (FK,XX,FUN)
      DX01=FUN*(XX)**2
      DX02=FUN*XX
      DX03=DX01*(XX**2)
C
C    Normalization integrals
      RETURN
20      I=INT((XX*2.0D0/XDEQ)+.51D0)
      IF (I.GT.NPOTMX) THEN
          WRITE (6,25) I,NPOTMX
25     FORMAT (' The array dimension for the POT array is ',
     1    'too small.'/
     2    ' An array dimension of at least',I8,' is required ',
     3    'whereas POT is dimensioned at only',I8,'.')
          STOP 'RAVDEQ aborted -- POT array index too big'
      ENDIF
      IF (N-NST) 30,30,40
C
C    First phase shift, obtain the potential.
30      CALL CHGDE2 (FK,XX,FUN)
      DX01=FUN*(XX**2)
      DX02=FUN*XX
      CD=(-GAM/FFI)*((X01/XX)+FFJ-X02)
      POT(I)=CD
      DX03=FN/XX*X03-(1.0D0-CD)*X04
      DX04=-FN/XX*X04+(1.0D0-CD)*X03
      RETURN
40      CD=POT(I)
      DX01=FN/XX*X01-(1.0D0-CD)*X02
      DX02=-FN/XX*X02+(1.0D0-CD)*X01
      RETURN
      END

      SUBROUTINE COULS (GAMP,X0,NF)
C    Coulomb function subroutine
C
C    Standard input data:
C          GAM     Atomic number times the fine structure constant
C            This must be identical to the GAMMA read into the
C            Ravenhall program.  This number is usually given to
C            six places.
C          X0      The fitting-on radius
C          NF      Number of Coulomb functions calculated
C          NPRINT(7) Set to non-zero value to write intermediate output.
C
C    Default options, and how to change them:
C          FCT1    (default=0.75) This is the factor multiplying N to
C            obtain the value of X, the independant variable at
C            which the series solution is calculated.  In order
C            to change this number, one may read it in by removing
C            the C* below.  A recompilation is necessary.  A card
C            giving the desired value of FCT1 must then be added
C            to the data set.
C          ZV15(1) (default=0.01) the step size
C          NZ002(1)(default=50) print frequency.  If either of these
C            values must be changed, one may read in the new
C            values by removing the C] below.  A recompilation
C            is again necessary.  New values must then be
C            included in the data set.  A maximum of 10
C            different regions may be specified
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Include parameter file to define array sizes
      INCLUDE 'ravpro.inc'
      DIMENSION FFRN(NFMX), GGRN(NFMX), FFIN(NFMX), GGIN(NFMX)
      DIMENSION NPRINT(10), NPUNCH(5)
      DIMENSION ZV10(4), ZV11(4), ZV12(3), ZV13(4)
      DIMENSION ZV14(11), ZV15(10)
      DIMENSION NZ001(10), NZ002(11)
      DOUBLE COMPLEX AA,BB,CC
      DOUBLE COMPLEX CDGAM
      EXTERNAL CFDEQF
      EQUIVALENCE (X01,ZV10(1)),(X02,ZV10(2)),(DX01,ZV11(1))
      EQUIVALENCE (DX02,ZV11(2)),(XX,ZV12(1))
      COMMON /DEQ/ ZV10,ZV11,ZV12,ZV13,ZV14,ZV15,
     1 GAM,FN,FFI,FFJ,FUN,CD,XDEQ,MDEQ,N,NST,FK
      COMMON /COUL/ FFRN,GGRN,FFIN,GGIN
      COMMON /IO/ NPRINT,NPUNCH
      COMMON /CONST/ XMN,HBARC,ALPHA,XME,PI,DTR,RTD
      XXP1(X)=X*X+X
C
C    Set up intial default conditions
      GAM=GAMP
      FCT1=0.75D0
      ZVSAV=ZV15(1)
      ZV15(1)=0.009999999999946000D0
      ZV15(2)=0.0D0
      NZ002(1)=50
      ZV14(1)=0.0D0
      ZV14(2)=X0
C
C    Read in FCT1
C*      READ (7,400) FCT1
C*400     FORMAT(D10.5)
C
C    Read in integration variables
C          IZ2     Number of regions
C          ZV14(i) Starting value of region i
C          ZV15(i) Step size in region i
C    NZ002(i) Print frequency in region i
C]      READ (7,410) IZ2,(ZV14(I),ZV15(I),NZ002(I),I=1,IZ2)
C]410     FORMAT(I5/(2D10.5,I5))
C]      J=IZ2-1
C]      DO 420 I=1,J
C]      IF(DABS(ZV15(I))-1.0D-8)430,430,440
C]430      WRITE(8,435)
C]435     FORMAT('0Zero step size error')
C]      CALL EXIT
C]440      CONTINUE
C
C    Write out the input
      IF(NPRINT(1).EQ.0)GO TO 8015
      WRITE (6,8014) GAM,FCT1,ZV14(2),ZV15(1),NF,X0
8014     FORMAT ('1',T11,'- * - R A V E N H A L L  P R O G R A M  ',
     1 'C O U L O M B  F U N C T I O N S - * -'/
     2 '0Input for calculation of the Coulomb functions for ',
     3 'phase shift analysis'/
     4 ' GAMMA = ',1PE23.16/
     5 ' Integrate Coulomb equations from ',E23.16,' times N to r =',
     6 E23.16,' in steps of ',E23.16/
     7 1X,I3,' Phase shifts will be possible'/
     8 ' Fitting on radius = ',E23.16)
8015      DO 3240 N=1,NF
      NCD=1
C
C    Calculation of RHO
      FN=N
      RHO=DSQRT(FN*FN-GAM*GAM)
      FCT=FCT1*FN
C
C    Series solution (first do the regular function, then do the
C    irregular function).
      IF (FCT-X0) 1004,1009,1009
1004      XSTO=X0
      X0=FCT
      GO TO 1010
1009      XSTO=X0
1010      T=RHO
      M=2
1020      AM1=1.0D0
      BM1=GAM/(T+FN)
      G=0.0D0
      F=0.0D0
      FPRE=F
      GPRE=G
      MC=1
1070      M1=MC-1
      IF (M1) 1072,1072,1080
1072      G=1.0D0
      F=BM1
      GO TO 1100
1080      GPRE=G
      FPRE=F
      G=G+AM1
      F=BM1+F
      IF (DABS(GPRE)-DABS(G)) 1081,1081,1082
1081      GSAV=G
1082      IF (DABS(FPRE)-DABS(F)) 1083,1083,1084
1083      FSAV=F
1084      CONTINUE
1100      FM=MC
      AM=AM1
      BM=BM1
      IF (T) 1160,1130,1130
1130      AM1=X0*((-GAM*AM-(T+FM+FN)*BM)/FM)/(FM+2.0D0*T)
      BM1=X0*((-GAM*BM+(T+FM-FN)*AM)/FM)/(FM+2.0D0*T)
      GO TO 1190
1160      N2=2*N
      IF (MC-N2) 1130,1180,1130
1180      AM1=X0*(-GAM*AM-(T+FM+FN)*BM)*((FN+RHO)/(4.0D0*FN*GAM**2))
      BM1=X0*(-GAM*BM+(T+FM-FN)*AM)*((FN+RHO)/(4.0D0*FN*GAM**2))
1190      T1=DABS(G)
      T2=DABS(F)
      T11=DABS(AM1)
      T22=DABS(BM1)
      T111=T11/T1
      T222=T22/T2
      IF (T111-1.0D-12) 1205,1205,1210
1205      IF (T222-1.0D-12) 1230,1230,1210
1210      MC=MC+1
      GO TO 1070
1230      IF (T) 1280,1280,1240
1240      FR=F
      GR=G
      IF(NPRINT(7).NE.0)WRITE(8,1085)GSAV,FSAV
1085     FORMAT ('0Max term GR = ',1PE23.16,', Max term FR = ',E23.16)
      T=-T
      GO TO 1020
1280      FI=F
      GI=G
      T=-T
      IF(NPRINT(7).EQ.0)GO TO 1340
      WRITE(8,1085)GSAV,FSAV
      WRITE (6,1320) X0,N,RHO,GAM,FR,FI,GR,GI,T
1320     FORMAT ('0X0 = ',1PE23.16,', N = ',I3,', RHO = ',E23.16,
     1 ', GAMMA = ',E23.16/
     2 ' FR = ',E23.16,', FI = ',E23.16,', GR = ',E23.16,', GI = ',
     3 E23.16/
     4 ' T = ',E23.16)
C
C    Calculate GAMMA functions
C    Normalized Coulomb function
1340      TT=-(GAM*GAM)/(RHO+FN)
      RA1=TT
      RA2=TT+TT
      AA=DCMPLX(RA1,GAM)
      BB=DCMPLX(RA2,0.0D0)
      CC=CDGAM(AA,0,NTEST,0)
      S1=1.0D0/DCABS(CC)
      S2=1.0D0/DBLE(CDGAM(BB,0,NTEST,0))
      TERM=2.0D0*X0*DCABS(AA)/XXP1(RA2)
      IF(N.LT.2)GO TO 1542
      NTF=N-1
      DO 1541 I=1,NTF
1541      TERM=TERM*2.0D0*X0*DCABS(AA+I)/XXP1(RA2+I+I)
1542      SERS=S2*TERM/(2.0D0*RHO*S1)
      FNR=DEXP(0.5D0*PI*GAM)*DSQRT(FN*(FN+RHO)/2.0D0)
      A1=(2.0D0*X0)**TT
      FNR=FNR*A1
      FNR=FNR*SERS
      IF(NPRINT(7).NE.0)WRITE(8,1573)CC,FNR,SERS,S1,S2,TERM
1573     FORMAT ('0',4(1X,1PE23.16)/1X,3(1X,E23.16))
C
C    Normalized functions
      T=-T
      A1=PI*GAM
      A2=DSIN(PI*(FN-RHO))
      CHCK=-DABS(GAM)*DEXP(A1)*A2/((RHO+RHO)*DCOSH(A1)*DSQRT(A2*A2/
     1 (1.0D0-A2*A2)+DTANH(A1)**2))
      FNI=CHCK/FNR
      FRN=FR*FNR
      FIN=FI*FNI
      GRN=GR*FNR
      GIN=GI*FNI
      IF(NPRINT(7).NE.0)WRITE(8,1822)FRN,FIN,GRN,GIN
1822     FORMAT ('0Normalized F and G'/' FRN = ',1PE23.16,', FIR = ',
     1 E23.16,', GRN = ',E23.16,', GIN = ',E23.16)
      X0=XSTO
      IF (FCT-X0) 1862,3030,3030
C
C    Initial conditions for integration
1862      X01=GRN
      X02=FRN
      ZV14(1)=FCT
      NZ001(6)=002
      M=3
      NCD=1
      GO TO 9030
9015      IF(NCD-2)9040,1980,3030
1980      CGI=X01
      CFI=X02
      M=2
      NCD=3
      FRN=CFR
      FIN=CFI
      GRN=CGR
      GIN=CGI
      GO TO 3030
9030      CONTINUE
C
C    Integrate
      XX=ZV14(1)
      ZV12(2)=ZV15(1)
      ZV12(3)=ZV14(2)
      NZ001(9)=NZ002(1)
C
C    Call Runga-Kutta-Gill Routine
C          ZV12    Independent variable array
C          ZV10    Array of solutions (see equivalence statements)
C          ZV11    Array of derivatives (see equivalence statements)
C          ZV13    Array used by DEQ
C          NZ001(6)Number of equations
C          IFUN    Control return
C          X01     Solutions
C          X02     Solutions
C          DX01    Derivatives
C          DX02    Derivatives
C          XX      Independent variable (not suitabely named)
      IF (ZV12(1).LT.ZV12(3)) GO TO 5100
      WRITE (6,5101) ZV12(1),ZV12(3)
5101     FORMAT (' Bad DEQ call in COULS: ',2(1X,1PE23.16))
      STOP 'COULS aborted -- bad DEQ call'
5100      CALL INDEQ2(NZ001(9))
      NZ005=1
9011      CALL DEQ2(IRET,ZV12,ZV10,ZV11,ZV13,NZ001(6),CFDEQF)
      IF (IRET.EQ.0) GO TO 9013
      IF (NPRINT(7).EQ.0) GO TO 9011
      WRITE(8,9016)XX,DX01,DX02,X01,X02
9016     FORMAT ('0XX = ',1PE23.16,', GDOT = ',E23.16,', FDOT = ',E23.16/
     1 ' G = ',E23.16,', F = ',E23.16)
      GO TO 9011
9013      NZ005=NZ005+1
      ZV12(2)=ZV15(NZ005)
      ZV12(3)=ZV14(NZ005+1)
      NZ001(9)=NZ002(NZ005)
      IF (ZV15(NZ005)) 9011,9035,9011
9035      IF (FCT-X0) 9015,3030,3030
9040      CGR=X01
      CFR=X02
C
C    Integration finished for this case.  Now do the other function.
      X01=GIN
      X02=FIN
      ZV14(1)=FCT
      NCD=2
      GO TO 9030
C
C    Write out the normalized Coulomb functions.
C    The resulting Coulomb functions are in FRN, FIN, GRN, and GIN.
3030      IF (NPRINT(1).NE.0) WRITE (6,2117) FRN,GRN,FIN,GIN,FCT,X0
2117     FORMAT ('0FRN = ',1PE23.16,', GRN = ',E23.16/
     1 ' FIN = ',E23.16,', GIN = ',E23.16/
     2 ' FCT = ',E23.16,', X0 = ',E23.16)
C
C    Punch Coulomb functions.
      IF (NPUNCH(1).NE.0) WRITE (9,2217) FRN,GRN,N
      IF (NPUNCH(1).NE.0) WRITE (9,2217) FIN,GIN,N
2217     FORMAT (2E23.16,23X,I8)
      FFRN(N)=FRN
      GGRN(N)=GRN
      FFIN(N)=FIN
      GGIN(N)=GIN
C
C    Now do the next value of angular momentum
3240      CONTINUE
C
C    All N values done. Go to the next case.
      ZV15(1)=ZVSAV
      RETURN
      END

      FUNCTION CDGAM(AA,N,NTEST,NFACT)
C    Calculation of a GAMMA function of a complex argument.
C    Method taken from Davis Tables of Higher Mathematical Functions.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX AA,RD,RN,CDGAM
      DOUBLE PRECISION B(25)
      DATA B/-.4227843350984671D0,-.2330937364217867D0,
     1 .1910911013876915D0,-.0245524900054D0,-.0176452445501443D0,
     2 .0080232730222673D0,-.0008043297756043D0,-.0003608378162548D0,
     3 .0001455961421399D0,-.0000175458597517D0,-.258899502902D-5,
     4 .13385015469D-5,-.2054743149D-6,-.1595268D-9,.62756218D-8,
     5 -.12736143D-8,.923397D-10,.120030D-10,-.42207D-11,.5240D-12,
     6 -.139D-13,-.67D-14,.128D-14,.12D-15,-.2D-16/
      NN=N
      RN=DCMPLX(0.0D0,0.0D0)
      DO 300 I=1,22
300      RN=(RN+B(23-I))*AA
      RN=RN+1.0D0
      NFACT=NFACT
      NTEST=NTEST
C
C    Calculation of the coefficient for positive RHO
      IF(NN)350,350,390
350      CDGAM=1.0D0/(AA*(AA+1.0D0)*RN)
      RETURN
390      IF(NN-2)400,510,430
400      CDGAM=1.0D0/((AA+1.0D0)*RN)
      RETURN
430      IF(NN-3)450,435,450
435      CDGAM=(AA+2.0D0)/RN
      RETURN
450      RD=AA+2.0D0
      DO 480 J=4,NN
480      RD=RD*(AA+(J-1))
      CDGAM=RD/RN
      RETURN
510      CDGAM=1.0D0/RN
      RETURN
      END

      FUNCTION ERF(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(27),AA(18)
      DATA A(27)   / 3.887303655222904D0      /
      DATA A(26)   /-1.381631420019799D0      /
      DATA A(25)    / .647316404854584D0       /
      DATA A(24)    /-.305931024422036D0       /
      DATA A(23)    /.138679747202030D0        /
      DATA A(22)    /-.059247456591259D0       /
      DATA A(21)    / .236917518249282D-01   /
      DATA A(20)    /-.884736263524045D-02   /
      DATA A(19)    / .308566171136092D-02   /
      DATA A(18)    /-.100638635123798D-02   /
      DATA A(17)   /.307546328843079D-03    /
      DATA A(16)   /-.882619837553631D-04   /
      DATA A(15)   / .238450961660726D-04   /
      DATA A(14)   /-.607910028505827D-05   /
      DATA A(13)   / .146597217338083D-05   /
      DATA A(12)   /-.033515993427206D-05   /
      DATA A(11)   / .007280579544232D-05   /
      DATA A(10)   /-.001505791176668D-05   /
      DATA A(9)   / .000297094742055D-05   /
      DATA A(8)   /-.000056021273938D-05   /
      DATA A(7)   / .000010113162390D-05   /
      DATA A(6)   /-.1750650485D-10        /
      DATA A(5)   /.0291038139D-10        /
      DATA A(4)   /-.0046532645D-10        /
      DATA A(3)   / .0007164815D-10        /
      DATA A(2)   /-.0001063749D-10        /
      DATA A(1)   / .0000152467D-10        /
      DATA AA(18)  / 1.97070527225754D0       /
      DATA AA(17)   /-.143397402717750D-01   /
      DATA AA(16)   / .297361692202619D-03   /
      DATA AA(15)   /-.980351604336237D-05   /
      DATA AA(14)   / .043313342034728D-05   /
      DATA AA(13)   /-.002362150026241D-05   /
      DATA AA(12)   / .000151549676581D-05   /
      DATA AA(11)   /-.000011084939856D-05   /
      DATA AA(10)   / .0904259014D-10        /
      DATA AA(9)   /-.0080947054D-10        /
      DATA AA(8)  / .0007853856D-10        /
      DATA AA(7)  /-.0000817918D-10        /
      DATA AA(6)  / .90715D-15             /
      DATA AA(5)  /-.10646D-15             /
      DATA AA(4)  / .01315D-15             /
      DATA AA(3)  /-.00170D-15             /
      DATA AA(2)  / .00023D-15             /
      DATA AA(1)  /-.00003D-15             /
      DATA C /.282094792D0/
      S=DABS(X)
      IF(S.GT.4.D0) GO TO 20
      COEFF=.25D0*X*X-2.D0
C
C    Use recurrence relation to compute the Chebyshev polynomials for
C    the argument X/4
      BJ=0.D0
      BJ1=0.D0
      DO 10 J=1,27
      BJ2=BJ1
      BJ1=BJ
10      BJ=COEFF*BJ1-BJ2+A(J)
      ERF=.125D0*X*(BJ-BJ2)
      RETURN
C
20      Y=1.D0/S
      COEFF=64.D0*Y*Y-2.D0
      BJ=0.D0
      BJ1=0.D0
      DO 30 J=1,18
      BJ2=BJ1
      BJ1=BJ
30      BJ=COEFF*BJ1-BJ2+AA(J)
      ERF=DSIGN(1.D0-(BJ-BJ2)*C*Y*DEXP(-S*S),X)
      RETURN
      END

      FUNCTION CDEXPI(X)
C    Complex exponential of an imaginary argument.  Input argument,
C    X, is double precision.  Output argumen is the double complex
C    result (cos(x), sin(x))
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX CDEXPI
      CDEXPI=DCMPLX(DCOS(X),DSIN(X))
      RETURN
      END

      FUNCTION CDIVIDE(X,Y)
C    This routine is meant to divide two double complex numbers
C    without the problems inherent when both numbers are about
C    the square root of the largest real in magnetude.  This
C    will not be necessary on machines with large possible
C    exponents, but if the largest possible exponent in the
C    floating point format is 77 and smaller (as on the NORD-500)
C    then this routine is required.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX CDIVIDE,X,Y
      CDIVIDE=X
      XSIZ=DCABS(X)
      IF (XSIZ.EQ.0.0D0) RETURN
      XREL=DBLE(X)/XSIZ
      XIMA=DIMAG(X)/XSIZ
      YSIZ=DCABS(Y)
      IF (YSIZ.EQ.0.0D0) STOP 'CDIVIDE aborted -- division by zero'
      YREL=DBLE(Y)/YSIZ
      YIMA=DIMAG(Y)/YSIZ
      ZSIZ=XSIZ/YSIZ
      ZREL=ZSIZ*(XREL*YREL+XIMA*YIMA)
      ZIMA=ZSIZ*(XIMA*YREL-XREL*YIMA)
      CDIVIDE=DCMPLX(ZREL,ZIMA)
      RETURN
      END

      FUNCTION DCABS(X)
C    This routine returns the magnetude of a complex number even
C    under extremes of the exponent value of the magnetude.
C    This is because the following algorythm is used:
C    MAG=MAX(r,i)*SQRT((r/MAX(r,i))**2+(i/MAX(r,i))**2)
C          If MAX(r,i) is greater than 1
C    MAG=MIN(r,i)*SQRT((r/MIN(r,i))**2+(i/MIN(r,i))**2)
C          If MAX(r,i) is less than 1
C    where r and i are absolute values of the real and imaginary
C    parts and MAG=SQRT(r**2+I**2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX X
      REL=DABS(DBLE(X))
      DIM=DABS(DIMAG(X))
      D=DMAX1(REL,DIM)
      DCABS=0.0D0
      IF (D.EQ.0.0D0) RETURN
      IF (D.GT.1.0D0) GO TO 1
      C=DMIN1(REL,DIM)
      IF (C.NE.0.0D0) D=C
1      DCABS=D*DSQRT((REL/D)**2+(DIM/D)**2)
      RETURN
      END

      BLOCK DATA RAVBD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C    Include parameter file to define array sizes
      INCLUDE 'ravpro.inc'
      COMMON /CONST/ XMN,HBARC,ALPHA,XME,PI,DTR,RAD
C    Constants from RMP 52, No. 2 (April 1980) as reprinted in the
C    CERN Particle Properties Data Booklet.
C    One atomic mass unit in MeV
      DATA XMN/931.5016D0/
C    H-bar times c in MeV-fm
      DATA HBARC/197.32858D0/
C    The fine structure constant
      DATA ALPHA/7.29735D-3/
C    The mass of the electron in MeV
      DATA XME/.5110034D0/
C    The mathematical constant PI
      DATA PI/3.1415926535897D0/
C    PI/180
      DATA DTR/.0174533D0/
C    180/PI
      DATA RAD/57.29578D0/
C
      COMMON /POT/ POT(NPOTMX)
      DATA POT/NPOTMX*0.0D0/
      END
