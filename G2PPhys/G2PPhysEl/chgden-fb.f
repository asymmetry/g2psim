      SUBROUTINE CHGDE1(A1,A2,A3,A4)
C    This is a straight Fourier-Bessel series charge distribution.
C    The COEF's multiply j0 (not like Friar and Negele but like SESFIT)
C    Parameter 1 is a charge density normalization parameter.
C    Parameter 2 is the cutoff radius (which should not be fit).
C    Parameters 3, 4,... are coefficients 2, 3, ...  Coefficient 1 is
C    always 1.  The model number should be the number of coefficients
C    to use (including the first).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION COEF(48)
      COMMON /CHARGE/ MODELN,RHON,RMAX,COEF
      DATA PI/3.1415926535897D0/
      DATA RMIN/.1D-6/
C
      NCOEF=MODELN-1
100      A2=1.D0
      A3=0.D0
      A4=-1.D0
      DO 205 I=1,NCOEF
          A2=A2+COEF(I)
          A4=A4-COEF(I)*DFLOAT(I+1)**2
205      CONTINUE
      A2=A2*RHON/A1**3
      RHO0=A2
      A4=A4*RHON*(PI/RMAX)**2/A1**5
      RETURN
C
      ENTRY CHGDE2(A1,A2,A3)
      R=A2/A1
      A3=0.D0
      IF (R.GT.RMAX) RETURN
      IF (R.LT.RMIN) GO TO 3
      A3=DSIN(PI*R/RMAX)*RMAX/(PI*R)
      NCOEF=MODELN-1
      DO 1 I=1,NCOEF
          ARG=(I+1)*PI*R/RMAX
          A3=A3+COEF(I)*DSIN(ARG)/ARG
1      CONTINUE
      A3=A3*RHON/(A1**3)
      RETURN
3      A3=RHO0
      RETURN
C
      ENTRY CHGDE3(NPAR)
      NPAR=MODELN+1
      RETURN
      END
