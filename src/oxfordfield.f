C======================================================================
      SUBROUTINE sda_ptf_solid_nh3(X,B)
C----------------------------------------------------------------------
C
C     Purpose and Methods: Calculate magnetic field of Helmholtz coils
C
C     Library belongs: libsda.a
C
C     Called by: sda_minit
C
C     Author    U. Hartfiel     (1982)
C     Modified: M.Guckes    July 1988
C               B.Niczyporuk Sep.1998 (Real magnet geometry utilized + ...)
C
C--  COIL PARAMETERS
C--  ===============
C
C--  N       Number of coils
C--  H0      Field in the origin
C--  RJ      Current density (to calculate the magnetic field)
C--  RI      Inner coil radius
C--  RA      Outer coil radius
C
C--  ZAU     Axial distance of the gap to R-axis at RI
C--  ZAO     Axial distance of the gap to R-axis at RA
C--  ZBU     Axial distance of the outer edge to R-axis at RI
C--  ZBO     Axial distance of the outer edge to R-axis at RA
C
C          
C          A R
C          I
C          I               ZA0      ZBO 
C       RA +	             ________
C          I                | 	     |
C          I		        |	     | 
C          I        	    |	     |
C       RI +                |________|
C          I      ____     ZAU      ZBU 
C          I     |    |         
C          I     |____|
C          I      
C       ---0-----+----+-----+--------+---->  Z
C          I 
C
C
C----------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------
C
      SAVE
C
C Subroutine variables, By Jixie:  units: position in cm, field in kG
      REAL X(3), B(3)
C
C Local variables
C
C Number of loops
      INTEGER i, NLOOP
      PARAMETER (NLOOP = 14)
C
      REAL RI(NLOOP), ZAU(NLOOP), ZBU(NLOOP), RA(NLOOP),
     1     ZAO(NLOOP), ZBO(NLOOP), RJ(NLOOP), RJQ1(NLOOP), TURNS(NLOOP)
C
      REAL H0, S5, S6, D5, R, Z, CURm
C
      LOGICAL LFIRST
C
C Max field [KG]
      DATA H0 /50.0/
C
C Maximum Current [A]
c*    DATA CURm/123.646/
      DATA CURm/106.0958/
C
C Oxford coils geometry
      DATA RI /6.478,   6.478, 10.963,  10.963, 11.648,  11.648,
     .        14.451,  14.451, 18.936,  18.936, 19.559,  19.559,
     .        22.424,  22.424/
      DATA RA /8.487,   8.487, 11.648,  11.648, 13.877,  13.877,
     .        17.465,  17.465, 19.510,  19.510, 21.908,  21.908,
     .        25.430,  25.430/
C
      DATA ZAU/6.042,  -9.047,  6.967, -13.935,  6.967, -13.935,
     .         8.709, -14.930, 10.949, -16.920, 10.949, -16.920,
     .        12.690, -16.323/
      DATA ZBU/9.047,  -6.042, 13.935,  -6.967, 13.935,  -6.967,
     .        14.930,  -8.709, 16.920, -10.949, 16.920, -10.949,
     .        16.323, -12.690/
C
      DATA ZAO/6.042,  -9.047,  6.967, -13.935,  6.967, -13.935,
     .         8.709, -14.930, 10.949, -16.920, 10.949, -16.920,
     .        12.690, -16.323/
      DATA ZBO/9.047,  -6.042, 13.935,  -6.967, 13.935,  -6.967,
     .        14.930,  -8.709, 16.920, -10.949, 16.920, -10.949,
     .        16.323, -12.690/
C
C Number of turns per coil
      DATA TURNS/1118.0, 1118.0, 616.0, 616.0, 2003.0, 2003.0,
     .           3515.0, 3515.0, 638.0, 638.0, 3534.0, 3534.0,
     .           3853.0, 3853.0/
C                    
      DATA LFIRST /.TRUE./
C
      IF (LFIRST) THEN
         LFIRST = .FALSE.
C
C Calculate Current Density [A/CM**2] for each coil
         DO i = 1,NLOOP
           RJ(i) = TURNS(i)*CURm/( (RA(i) - RI(i))*(ZBU(i) - ZAU(i)) )
         ENDDO
C
C By Chao: special treat for the 1st pair of coils (inversed wired)
         RJ(1) = -1 * RJ(1)
         RJ(2) = -1 * RJ(2)
C
C  Integration over D(I) and all coils
         CALL DUSP(NLOOP,0.,0.,RJ,RI,ZAU,ZBU,RA,ZAO,ZBO,H0,S5,S6,D5)
         DO i =1,NLOOP
           RJQ1(i) = RJ(i)*(H0/S6)
         ENDDO
C   By Jixie: There is some problem in runing this IO in windows, 
c   so I just comment out the following 4 lines
C      WRITE(6,101) RJ,RJQ1
C 101     FORMAT(' RJ  =',8F9.2/' RJQ1=',8F9.2)
C         WRITE(6,102) H0,S5,S6,D5
C 102     FORMAT(' H0,S5,S6,D5 =',4F10.4)  
      ENDIF
C
C
      R = SQRT(X(1)**2+X(2)**2)
      Z = X(3)
C      print *,r,z
C
      CALL DUSP(NLOOP,R,Z,RJ,RI,ZAU,ZBU,RA,ZAO,ZBO,H0,S5,S6,D5)
C
      IF (R.EQ.0) R=1.E-5
      B(1) = X(1)*(S5/R)
      B(2) = X(2)*(S5/R)
      B(3) = S6
C
      RETURN
      END
