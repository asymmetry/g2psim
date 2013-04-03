C     This file is included in RAVPRO in order to define
C     array sizes.
C
C     The size of the POT array.  This should be greater than 2 times
C     the maximum fitting-on radius divided by the minimum step size.
C     The value originally in RAVPRO was 12000.
      PARAMETER (NPOTMX=50000)
C
C     The maximum number of phase shifts (ie maximum value of NF).  This
C     can be estimated at about 3 times the maximum fitting-on radius,
C     or the maximum energy divided by 10.  The value originally in RAVPRO
C     was 60.
      PARAMETER (NFMX=60)
