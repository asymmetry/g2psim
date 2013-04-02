      SUBROUTINE INDEQ2(ISTEP)
C    This routine initializes the routine DEQ2 which uses the
C    Runge-Kutta-Gill method (Gill, Proc. Campridge Phil. Soc.,
C    47 (1950)) to solve a set of coupled differential equations.
C    The argument ISTEP is the number of steps between print-out
C    returns (IRET=1) when DEQ2 is running. Set IRET negative or
C    0 if no print-out returns are desired.  Set it to 1 to do a
C    print-out return for every step.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /DEQ2C/ I,NSTEP
      IF (ISTEP.GT.0) THEN
          NSTEP=ISTEP
          I=ISTEP
      ELSE
          NSTEP=-1
          I=0
      ENDIF
      RETURN
      END
