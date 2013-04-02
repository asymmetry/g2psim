      SUBROUTINE DEQ2(IRET,X,Y,F,Q,NEQ,DIRIV)
C    This routine uses the Runge-Kutta-Gill method (Gill, Proc.
C    Campridge Phil. Soc., 47 (1950)) to solve a set of coupled
C    differential equations of the form:
C          d(Y1(x))/dx=F1(Y1(x),Y2(x),....Yi(x))
C          d(Y2(x))/dx=F2(Y1(x),Y2(x),....Yi(x))
C        .
C        .
C        .
C          d(Yi(x))/dx=Fi(Y1(x),Y2(x),....Yi(x))
C    It must be initialized by a call to INDEQ2.  IRET is set to
C    0 on return if the final value of the independent variable
C    is reached.  IRET is 1 for an intermediate print-out return.
C    X is an array of dimension 3 which contains the initial value,
C    the step, and the final value of the independent variable.
C    X(1) is destroyed.  The other 2 remain intact.  F is an array
C    of at least NEQ elements.  Each element is the value of the
C    right hand side of the corresponding equation at x=X(1).
C    The F array must be set by the subroutine DIRIV. Y is an
C    array of dimension at least NEQ in which DEQ2 returns the
C    solution at X(1) for IRET=1 or X(3) for IRET=0. Y must be
C    initialized in the calling program to its value at X(1).
C    Q is a scratch array of dimension at least NEQ.  DIRIV is a
C    subroutine which evaluates F given X(1) and the elements of
C    the array Y.  It must be declared EXTERNAL in the calling program.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL DIRIV
      DOUBLE PRECISION X(3),Y(NEQ),F(NEQ),Q(NEQ)
      COMMON /DEQ2C/ I,NSTEP
      DATA SQT2/1.414213562373095D0/
C
      H=X(2)
10      IF (I.EQ.NSTEP) THEN
          I=0
          IRET=1
          RETURN
      ENDIF
      I=I+1
      CALL DIRIV
      DO 30 J=1,NEQ
          Q(J)=H*F(J)
          Y(J)=Y(J)+.5D0*Q(J)
30      CONTINUE
      X(1)=X(1)+.5D0*H
      CALL DIRIV
      DO 40 J=1,NEQ
          XK=H*F(J)
          Y(J)=Y(J)+(1.0D0-1.0D0/SQT2)*(XK-Q(J))
          Q(J)=(2.0D0-SQT2)*XK+(-2.0D0+3.0D0/SQT2)*Q(J)
40      CONTINUE
      CALL DIRIV
      DO 50 J=1,NEQ
          XK=H*F(J)
          Y(J)=Y(J)+(1.0D0+1.0D0/SQT2)*(XK-Q(J))
          Q(J)=(2.0D0+SQT2)*XK+(-2.0D0-3.0D0/SQT2)*Q(J)
50      CONTINUE
      X(1)=X(1)+.5D0*H
      CALL DIRIV
      DO 60 J=1,NEQ
          Y(J)=Y(J)+H*F(J)/6.0D0-Q(J)/3.0D0
60      CONTINUE
      IF (2.D0*DABS(X(1)-X(3)).GE.DABS(H)) GO TO 10
      IRET=0
      RETURN
      END
