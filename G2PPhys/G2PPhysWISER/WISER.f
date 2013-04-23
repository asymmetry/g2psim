********************************************************************8
*     Chao Gu 2013.04
*     Modified from X. Qian's modified version at 2005
********************************************************************

      SUBROUTINE wiser(Z,N,PART,E,P,TH,radlen,XSEC)
! Z, N are target info
! PART is the particle whose cross section is to be calcualted
! E is beam energy in GeV
! P is scattered particle momentum in GeV
! TH is scattering angle in rad
! radlen is total radiation length
! XSEC is the output cross section in nbarn/GeV/sr, per nucleon
! IMPORTANT: to get the cross section per nuclei, times this by (A)**0.8

      IMPLICIT double precision (A-H,O-Z)
      INTEGER Z,N,PART

      PI=3.14159265358979323846D0

      E1=E*1000.
      P1=P*1000.
      TH1=TH*180./PI
      IA=Z+N
      radlen1=radlen*100 ! in % now

      Call WISER_ALL_SIG(E1,P1,TH1,radlen1,PART,TOTAL)

      XSEC=TOTAL ! TOTAL in nanobarn/GeV*str
      XSEC=XSEC*1.0D-6 ! ub/MeV*str

      END 

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C Xiaochao Zheng, this subroutine added July 07, 2000

      double precision FUNCTION QUADMO(FUNCT,PLOWER,PUPPER,EPSLON,NLVL)
      double precision FUNCT,PLOWER,PUPPER,EPSLON
      INTEGER NLVL                                                  
      INTEGER LEVEL,MINLVL/3/,MAXLVL/24/,IRETRN(50),I                 
      double precision VALINT(50,2), VMX(50), RX(50), FMX(50), FRX(50),
     1   FMRX(50), ESTRX(50), EPSX(50)                                 
      double precision  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR,
     1   ESTINT, VL, AREA, ABAREA, VM, COEF, ROMBRG, EPS 

!      IMPLICIT double precision (A-H,O-Z)
         LEVEL = 0                                                     
         NLVL = 0                                                      
         ABAREA = 0.0                                                  
         VL = PLOWER                                                     
         R = PUPPER                                                     
         FL = FUNCT(VL)                                                 
         FM = FUNCT(0.5*(VL+R))                                         
         FR = FUNCT(R)                                                 
         EST = 0.0                                                     
         EPS = EPSLON                                                  
  100 LEVEL = LEVEL+1                                                  
      VM = 0.5*(VL+R)                                                    
      COEF = R-VL                                                       
      IF(COEF.NE.0) GO TO 150                                          
         ROMBRG = EST                                                  
         GO TO 300                                                     
  150 FML = FUNCT(0.5*(VL+VM))                                           
      FMR = FUNCT(0.5*(VM+R))                                           
      ESTL = (FL+4.0*FML+FM)*COEF                                      
      ESTR = (FM+4.0*FMR+FR)*COEF                                      
      ESTINT = ESTL+ESTR                                               
      AREA=ABS(ESTL)+ABS(ESTR)                                       
      ABAREA=AREA+ABAREA-ABS(EST)                                     
      IF(LEVEL.NE.MAXLVL) GO TO 200                                    
         NLVL = NLVL+1                                                 
         ROMBRG = ESTINT                                               
         GO TO 300                                                     
  200 IF((ABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.
     1         (LEVEL.LT.MINLVL))  GO TO 400                           
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                            
  300    LEVEL = LEVEL-1                                               
         I = IRETRN(LEVEL)                                              
         VALINT(LEVEL, I) = ROMBRG                                     
         GO TO (500, 600), I                                           
  400    IRETRN(LEVEL) = 1                                              
         VMX(LEVEL) = VM                                                 
         RX(LEVEL) = R                                                 
         FMX(LEVEL) = FM                                               
         FMRX(LEVEL) = FMR                                             
         FRX(LEVEL) = FR                                               
         ESTRX(LEVEL) = ESTR                                           
         EPSX(LEVEL) = EPS                                             
         EPS = EPS/1.4                                                 
         R = VM                                                         
         FR = FM                                                       
         FM = FML                                                      
         EST = ESTL                                                    
         GO TO 100                                                     
  500    IRETRN(LEVEL) = 2                                              
         VL = VMX(LEVEL)                                                 
         R = RX(LEVEL)                                                 
         FL = FMX(LEVEL)                                               
         FM = FMRX(LEVEL)                                              
         FR = FRX(LEVEL)                                               
         EST = ESTRX(LEVEL)                                            
         EPS = EPSX(LEVEL)                                             
         GO TO 100                                                     
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                         
      IF(LEVEL.GT.1) GO TO 300                                         
      QUADMO = ROMBRG /12.0D0                                          
      RETURN                                                           
      END                                                              

      Subroutine WISER_ALL_SIG(E0MM,PMM,THETA_DEG,RAD_LEN,TYPE,SIGMA)

!------------------------------------------------------------------------------
! Calculate pi,K,p  cross section for electron beam on a proton target
! IntegrateQs over function WISER_FIT using integration routine QUADMO
! E0         is electron beam energy, OR max of Brem spectra
! P,E       is scattered particle  momentum,energy
! THETA_DEG  is kaon angle in degrees
! RAD_LEN (%)is the radiation length of target, including internal
!                (typically 5%)
!               = .5 *(target radiation length in %) +5.
!       ***  =100. IF BREMSTRULUNG PHOTON BEAM OF 1 EQUIVIVENT QUANTA
***
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
! SIGMA      is output cross section in nanobars/GeV-str
!------------------------------------------------------------------------------

      IMPLICIT NONE       
      double precision E0MM,PMM,E0,P,THETA_DEG,RAD_LEN,SIGMA
      INTEGER TYPE
      COMMON/WISER_ALL/ E,P_COM,COST,P_T,TYPE_COM,PARTICLE,M_X,U_MAN
      double precision E,P_COM,COST,P_T,M_X,U_MAN
C      double precision E1,P_COM1,COST1,P_T1,M_X1,U_MAN1
      INTEGER TYPE_COM,PARTICLE
!  Wiser's fit    pi+     pi-    k+     k-     p+      p-   
      double precision A5(6)/-5.49,  -5.23, -5.91, -4.45, -6.77,  -6.53/
      double precision A6(6)/-1.73,  -1.82, -1.74, -3.23,  1.90,  -2.45/
      double precision MASS2(3)/.019488, .2437, .8804/
      double precision MASS(3)/.1396, .4973, .9383/ 
      double precision MP/.9383/,  MP2/.8804/, RADDEG/.0174533/
      double precision  M_L,SIG_E
      double precision E_GAMMA_MIN,WISER_ALL_FIT,QUADMO,E08
      double precision EPSILON/.003/
      EXTERNAL WISER_ALL_FIT                        
      INTEGER N,CHARGE

      P=PMM/1000.
      E0=E0MM/1000.

      P_COM = P
      TYPE_COM = TYPE
      PARTICLE = (TYPE+1)/2       ! 1= pi, 2= K, 3 =P
      CHARGE = TYPE -2*PARTICLE +2  ! 1 for + charge, 2 for - charge
      E08 =E0
                
      E =SQRT(MASS2(PARTICLE) + P**2)

      COST = COS(RADDEG * THETA_DEG)
      P_T = P * SIN(RADDEG * THETA_DEG)
      IF(TYPE.LE.4) THEN  !mesons
       IF(CHARGE.EQ.1) THEN   ! K+ n final state
        M_X = MP
       ELSE   ! K- K+ P final state
        M_X = MP+ MASS(PARTICLE)
       ENDIF
      ELSE  ! baryons 
       IF(CHARGE.EQ.1) THEN   ! pi p  final state
        M_X = MASS(1)  ! pion mass
       ELSE   ! P P-bar  P final state
        M_X = 2.*MP
       ENDIF
      ENDIF
      E_GAMMA_MIN = (M_X**2 -MASS2(PARTICLE ) -MP2+2.*MP*E)/
     >  (2.*(MP -E +P*COST))
!      WRITE(10,'(''E_GAMMA_MIN='',F10.2,''  p_t='',F8.2)')
!     >     E_GAMMA_MIN,P_T
!      E_GAMMA_MIN = MP *(E + MASS(PARTILCE))/(MP -P*(1.-COST))
      
*      print *,E_GAMMA_MIN

      
      IF(E_GAMMA_MIN.GT..1) THEN !Kinematically allowed?
       M_L = SQRT(P_T**2 + MASS2(PARTICLE))    

       IF(TYPE.NE.5) THEN  ! everything but proton
c          print *,'calling quadmo:',E_GAMMA_MIN,E08,EPSILON,N
c          print *,'quadmo=',
c     >         QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N),
c     >         EXP(A5(TYPE) *M_L), EXP(A6(TYPE) *P_T**2/E)
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) *EXP(A6(TYPE) *P_T**2/E)
       ELSE ! proton

        U_MAN = ABS(MP2 + MASS2(PARTICLE) -2.*MP*E)
c        print *,'quadmo=',
c     >       QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N),
c     >       EXP(A5(TYPE) *M_L) 
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) 
       ENDIF
       SIGMA = P**2/E * 1000. * RAD_LEN/100. *SIG_E 
      ELSE   ! Kinematically forbidden
       SIGMA = 0.
      ENDIF
c      print *,'SIGMA=',SIGMA,SIG_E,RAD_LEN,P,E

      RETURN
      END


      double precision FUNCTION WISER_ALL_FIT(E_GAMMA)

!---------------------------------------------------------
! Calculates  pi, k, p  cross section for gamma + p -> k
!  It is already divided by E_GAMMA, the bremstrulung spectra
! David Wiser's fit from Thesis, eq. IV-A-2 and Table III.
! Can be called from WISER_SIG using integration routine QUADMO
! E,P are KAON energy and momentum
! P_t is KAON transverse momentum
! P_CM is KAON center of mass momentum
! P_CM_L is KAON center of mass longitudinal momentum
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
! E_GAMMA is photon energy.
!             Steve Rock 2/21/96
!---------------------------------------------------------
                           
      IMPLICIT NONE       
      COMMON/WISER_ALL/ E,P,COST,P_T,TYPE,PARTICLE,M_X,U_MAN

      double precision  E,P,COST,P_T,M_X,U_MAN
      INTEGER  TYPE  !  1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
      INTEGER PARTICLE   ! 1= pi, 2= K, 3 =P
!  Wiser's fit    pi+     pi-    k+     k-     p+       p- 
      double precision A1(6)/566.,  486.,  368., 18.2, 1.33E5,  1.63E3 / 
      double precision A2(6)/829.,  115.,  1.91, 307., 5.69E4, -4.30E3 / 
      double precision A3(6)/1.79,  1.77,  1.91, 0.98, 1.41,    1.79 / 
      double precision A4(6)/2.10,  2.18,  1.15, 1.83,  .72,    2.24 /
      double precision A6/1.90/,A7/-.0117/ !proton only
      double precision MASS2(3)/.019488, .2437, .8804/
      double precision MASS(3)/.1396, .4973, .9383/ 
      double precision MP2/.8804/,MP/.9383/, RADDEG/.0174533/
      double precision X_R,S,B_CM, GAM_CM,  P_CM
      double precision P_CM_MAX, P_CM_L
      double precision E_GAMMA
                                            

!Mandlestam variables                                                
      S = MP2 + 2.* E_GAMMA * MP    

!Go to Center of Mass to get X_R
      B_CM = E_GAMMA/(E_GAMMA+MP)
      GAM_CM = 1./SQRT(1.-B_CM**2)
      P_CM_L = -GAM_CM *B_CM *E + 
     >          GAM_CM * P * COST
      P_CM = SQRT(P_CM_L**2 + P_T**2)  


      P_CM_MAX =SQRT (S +(M_X**2-MASS2(PARTICLE))**2/S 
     >    -2.*(M_X**2 +MASS2(PARTICLE)) )/2.
      X_R =  P_CM/P_CM_MAX   
       IF(X_R.GT.1.) THEN  ! Out of kinematic range
        WISER_ALL_FIT = 0.
       ELSEIF(TYPE.NE.5) THEN  ! not the proton
        WISER_ALL_FIT = (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)/E_GAMMA  
       ELSE ! special formula for proton
        WISER_ALL_FIT = ( (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)          /
     >   (1.+U_MAN)**(A6+A7*S) )/E_GAMMA  
       ENDIF
      
      RETURN
      END


