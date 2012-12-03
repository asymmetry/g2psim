      function x_r12p5_dent(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 49)
      data ncoeff/ 48/
      data avdat/ -0.6064667E+01/
      data xmin/
     1 -0.99947E-01,-0.49713E-01,-0.94924E-01,-0.31973E-01,-0.49984E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99987E-01, 0.55685E-01, 0.85487E-01, 0.29976E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.41662822E-02, 0.10531462E+00, 0.17303853E-04,-0.50679576E-01,
     +  0.60727708E-02, 0.30205394E-02, 0.17367827E-03,-0.97579695E-03,
     + -0.30769382E-02,-0.42611460E-03,-0.18279970E-02,-0.28059838E-03,
     +  0.27733983E-03,-0.18383824E-02,-0.27542058E-03,-0.55970112E-03,
     + -0.19668061E-02,-0.33105069E-03,-0.12086581E-02,-0.15728103E-02,
     +  0.88135610E-04, 0.31611108E-03,-0.14528136E-03, 0.33993373E-03,
     +  0.18326941E-03,-0.55294891E-03, 0.21521666E-03,-0.10267943E-02,
     + -0.29861436E-02,-0.14341164E-02,-0.78138371E-04,-0.20009637E-03,
     + -0.47101628E-04,-0.12694606E-03,-0.12780273E-03,-0.26328615E-04,
     + -0.97806146E-03,-0.15510630E-02,-0.88445295E-03,-0.86898392E-03,
     + -0.53307286E-03,-0.67665044E-03, 0.20176822E-04,-0.67487723E-04,
     + -0.16331709E-04, 0.13186960E-03,-0.13319646E-03, 0.65768370E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_r12p5_dent=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)*x11            *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x51
     8  +coeff(  8)*x11*x22            
      x_r12p5_dent=x_r12p5_dent
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)*x11            *x52
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)    *x21    *x44    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)*x12*x21            
     8  +coeff( 17)*x11    *x31*x41    
      x_r12p5_dent=x_r12p5_dent
     9  +coeff( 18)*x13                
     1  +coeff( 19)*x11    *x32        
     2  +coeff( 20)*x11*x22    *x42    
     3  +coeff( 21)*x11        *x44    
     4  +coeff( 22)    *x21    *x41    
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)    *x21*x31        
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)    *x23            
      x_r12p5_dent=x_r12p5_dent
     9  +coeff( 27)*x11    *x31        
     1  +coeff( 28)*x11        *x42    
     2  +coeff( 29)*x11*x22*x31*x41    
     3  +coeff( 30)*x11*x22*x32        
     4  +coeff( 31)    *x22            
     5  +coeff( 32)        *x31*x41    
     6  +coeff( 33)*x12                
     7  +coeff( 34)        *x32        
     8  +coeff( 35)    *x22    *x42    
      x_r12p5_dent=x_r12p5_dent
     9  +coeff( 36)            *x44    
     1  +coeff( 37)    *x23*x31*x41    
     2  +coeff( 38)*x12*x21*x31*x41    
     3  +coeff( 39)*x12*x21*x32        
     4  +coeff( 40)*x12*x23    *x42    
     5  +coeff( 41)    *x23*x32*x42    
     6  +coeff( 42)*x12*x21    *x44    
     7  +coeff( 43)        *x31        
     8  +coeff( 44)            *x42    
      x_r12p5_dent=x_r12p5_dent
     9  +coeff( 45)                *x52
     1  +coeff( 46)*x11*x22    *x41    
     2  +coeff( 47)    *x22*x31*x41    
     3  +coeff( 48)    *x21*x31*x41*x51
c
      return
      end
      function y_r12p5_dent(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 32)
      data ncoeff/ 31/
      data avdat/ -0.2526331E-02/
      data xmin/
     1 -0.99947E-01,-0.49713E-01,-0.94924E-01,-0.31973E-01,-0.49984E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99987E-01, 0.55685E-01, 0.85487E-01, 0.29976E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.58295676E-02, 0.35715722E-04, 0.14128375E+00, 0.74085198E-01,
     +  0.10906773E-01, 0.75318869E-02,-0.74728049E-03,-0.52880039E-02,
     +  0.37876570E-02,-0.19549190E-02,-0.64735295E-03,-0.24689070E-02,
     + -0.12205956E-02,-0.92498753E-02, 0.25475046E-03,-0.42712479E-03,
     + -0.46904647E-03,-0.72990649E-03,-0.32215442E-02,-0.67934883E-03,
     +  0.13349038E-03,-0.30652122E-02,-0.99690421E-03,-0.11927497E-02,
     + -0.31696856E-02,-0.49222563E-03,-0.15605809E-02,-0.11635218E-02,
     + -0.15283386E-02, 0.62234760E-02, 0.23228386E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_r12p5_dent=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)        *x31        
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)        *x31    *x51
     7  +coeff(  7)                *x51
     8  +coeff(  8)*x12*x22*x31        
      y_r12p5_dent=y_r12p5_dent
     9  +coeff(  9)*x12*x24    *x41    
     1  +coeff( 10)    *x22*x31        
     2  +coeff( 11)*x11*x21*x31        
     3  +coeff( 12)    *x24    *x41    
     4  +coeff( 13)    *x24*x31        
     5  +coeff( 14)*x12*x22    *x41    
     6  +coeff( 15)    *x22            
     7  +coeff( 16)    *x21    *x41    
     8  +coeff( 17)    *x21*x31        
      y_r12p5_dent=y_r12p5_dent
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)        *x31*x42    
     2  +coeff( 20)        *x31    *x52
     3  +coeff( 21)*x12                
     4  +coeff( 22)        *x32*x41    
     5  +coeff( 23)        *x33        
     6  +coeff( 24)*x12    *x31        
     7  +coeff( 25)    *x22    *x43    
     8  +coeff( 26)            *x45    
      y_r12p5_dent=y_r12p5_dent
     9  +coeff( 27)*x12        *x43    
     1  +coeff( 28)*x13*x21    *x41    
     2  +coeff( 29)*x13*x21*x31        
     3  +coeff( 30)*x12*x22    *x43    
     4  +coeff( 31)    *x24    *x45    
c
      return
      end
      function t_r12p5_dent(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/  0.1310190E+01/
      data xmin/
     1 -0.99947E-01,-0.49713E-01,-0.94924E-01,-0.31973E-01,-0.49984E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99987E-01, 0.55685E-01, 0.85487E-01, 0.29976E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.24501546E-02,-0.44419128E-01, 0.65398850E-01, 0.43605166E-02,
     +  0.14145293E-02, 0.80602862E-04,-0.26640936E-02,-0.15590950E-02,
     + -0.32423344E-03, 0.14860311E-03,-0.16693911E-03,-0.40861036E-03,
     + -0.75910857E-03,-0.14036569E-02,-0.15160076E-02,-0.74190582E-03,
     + -0.19089725E-03,-0.88271692E-04,-0.44840269E-03, 0.92822747E-05,
     + -0.77769007E-04, 0.14146890E-03, 0.25976932E-03, 0.12637513E-03,
     +  0.24992594E-03,-0.17275702E-03,-0.58023341E-03,-0.47440629E-03,
     + -0.86435099E-03,-0.58532547E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_r12p5_dent=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x11            *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x51
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)    *x21    *x42    
      t_r12p5_dent=t_r12p5_dent
     9  +coeff(  9)*x11            *x52
     1  +coeff( 10)*x12*x23            
     2  +coeff( 11)    *x23*x32        
     3  +coeff( 12)    *x23*x31*x41    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)*x11    *x31*x41    
     7  +coeff( 16)*x11        *x42    
     8  +coeff( 17)    *x21        *x52
      t_r12p5_dent=t_r12p5_dent
     9  +coeff( 18)*x13    *x32        
     1  +coeff( 19)*x12*x21*x32        
     2  +coeff( 20)*x13    *x31*x41    
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)*x11    *x31        
     5  +coeff( 23)    *x21*x31        
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)    *x21    *x41    
     8  +coeff( 26)*x13                
      t_r12p5_dent=t_r12p5_dent
     9  +coeff( 27)*x12*x21            
     1  +coeff( 28)    *x23            
     2  +coeff( 29)*x11    *x32        
     3  +coeff( 30)*x12*x21*x31*x41    
c
      return
      end
      function p_r12p5_dent(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.4289745E-02/
      data xmin/
     1 -0.99947E-01,-0.49713E-01,-0.94924E-01,-0.31973E-01,-0.49984E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99987E-01, 0.55685E-01, 0.85487E-01, 0.29976E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.41782409E-02, 0.21590858E-02,-0.11584162E+00,-0.73559783E-01,
     + -0.30013179E-01,-0.19284464E-01, 0.12350717E-01, 0.13807215E-01,
     +  0.20899454E-01, 0.13482396E-01,-0.10603587E-02,-0.15012794E-02,
     +  0.11428604E-01, 0.78551462E-02,-0.92866467E-02,-0.68537551E-02,
     + -0.40260130E-02,-0.35378449E-02,-0.42975741E-02,-0.31908497E-02,
     + -0.86087745E-03, 0.27640522E-02, 0.22847191E-03, 0.67099318E-03,
     +  0.34176267E-03, 0.33432984E-03, 0.21658884E-02,-0.91158494E-03,
     + -0.88495942E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_r12p5_dent=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      p_r12p5_dent=p_r12p5_dent
     9  +coeff(  9)*x11    *x31        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)                *x51
     3  +coeff( 12)*x11                
     4  +coeff( 13)*x11*x21*x31        
     5  +coeff( 14)*x11*x21    *x41    
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)*x11    *x31    *x51
      p_r12p5_dent=p_r12p5_dent
     9  +coeff( 18)*x11        *x41*x51
     1  +coeff( 19)*x12    *x31        
     2  +coeff( 20)*x12        *x41    
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)    *x21    *x41*x51
     5  +coeff( 23)    *x23*x31    *x51
     6  +coeff( 24)    *x22            
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)*x12                
      p_r12p5_dent=p_r12p5_dent
     9  +coeff( 27)    *x21*x31    *x51
     1  +coeff( 28)        *x31    *x52
     2  +coeff( 29)            *x41*x52
c
      return
      end
      function l_r12p5_dent(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.5422001E-02/
      data xmin/
     1 -0.99947E-01,-0.49713E-01,-0.94924E-01,-0.31973E-01,-0.49984E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99987E-01, 0.55685E-01, 0.85487E-01, 0.29976E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.53719045E-02,-0.12681403E-02, 0.59499597E-03, 0.31066129E-05,
     + -0.11296324E-01, 0.50091357E-05,-0.73601538E-02,-0.55968571E-02,
     + -0.53124921E-04, 0.13315293E-01,-0.10660241E-05,-0.15047852E-05,
     + -0.62998715E-02,-0.46651101E-04, 0.11675342E-04, 0.10819367E-04,
     + -0.46124891E-03, 0.38712692E-05, 0.76009069E-05, 0.26166799E-05,
     + -0.70209318E-03,-0.26119785E-05, 0.10222280E-02, 0.59656768E-05,
     +  0.34030495E-03,-0.11044565E-04,-0.64759603E-04,-0.10054189E-03,
     + -0.89027199E-04,-0.25654135E-04, 0.19425543E-05, 0.48592642E-05,
     + -0.19195579E-04,-0.26733604E-04, 0.74435532E-03, 0.74782566E-03,
     + -0.34236722E-02, 0.98305603E-03, 0.31621643E-04, 0.18560946E-04,
     +  0.50090736E-03,-0.83304018E-04,-0.76010801E-04,-0.75378074E-06,
     + -0.35933757E-04, 0.38740662E-03, 0.60495472E-03, 0.31563305E-03,
     + -0.12654946E-03,-0.46048972E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_r12p5_dent=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      l_r12p5_dent=l_r12p5_dent
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)*x11    *x31        
     3  +coeff( 12)*x11        *x41    
     4  +coeff( 13)*x12                
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)        *x31*x42    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)    *x22        *x51
      l_r12p5_dent=l_r12p5_dent
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)*x11    *x32        
     2  +coeff( 20)*x11        *x42    
     3  +coeff( 21)*x11*x21        *x51
     4  +coeff( 22)*x11            *x52
     5  +coeff( 23)*x12            *x51
     6  +coeff( 24)*x13                
     7  +coeff( 25)    *x22*x32        
     8  +coeff( 26)        *x34        
      l_r12p5_dent=l_r12p5_dent
     9  +coeff( 27)        *x33*x41    
     1  +coeff( 28)        *x32*x42    
     2  +coeff( 29)        *x31*x43    
     3  +coeff( 30)            *x44    
     4  +coeff( 31)        *x34*x41    
     5  +coeff( 32)    *x22    *x43    
     6  +coeff( 33)*x13    *x32        
     7  +coeff( 34)    *x22*x34        
     8  +coeff( 35)            *x41    
      l_r12p5_dent=l_r12p5_dent
     9  +coeff( 36)*x11                
     1  +coeff( 37)        *x32        
     2  +coeff( 38)        *x31*x41*x51
     3  +coeff( 39)        *x34    *x51
     4  +coeff( 40)            *x44*x51
     5  +coeff( 41)        *x32    *x51
     6  +coeff( 42)        *x31    *x51
     7  +coeff( 43)            *x41*x51
     8  +coeff( 44)                *x52
      l_r12p5_dent=l_r12p5_dent
     9  +coeff( 45)*x11            *x51
     1  +coeff( 46)            *x42*x51
     2  +coeff( 47)    *x22*x31*x41    
     3  +coeff( 48)    *x22    *x42    
     4  +coeff( 49)*x12            *x52
     5  +coeff( 50)    *x22*x31        
c
      return
      end
      function x_r12p5_dext(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 68)
      data ncoeff/ 67/
      data avdat/  0.1672101E-01/
      data xmin/
     1 -0.99947E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99982E-01, 0.52704E-01, 0.83287E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.11130186E-02, 0.39714050E+00, 0.50781399E-03, 0.13449566E+00,
     + -0.27693546E+00, 0.37461691E-01, 0.93446365E-02, 0.14209104E-01,
     + -0.12934072E-01,-0.41329563E-02,-0.34151778E-02,-0.81338780E-02,
     +  0.28413257E-02,-0.14994295E-01,-0.91119278E-02,-0.21398361E-02,
     +  0.20023372E-02, 0.21940186E-02,-0.28628882E-02,-0.25508415E-02,
     + -0.52239629E-02,-0.57260348E-02,-0.36158122E-02, 0.40219846E-03,
     +  0.41809937E-03, 0.71716285E-03,-0.87142678E-03, 0.28834799E-02,
     + -0.17441099E-02,-0.11426685E-02,-0.89241946E-02, 0.56219922E-03,
     + -0.76720794E-03,-0.82000968E-03, 0.60604868E-03,-0.28105748E-02,
     + -0.10340179E-02,-0.15054567E-01,-0.90312287E-02, 0.87611079E-04,
     +  0.53566479E-03, 0.33413039E-02, 0.16555317E-02,-0.97707182E-03,
     + -0.31925479E-02, 0.13486426E-02,-0.97051524E-02,-0.53939242E-02,
     + -0.14331508E-01, 0.91894745E-03,-0.98528946E-02,-0.65581608E-02,
     + -0.27180356E-02,-0.17194634E-02, 0.13283133E-02, 0.25523200E-02,
     + -0.39794608E-02,-0.76600751E-02, 0.61242017E-02, 0.12213603E-03,
     +  0.15113360E-03,-0.24949785E-02, 0.43883518E-03, 0.28910700E-03,
     + -0.10184115E-02,-0.31644365E-03,-0.14641099E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_r12p5_dext=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)    *x22            
      x_r12p5_dext=x_r12p5_dext
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x12                
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)*x12*x23            
     8  +coeff( 17)    *x21    *x41    
      x_r12p5_dext=x_r12p5_dext
     9  +coeff( 18)    *x21*x31        
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)        *x32        
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)*x11    *x32        
     6  +coeff( 24)*x11        *x44    
     7  +coeff( 25)        *x31        
     8  +coeff( 26)    *x23            
      x_r12p5_dext=x_r12p5_dext
     9  +coeff( 27)    *x21        *x52
     1  +coeff( 28)*x11*x21        *x51
     2  +coeff( 29)*x12*x21            
     3  +coeff( 30)*x12            *x51
     4  +coeff( 31)*x11*x22    *x42    
     5  +coeff( 32)*x11        *x41    
     6  +coeff( 33)    *x22        *x51
     7  +coeff( 34)            *x42*x51
     8  +coeff( 35)*x11    *x31        
      x_r12p5_dext=x_r12p5_dext
     9  +coeff( 36)*x11        *x42    
     1  +coeff( 37)*x13                
     2  +coeff( 38)*x11*x22*x31*x41    
     3  +coeff( 39)*x11*x22*x32        
     4  +coeff( 40)            *x41*x51
     5  +coeff( 41)        *x32    *x51
     6  +coeff( 42)*x11*x23            
     7  +coeff( 43)*x11*x22    *x41    
     8  +coeff( 44)*x11*x21    *x42    
      x_r12p5_dext=x_r12p5_dext
     9  +coeff( 45)*x12*x22            
     1  +coeff( 46)*x11*x22*x31        
     2  +coeff( 47)    *x21*x31*x43    
     3  +coeff( 48)*x12*x21    *x42    
     4  +coeff( 49)    *x21*x32*x42    
     5  +coeff( 50)*x11*x25            
     6  +coeff( 51)*x12*x21*x31*x41    
     7  +coeff( 52)    *x21*x33*x41    
     8  +coeff( 53)    *x25    *x42    
      x_r12p5_dext=x_r12p5_dext
     9  +coeff( 54)    *x23    *x44    
     1  +coeff( 55)    *x22*x32*x42    
     2  +coeff( 56)*x12*x21    *x43    
     3  +coeff( 57)*x12*x21*x32        
     4  +coeff( 58)    *x25*x31*x41    
     5  +coeff( 59)*x12*x21*x31*x42    
     6  +coeff( 60)                *x53
     7  +coeff( 61)*x11            *x52
     8  +coeff( 62)    *x24            
      x_r12p5_dext=x_r12p5_dext
     9  +coeff( 63)    *x23    *x41    
     1  +coeff( 64)    *x23        *x51
     2  +coeff( 65)    *x22*x31*x41    
     3  +coeff( 66)    *x22*x31    *x51
     4  +coeff( 67)    *x21    *x44    
c
      return
      end
      function y_r12p5_dext(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 55)
      data ncoeff/ 54/
      data avdat/ -0.4977679E-03/
      data xmin/
     1 -0.99947E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99982E-01, 0.52704E-01, 0.83287E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.41405382E-02, 0.40936666E-02, 0.10197769E+00,-0.39894623E-02,
     + -0.12019285E+00,-0.56514811E-01, 0.45625236E-01,-0.31233473E-01,
     +  0.36319368E-01, 0.31980135E-01, 0.18967435E-01,-0.31675950E-01,
     +  0.25578314E-01,-0.24223609E-02,-0.20738158E-01, 0.14146111E-01,
     + -0.10319786E-01,-0.73507163E-02,-0.40091095E-02,-0.30211019E-02,
     + -0.12418177E-01, 0.24004157E-02,-0.17912155E-02,-0.53562350E-02,
     + -0.37028079E-02,-0.31448784E-02, 0.82275900E-03,-0.12605377E-01,
     + -0.24525342E-02,-0.37534409E-02, 0.29278209E-03, 0.13669218E-02,
     +  0.42874570E-03, 0.30066769E-03, 0.24101313E-02, 0.11394558E-02,
     + -0.15470157E-02, 0.19689652E-02, 0.11531097E-02,-0.85609490E-02,
     + -0.11473923E-01, 0.16853251E-03,-0.53675869E-03, 0.12033522E-02,
     + -0.20526701E-02,-0.14880530E-02, 0.15170510E-02, 0.29557762E-02,
     +  0.71020536E-02,-0.33838104E-02,-0.36918242E-02,-0.54468559E-02,
     + -0.49216594E-02, 0.23365740E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_r12p5_dext=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)        *x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x21*x31        
      y_r12p5_dext=y_r12p5_dext
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)*x11    *x31        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)*x11*x21    *x41    
     5  +coeff( 14)*x11                
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)*x11*x21*x31        
     8  +coeff( 17)*x12        *x41    
      y_r12p5_dext=y_r12p5_dext
     9  +coeff( 18)*x12    *x31        
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)        *x31*x42    
     4  +coeff( 22)    *x22            
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)            *x43    
     7  +coeff( 25)    *x21*x31    *x51
     8  +coeff( 26)        *x31    *x52
      y_r12p5_dext=y_r12p5_dext
     9  +coeff( 27)*x12                
     1  +coeff( 28)        *x32*x41    
     2  +coeff( 29)*x11        *x41*x51
     3  +coeff( 30)        *x33        
     4  +coeff( 31)        *x33*x41    
     5  +coeff( 32)            *x42    
     6  +coeff( 33)    *x21        *x51
     7  +coeff( 34)                *x52
     8  +coeff( 35)        *x31*x41    
      y_r12p5_dext=y_r12p5_dext
     9  +coeff( 36)        *x32        
     1  +coeff( 37)    *x22    *x41*x51
     2  +coeff( 38)*x12        *x41*x51
     3  +coeff( 39)*x12    *x31    *x51
     4  +coeff( 40)*x12*x22    *x41    
     5  +coeff( 41)*x12*x22*x31        
     6  +coeff( 42)*x11            *x51
     7  +coeff( 43)*x11    *x31    *x51
     8  +coeff( 44)        *x32*x41*x51
      y_r12p5_dext=y_r12p5_dext
     9  +coeff( 45)*x11*x21    *x41*x51
     1  +coeff( 46)*x11*x21*x31    *x51
     2  +coeff( 47)*x12*x22            
     3  +coeff( 48)*x11*x21    *x43    
     4  +coeff( 49)*x11*x21*x31*x42    
     5  +coeff( 50)*x12    *x31*x42    
     6  +coeff( 51)*x13*x21    *x41    
     7  +coeff( 52)*x13*x21*x31        
     8  +coeff( 53)    *x24*x31*x42    
      y_r12p5_dext=y_r12p5_dext
     9  +coeff( 54)*x12*x21    *x41*x52
c
      return
      end
      function t_r12p5_dext(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/  0.5274432E+00/
      data xmin/
     1 -0.99947E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99982E-01, 0.52704E-01, 0.83287E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.39721131E-02, 0.45008045E-01,-0.85493378E-01, 0.22742946E-01,
     +  0.76307096E-02,-0.50473628E-02,-0.53054374E-02,-0.27935749E-02,
     +  0.27926662E-02,-0.18097078E-02,-0.24349659E-02,-0.12650673E-02,
     + -0.12501525E-02, 0.19525759E-02, 0.17465726E-02, 0.19936387E-02,
     +  0.27886287E-02,-0.14241460E-03, 0.12701957E-02,-0.88807353E-03,
     +  0.39142934E-02, 0.21455449E-03,-0.29791132E-03, 0.68198651E-03,
     + -0.71158411E-03, 0.81687403E-03,-0.72254872E-04,-0.27038937E-03,
     + -0.17829755E-03,-0.17884291E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_r12p5_dext=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11*x21            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)*x12                
      t_r12p5_dext=t_r12p5_dext
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x32        
     3  +coeff( 12)                *x52
     4  +coeff( 13)            *x42    
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)*x11    *x32        
     7  +coeff( 16)*x11    *x31*x41    
     8  +coeff( 17)    *x21    *x42    
      t_r12p5_dext=t_r12p5_dext
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)        *x32    *x51
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)    *x23*x31*x41    
     4  +coeff( 22)        *x31        
     5  +coeff( 23)    *x21    *x41    
     6  +coeff( 24)*x12            *x51
     7  +coeff( 25)*x11*x21        *x51
     8  +coeff( 26)        *x31*x41*x51
      t_r12p5_dext=t_r12p5_dext
     9  +coeff( 27)            *x41    
     1  +coeff( 28)*x11    *x31        
     2  +coeff( 29)*x11        *x41    
     3  +coeff( 30)        *x31    *x51
c
      return
      end
      function p_r12p5_dext(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 36)
      data ncoeff/ 35/
      data avdat/  0.3024764E-03/
      data xmin/
     1 -0.99947E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99982E-01, 0.52704E-01, 0.83287E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19089943E-02, 0.85264439E-03,-0.38246524E-01, 0.59923758E-02,
     + -0.57363762E-02,-0.12691601E-01, 0.64168023E-02, 0.10419062E-01,
     +  0.40504495E-02, 0.66471086E-02,-0.80682151E-02, 0.53969965E-04,
     + -0.91643642E-05, 0.73573110E-02,-0.30599646E-02,-0.78091741E-03,
     + -0.51000173E-03,-0.24311848E-02,-0.22671756E-02,-0.29234614E-02,
     + -0.85600402E-03, 0.47508421E-03,-0.30929199E-03,-0.21940530E-02,
     + -0.11218145E-02,-0.35560626E-03,-0.73987636E-03,-0.10606889E-02,
     +  0.61437158E-05,-0.59008942E-03, 0.35398162E-03, 0.22623004E-03,
     +  0.11003536E-03, 0.14063156E-03,-0.59127103E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_r12p5_dext=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      p_r12p5_dext=p_r12p5_dext
     9  +coeff(  9)*x11    *x31        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)                *x53
     5  +coeff( 14)*x11*x21    *x41    
     6  +coeff( 15)*x12        *x41    
     7  +coeff( 16)                *x51
     8  +coeff( 17)*x11                
      p_r12p5_dext=p_r12p5_dext
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)        *x31*x42    
     3  +coeff( 21)*x11    *x31    *x51
     4  +coeff( 22)    *x22            
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)        *x32*x41    
     7  +coeff( 25)            *x43    
     8  +coeff( 26)            *x41*x52
      p_r12p5_dext=p_r12p5_dext
     9  +coeff( 27)*x12    *x31        
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)*x12        *x42    
     3  +coeff( 30)    *x22*x33        
     4  +coeff( 31)        *x31*x41    
     5  +coeff( 32)            *x42    
     6  +coeff( 33)    *x21        *x51
     7  +coeff( 34)*x12                
     8  +coeff( 35)        *x33        
c
      return
      end
      function l_r12p5_dext(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.3434530E-01/
      data xmin/
     1 -0.99947E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99982E-01, 0.52704E-01, 0.83287E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13294291E-01,-0.55474007E+00, 0.78331126E-03, 0.61478460E-03,
     + -0.98333530E-01, 0.35459375E+00,-0.32980647E-01,-0.24037866E-01,
     + -0.38608804E-01, 0.44352792E-01,-0.15012405E-01,-0.41072289E-02,
     + -0.55123544E-02, 0.24203338E-01, 0.13590524E-01, 0.14412630E-01,
     + -0.22150269E-02,-0.20687983E-02,-0.26981251E-02, 0.15401280E-02,
     + -0.37156856E-02, 0.11442479E-01,-0.18023783E-02, 0.55833552E-02,
     +  0.10152202E-01, 0.45684515E-02,-0.41054352E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_r12p5_dext=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)*x11            *x51
      l_r12p5_dext=l_r12p5_dext
     9  +coeff(  9)    *x22            
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)*x12                
     3  +coeff( 12)        *x32        
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x21*x34        
      l_r12p5_dext=l_r12p5_dext
     9  +coeff( 18)    *x21    *x41    
     1  +coeff( 19)            *x42    
     2  +coeff( 20)                *x52
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x21*x32        
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)*x11    *x32        
     7  +coeff( 25)*x11    *x31*x41    
     8  +coeff( 26)*x11        *x42    
      l_r12p5_dext=l_r12p5_dext
     9  +coeff( 27)*x11*x21        *x51
c
      return
      end
      function x_r12p5_fp  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1823898E-01/
      data xmin/
     1 -0.99237E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49017E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99987E-01, 0.52704E-01, 0.94327E-01, 0.31909E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28344881E-01, 0.11976781E+00,-0.27404734E+00, 0.70601147E+00,
     +  0.58173854E-02, 0.59324998E-01,-0.68592136E-02,-0.58155622E-01,
     +  0.21748831E-01,-0.12026033E-01,-0.22558395E-01,-0.10719726E-01,
     + -0.22506585E-01, 0.88607455E-02, 0.47850520E-02, 0.52766590E-02,
     + -0.12895361E-01,-0.49434723E-02,-0.31570750E-02, 0.31364136E-02,
     + -0.72184377E-02, 0.84467251E-02,-0.16757157E-01,-0.36273573E-02,
     +  0.19647735E-02, 0.12649618E-02,-0.36661737E-02,-0.15786337E-01,
     + -0.17840832E-01, 0.38611181E-01,-0.13397729E-01,-0.72459192E-02,
     +  0.84727379E-02,-0.35464261E-01, 0.14242684E-01, 0.47280970E-02,
     +  0.25761961E-02, 0.55435364E-03, 0.27537073E-02,-0.12419029E-01,
     +  0.13928855E-01, 0.18114801E-01, 0.15488710E-01, 0.87623624E-02,
     + -0.23044704E-02,-0.80273105E-02, 0.63768262E-02,-0.84787868E-02,
     +  0.70622708E-02,-0.97154971E-03, 0.48270710E-02,-0.16925791E-01,
     + -0.18604904E-01,-0.13007276E-03,-0.98440680E-03,-0.78448644E-02,
     + -0.45952175E-03,-0.61833053E-02,-0.19819669E-02,-0.13088265E-01,
     + -0.55345711E-02,-0.20972444E-02, 0.86599393E-02,-0.10012769E-02,
     + -0.16573925E-01, 0.57622994E-03, 0.21931180E-02,-0.12540733E-01,
     + -0.41019046E-02, 0.75385743E-02,-0.18166417E-02,-0.27441055E-01,
     + -0.11513591E-01,-0.66267503E-02, 0.56208195E-02,-0.38546186E-02,
     +  0.11301106E-01,-0.12303248E-01, 0.59906091E-02,-0.14042138E-02,
     + -0.77787312E-02, 0.67390562E-02,-0.10695273E-01, 0.12747862E-02,
     +  0.28236448E-02, 0.24084421E-02,-0.31349966E-02, 0.89909119E-03,
     + -0.11745433E-01,-0.16347173E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x14 = x13*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x26 = x25*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      x_r12p5_fp  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11*x21            
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)                *x52
      x_r12p5_fp  =x_r12p5_fp  
     9  +coeff(  9)    *x22            
     1  +coeff( 10)*x12                
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11        *x42    
     6  +coeff( 15)*x11            *x52
     7  +coeff( 16)                *x53
     8  +coeff( 17)    *x21*x34        
      x_r12p5_fp  =x_r12p5_fp  
     9  +coeff( 18)        *x34    *x51
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11    *x32        
     5  +coeff( 23)*x11*x21    *x42    
     6  +coeff( 24)    *x24    *x42    
     7  +coeff( 25)*x13    *x31*x41    
     8  +coeff( 26)    *x23            
      x_r12p5_fp  =x_r12p5_fp  
     9  +coeff( 27)        *x32        
     1  +coeff( 28)    *x24            
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)*x11    *x31*x41    
     5  +coeff( 32)        *x32    *x51
     6  +coeff( 33)        *x31*x41*x51
     7  +coeff( 34)*x12*x22            
     8  +coeff( 35)*x13*x21            
      x_r12p5_fp  =x_r12p5_fp  
     9  +coeff( 36)*x12        *x42    
     1  +coeff( 37)        *x31*x41*x52
     2  +coeff( 38)*x11*x21        *x51
     3  +coeff( 39)*x12            *x51
     4  +coeff( 40)    *x22*x31*x41    
     5  +coeff( 41)    *x22    *x42    
     6  +coeff( 42)*x11*x22        *x51
     7  +coeff( 43)*x11*x21*x31*x41    
     8  +coeff( 44)    *x21*x32    *x51
      x_r12p5_fp  =x_r12p5_fp  
     9  +coeff( 45)*x14                
     1  +coeff( 46)*x12    *x31*x41    
     2  +coeff( 47)*x11*x23        *x51
     3  +coeff( 48)*x11    *x32    *x51
     4  +coeff( 49)        *x32    *x52
     5  +coeff( 50)                *x54
     6  +coeff( 51)*x13*x22            
     7  +coeff( 52)*x11*x22*x32        
     8  +coeff( 53)*x12*x21*x31*x41    
      x_r12p5_fp  =x_r12p5_fp  
     9  +coeff( 54)*x12*x21            
     1  +coeff( 55)*x13                
     2  +coeff( 56)    *x23        *x51
     3  +coeff( 57)            *x42*x51
     4  +coeff( 58)    *x22        *x52
     5  +coeff( 59)*x11*x21*x32        
     6  +coeff( 60)*x12*x21        *x51
     7  +coeff( 61)    *x21*x31*x41*x51
     8  +coeff( 62)    *x21    *x42*x51
      x_r12p5_fp  =x_r12p5_fp  
     9  +coeff( 63)*x11*x21        *x52
     1  +coeff( 64)    *x21        *x53
     2  +coeff( 65)    *x23*x31*x41    
     3  +coeff( 66)*x11            *x53
     4  +coeff( 67)*x11*x22*x31*x41    
     5  +coeff( 68)*x11*x22    *x42    
     6  +coeff( 69)    *x22*x31*x41*x51
     7  +coeff( 70)    *x22    *x42*x51
     8  +coeff( 71)*x12*x21    *x42    
      x_r12p5_fp  =x_r12p5_fp  
     9  +coeff( 72)    *x21*x32*x42    
     1  +coeff( 73)    *x21*x31*x43    
     2  +coeff( 74)*x11*x21    *x42*x51
     3  +coeff( 75)*x11    *x34        
     4  +coeff( 76)    *x26        *x51
     5  +coeff( 77)    *x22*x32*x42    
     6  +coeff( 78)*x14*x23            
     7  +coeff( 79)*x12    *x33*x41    
     8  +coeff( 80)            *x42*x54
      x_r12p5_fp  =x_r12p5_fp  
     9  +coeff( 81)*x14*x22        *x51
     1  +coeff( 82)*x11*x22    *x42*x52
     2  +coeff( 83)*x14*x22    *x42    
     3  +coeff( 84)        *x33*x41    
     4  +coeff( 85)*x13            *x51
     5  +coeff( 86)*x11    *x31*x41*x51
     6  +coeff( 87)*x12            *x52
     7  +coeff( 88)    *x22        *x53
     8  +coeff( 89)    *x21*x33*x41    
      x_r12p5_fp  =x_r12p5_fp  
     9  +coeff( 90)    *x21    *x44    
c
      return
      end
      function y_r12p5_fp  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 98)
      data ncoeff/ 97/
      data avdat/ -0.1257272E-03/
      data xmin/
     1 -0.99237E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49017E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99987E-01, 0.52704E-01, 0.81477E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.70413121E-03,-0.73400736E-01, 0.29802712E-01, 0.26972115E-01,
     + -0.72286953E-02, 0.35419851E-02,-0.38922828E-01,-0.16138865E-01,
     +  0.32166988E-03, 0.60319573E-04, 0.37077849E-02, 0.59886620E-03,
     +  0.18062619E-01,-0.13073058E-02, 0.30708977E-02, 0.66853210E-03,
     +  0.58234790E-02,-0.82049528E-02, 0.14738891E-01, 0.89289201E-02,
     + -0.14499552E-01,-0.10547183E-01, 0.86358860E-02,-0.14242555E-02,
     + -0.15623432E-02, 0.87999366E-02, 0.82377326E-02, 0.45457338E-02,
     + -0.17221252E-02,-0.13590853E-02, 0.11500124E-02, 0.19340317E-02,
     + -0.60023769E-03, 0.91006921E-03,-0.17638458E-02,-0.34050883E-02,
     +  0.33852304E-02,-0.21871813E-02, 0.39783902E-02, 0.37892533E-02,
     +  0.12400693E-02,-0.24279552E-02,-0.71900914E-03,-0.11860685E-02,
     + -0.11273727E-02, 0.22591525E-02,-0.90257038E-03,-0.22910668E-02,
     +  0.14243281E-02,-0.12095210E-02,-0.57707638E-02,-0.48988373E-02,
     +  0.30234784E-02,-0.15758061E-02, 0.56872312E-02,-0.13782980E-02,
     +  0.55335551E-02,-0.22842451E-02,-0.32461283E-02, 0.53585707E-02,
     +  0.12780557E-02,-0.31129974E-02,-0.39559314E-02,-0.14194679E-01,
     +  0.18228021E-01, 0.36498049E-03, 0.30711279E-03, 0.51861972E-03,
     + -0.74398343E-03, 0.13373510E-02,-0.13354648E-02,-0.20303439E-03,
     + -0.38983754E-03,-0.38378157E-02, 0.76289807E-03,-0.60620427E-03,
     +  0.22353502E-02,-0.11931690E-02, 0.15842323E-02, 0.18025823E-03,
     + -0.21151288E-02,-0.90624706E-03,-0.72886730E-02, 0.65376732E-03,
     +  0.35722980E-02, 0.30844202E-02, 0.15512578E-02,-0.18678886E-02,
     + -0.64660050E-03,-0.23619852E-02, 0.38324115E-02, 0.42381026E-02,
     + -0.19158658E-02,-0.37050480E-02,-0.24936658E-02,-0.15034527E-02,
     +  0.26307802E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
      x55 = x54*x5
c
c                  function
c
      y_r12p5_fp  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)            *x41*x51
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)*x11        *x41    
      y_r12p5_fp  =y_r12p5_fp  
     9  +coeff(  9)    *x22        *x51
     1  +coeff( 10)            *x42*x51
     2  +coeff( 11)            *x41*x52
     3  +coeff( 12)                *x53
     4  +coeff( 13)        *x31    *x52
     5  +coeff( 14)    *x21            
     6  +coeff( 15)                *x51
     7  +coeff( 16)*x11                
     8  +coeff( 17)    *x22    *x41    
      y_r12p5_fp  =y_r12p5_fp  
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)    *x22*x31        
     2  +coeff( 20)*x11        *x41*x51
     3  +coeff( 21)*x11*x21*x31        
     4  +coeff( 22)*x11    *x31    *x51
     5  +coeff( 23)*x12    *x31        
     6  +coeff( 24)    *x23*x31    *x51
     7  +coeff( 25)                *x52
     8  +coeff( 26)        *x31*x42    
      y_r12p5_fp  =y_r12p5_fp  
     9  +coeff( 27)    *x21*x31    *x51
     1  +coeff( 28)        *x33        
     2  +coeff( 29)        *x31    *x53
     3  +coeff( 30)    *x22            
     4  +coeff( 31)*x11*x21            
     5  +coeff( 32)            *x43    
     6  +coeff( 33)*x12                
     7  +coeff( 34)*x11*x21    *x41    
     8  +coeff( 35)    *x22    *x41*x51
      y_r12p5_fp  =y_r12p5_fp  
     9  +coeff( 36)            *x41*x53
     1  +coeff( 37)    *x21*x33        
     2  +coeff( 38)        *x33    *x51
     3  +coeff( 39)*x12*x21    *x41    
     4  +coeff( 40)        *x32*x43*x51
     5  +coeff( 41)*x13*x22    *x41    
     6  +coeff( 42)*x13*x22*x31        
     7  +coeff( 43)            *x42    
     8  +coeff( 44)        *x31*x41    
      y_r12p5_fp  =y_r12p5_fp  
     9  +coeff( 45)        *x32        
     1  +coeff( 46)        *x32*x41    
     2  +coeff( 47)    *x23*x31        
     3  +coeff( 48)    *x21*x31*x42    
     4  +coeff( 49)    *x21*x31    *x52
     5  +coeff( 50)        *x32*x42    
     6  +coeff( 51)*x11*x22    *x41    
     7  +coeff( 52)    *x21    *x41*x53
     8  +coeff( 53)        *x31    *x54
      y_r12p5_fp  =y_r12p5_fp  
     9  +coeff( 54)*x11    *x33        
     1  +coeff( 55)*x11        *x41*x53
     2  +coeff( 56)*x13        *x41    
     3  +coeff( 57)*x12*x22    *x41    
     4  +coeff( 58)*x11    *x33    *x51
     5  +coeff( 59)        *x32*x41*x53
     6  +coeff( 60)*x12*x22*x31        
     7  +coeff( 61)*x11    *x34        
     8  +coeff( 62)*x11*x23*x31    *x51
      y_r12p5_fp  =y_r12p5_fp  
     9  +coeff( 63)*x12    *x33    *x51
     1  +coeff( 64)*x11*x21*x33    *x52
     2  +coeff( 65)*x11*x21*x32*x41*x53
     3  +coeff( 66)*x11            *x51
     4  +coeff( 67)    *x21    *x42    
     5  +coeff( 68)        *x31*x41*x51
     6  +coeff( 69)    *x21*x32        
     7  +coeff( 70)    *x23    *x41    
     8  +coeff( 71)            *x43*x51
      y_r12p5_fp  =y_r12p5_fp  
     9  +coeff( 72)*x12        *x41    
     1  +coeff( 73)    *x21*x32    *x51
     2  +coeff( 74)        *x32*x41*x51
     3  +coeff( 75)        *x32    *x52
     4  +coeff( 76)*x11        *x43    
     5  +coeff( 77)*x11*x21*x31    *x51
     6  +coeff( 78)*x11    *x31*x41*x51
     7  +coeff( 79)*x11    *x31    *x52
     8  +coeff( 80)    *x21    *x43*x51
      y_r12p5_fp  =y_r12p5_fp  
     9  +coeff( 81)            *x41*x54
     1  +coeff( 82)*x11    *x32    *x51
     2  +coeff( 83)    *x21*x32*x41*x51
     3  +coeff( 84)*x12    *x31*x41    
     4  +coeff( 85)*x11*x21    *x41*x52
     5  +coeff( 86)    *x22*x33        
     6  +coeff( 87)        *x33    *x52
     7  +coeff( 88)*x11    *x31    *x53
     8  +coeff( 89)*x13*x21            
      y_r12p5_fp  =y_r12p5_fp  
     9  +coeff( 90)    *x24    *x41*x51
     1  +coeff( 91)*x11*x21*x32*x41    
     2  +coeff( 92)*x11    *x32*x41*x51
     3  +coeff( 93)    *x24*x31    *x51
     4  +coeff( 94)        *x31    *x55
     5  +coeff( 95)*x12        *x41*x52
     6  +coeff( 96)*x11*x21*x33        
     7  +coeff( 97)*x12*x21*x31    *x51
c
      return
      end
      function t_r12p5_fp  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/ -0.5195621E-02/
      data xmin/
     1 -0.99237E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49017E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99987E-01, 0.52704E-01, 0.81477E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.59245983E-02,-0.33418745E-01,-0.39874506E-02, 0.16452523E-03,
     +  0.38504781E-03, 0.12491446E+00, 0.58122114E-02,-0.91456366E-03,
     +  0.73195365E-02,-0.11924954E-01,-0.36180613E-02,-0.40109474E-02,
     + -0.11783703E-02,-0.36060736E-02, 0.20005826E-02,-0.16428060E-02,
     + -0.57809841E-03, 0.49336208E-03, 0.17573833E-02,-0.34633044E-02,
     +  0.19919903E-02,-0.21586744E-02, 0.61252067E-03,-0.30658140E-02,
     +  0.18727503E-02,-0.20981538E-02, 0.46539874E-03, 0.21451685E-03,
     +  0.21217811E-03,-0.68919692E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_r12p5_fp  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)        *x31        
     5  +coeff(  5)            *x41    
     6  +coeff(  6)                *x51
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)*x11            *x51
      t_r12p5_fp  =t_r12p5_fp  
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x12                
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)*x11            *x52
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)*x13        *x42    
      t_r12p5_fp  =t_r12p5_fp  
     9  +coeff( 18)    *x23            
     1  +coeff( 19)*x11    *x32        
     2  +coeff( 20)    *x21*x32        
     3  +coeff( 21)*x11        *x42    
     4  +coeff( 22)        *x32    *x51
     5  +coeff( 23)*x11*x23            
     6  +coeff( 24)*x11*x21    *x42    
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)*x13    *x31*x41    
      t_r12p5_fp  =t_r12p5_fp  
     9  +coeff( 27)    *x21*x31        
     1  +coeff( 28)    *x21    *x41    
     2  +coeff( 29)        *x31    *x51
     3  +coeff( 30)*x11*x22            
c
      return
      end
      function p_r12p5_fp  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 44)
      data ncoeff/ 43/
      data avdat/ -0.7074261E-03/
      data xmin/
     1 -0.99237E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49017E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99987E-01, 0.52704E-01, 0.81477E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16546138E-02, 0.55494159E-01,-0.45440771E-01, 0.11360578E-01,
     +  0.26708273E-01,-0.23103072E-03,-0.21930559E-01,-0.18971033E-01,
     + -0.65800021E-02,-0.15184822E-01, 0.49221508E-05,-0.25471491E-04,
     +  0.93644063E-04, 0.36363563E-04, 0.32800518E-02,-0.97988499E-03,
     + -0.18495584E-02, 0.22524297E-02, 0.10658820E-02, 0.10133443E-01,
     + -0.25641135E-03, 0.56515983E-02,-0.92887934E-02, 0.43206252E-02,
     + -0.12896291E-02, 0.14338053E-01, 0.64588846E-02,-0.57412274E-02,
     +  0.44588642E-02, 0.74209739E-03,-0.19544623E-04,-0.47796129E-03,
     +  0.21563764E-02, 0.50634714E-02, 0.23760998E-02, 0.29004605E-02,
     +  0.15157157E-02,-0.24124102E-02,-0.75521448E-03,-0.63239608E-03,
     + -0.98795246E-03,-0.67756092E-03,-0.40101839E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_r12p5_fp  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      p_r12p5_fp  =p_r12p5_fp  
     9  +coeff(  9)*x11    *x31        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)            *x42*x51
     5  +coeff( 14)                *x53
     6  +coeff( 15)*x11        *x41*x51
     7  +coeff( 16)    *x24    *x41    
     8  +coeff( 17)    *x21            
      p_r12p5_fp  =p_r12p5_fp  
     9  +coeff( 18)                *x51
     1  +coeff( 19)*x11                
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)        *x31    *x52
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)*x12        *x41    
     7  +coeff( 25)    *x22            
     8  +coeff( 26)    *x22    *x41    
      p_r12p5_fp  =p_r12p5_fp  
     9  +coeff( 27)        *x31*x42    
     1  +coeff( 28)*x11*x21*x31        
     2  +coeff( 29)*x12    *x31        
     3  +coeff( 30)*x11*x21            
     4  +coeff( 31)*x11            *x51
     5  +coeff( 32)*x12                
     6  +coeff( 33)        *x33        
     7  +coeff( 34)        *x32*x41    
     8  +coeff( 35)            *x43    
      p_r12p5_fp  =p_r12p5_fp  
     9  +coeff( 36)    *x21*x31    *x51
     1  +coeff( 37)            *x41*x52
     2  +coeff( 38)*x11    *x31    *x51
     3  +coeff( 39)        *x33*x41    
     4  +coeff( 40)        *x32        
     5  +coeff( 41)        *x31*x41    
     6  +coeff( 42)            *x42    
     7  +coeff( 43)                *x52
c
      return
      end
      function l_r12p5_fp  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/ -0.3774830E-01/
      data xmin/
     1 -0.99237E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49017E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99987E-01, 0.52704E-01, 0.81477E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.24858454E-01,-0.35668004E+00,-0.32989040E-01, 0.21483250E+00,
     + -0.40239185E-01,-0.33062227E-01, 0.45495387E-01,-0.66622859E-02,
     + -0.17825836E-01,-0.13883304E-01,-0.84385918E-02, 0.10594022E-01,
     +  0.65853847E-02, 0.20419473E-02,-0.36827922E-02, 0.67371614E-02,
     +  0.64424942E-02,-0.33326370E-02, 0.10941052E-01, 0.51717330E-02,
     +  0.90060476E-02,-0.77810264E-02,-0.56454964E-03, 0.65870197E-02,
     +  0.46483488E-03, 0.17870513E-02,-0.13690229E-02,-0.96746575E-03,
     + -0.30731547E-02, 0.57068011E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_r12p5_fp  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)*x11            *x51
      l_r12p5_fp  =l_r12p5_fp  
     9  +coeff(  9)*x12                
     1  +coeff( 10)        *x32        
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)                *x53
     5  +coeff( 14)        *x31        
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      l_r12p5_fp  =l_r12p5_fp  
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)*x11    *x32        
     3  +coeff( 21)*x11    *x31*x41    
     4  +coeff( 22)*x11*x21        *x51
     5  +coeff( 23)*x11            *x52
     6  +coeff( 24)*x12            *x51
     7  +coeff( 25)            *x41    
     8  +coeff( 26)        *x31*x41    
      l_r12p5_fp  =l_r12p5_fp  
     9  +coeff( 27)        *x31    *x51
     1  +coeff( 28)*x11        *x41    
     2  +coeff( 29)    *x23            
     3  +coeff( 30)    *x21*x31*x41    
c
      return
      end
      function x_r12p5_q1ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.3698746E-02/
      data xmin/
     1 -0.99959E-01,-0.49985E-01,-0.94921E-01,-0.31889E-01,-0.49871E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99982E-01, 0.56044E-01, 0.85487E-01, 0.29980E-01, 0.49753E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28088358E-02, 0.11379560E+00,-0.19078256E-06, 0.58088510E-04,
     +  0.79326451E-01, 0.10102723E-02, 0.99307415E-03,-0.10993716E-03,
     + -0.69966467E-04, 0.67561617E-04,-0.84206702E-04,-0.49470902E-04,
     + -0.51315077E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_r12p5_q1ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)*x11*x22            
      x_r12p5_q1ex=x_r12p5_q1ex
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)*x11    *x31*x41    
     2  +coeff( 11)    *x25            
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)*x11            *x52
c
      return
      end
      function y_r12p5_q1ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.3374772E-02/
      data xmin/
     1 -0.99959E-01,-0.49985E-01,-0.94921E-01,-0.31889E-01,-0.49871E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99982E-01, 0.56044E-01, 0.85487E-01, 0.29980E-01, 0.49753E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.52203899E-02, 0.91709085E-01, 0.77632663E-04, 0.11021557E+00,
     + -0.10375979E-02, 0.17532704E-04,-0.68410154E-03,-0.25807307E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      y_r12p5_q1ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)        *x31        
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x43*x51
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)*x11*x21*x31        
c
      return
      end
      function t_r12p5_q1ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/  0.1919874E-02/
      data xmin/
     1 -0.99959E-01,-0.49985E-01,-0.94921E-01,-0.31889E-01,-0.49871E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99982E-01, 0.56044E-01, 0.85487E-01, 0.29980E-01, 0.49753E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.15107320E-02,-0.40871304E-01, 0.81875976E-02, 0.19270645E-02,
     +  0.21358340E-02,-0.25546066E-02, 0.89828216E-04, 0.11946512E-03,
     + -0.15192934E-02,-0.19714760E-02,-0.87757735E-03,-0.11721627E-02,
     +  0.18329022E-05, 0.10802281E-03,-0.17551628E-03,-0.76851208E-03,
     + -0.11048208E-02,-0.56339044E-03,-0.11585184E-02, 0.46651490E-04,
     + -0.15372114E-04, 0.68698748E-04,-0.24540943E-05, 0.62251979E-05,
     + -0.81757673E-04,-0.60667357E-04, 0.19554803E-03, 0.22471670E-03,
     +  0.16168843E-03, 0.17992232E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      t_r12p5_q1ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x11            *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21*x31*x41    
     7  +coeff(  7)*x13    *x31*x41    
     8  +coeff(  8)                *x51
      t_r12p5_q1ex=t_r12p5_q1ex
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)*x11    *x31*x41    
     2  +coeff( 11)*x11        *x42    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x13    *x32        
     5  +coeff( 14)*x11*x22*x32        
     6  +coeff( 15)*x13                
     7  +coeff( 16)*x12*x21            
     8  +coeff( 17)*x11*x22            
      t_r12p5_q1ex=t_r12p5_q1ex
     9  +coeff( 18)    *x23            
     1  +coeff( 19)*x11    *x32        
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)*x11    *x33        
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)*x11    *x32*x41    
     6  +coeff( 24)*x11        *x43    
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)    *x22            
      t_r12p5_q1ex=t_r12p5_q1ex
     9  +coeff( 27)*x11    *x31        
     1  +coeff( 28)    *x21*x31        
     2  +coeff( 29)*x11        *x41    
     3  +coeff( 30)    *x21    *x41    
c
      return
      end
      function p_r12p5_q1ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  8)
      data ncoeff/  7/
      data avdat/ -0.1524253E-02/
      data xmin/
     1 -0.99959E-01,-0.49985E-01,-0.94921E-01,-0.31889E-01,-0.49871E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99982E-01, 0.56044E-01, 0.85487E-01, 0.29980E-01, 0.49753E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.26083640E-02,-0.60889863E-06, 0.42855073E-01, 0.61011359E-01,
     + -0.22571362E-02,-0.15906629E-02, 0.14312865E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      p_r12p5_q1ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x51
c
      return
      end
      function l_r12p5_q1ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/ -0.1452321E-02/
      data xmin/
     1 -0.99959E-01,-0.49985E-01,-0.94921E-01,-0.31889E-01,-0.49871E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99982E-01, 0.56044E-01, 0.85487E-01, 0.29980E-01, 0.49753E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14401743E-02,-0.31707296E-03, 0.16040257E-03, 0.80834764E-06,
     +  0.27811662E-04,-0.27758349E-02,-0.11676791E-06,-0.28830886E-03,
     +  0.24011005E-07,-0.10045568E-02,-0.17496753E-02, 0.48232079E-03,
     + -0.29741795E-03,-0.32237321E-06,-0.27096800E-06, 0.61203791E-04,
     + -0.25559510E-04, 0.70862232E-04, 0.34489545E-04, 0.26295074E-04,
     +  0.29737794E-04,-0.55063088E-05,-0.61495866E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_r12p5_q1ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      l_r12p5_q1ex=l_r12p5_q1ex
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)*x12                
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)        *x33        
     7  +coeff( 16)        *x31        
     8  +coeff( 17)    *x22        *x51
      l_r12p5_q1ex=l_r12p5_q1ex
     9  +coeff( 18)        *x31*x41*x51
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)*x12            *x51
     3  +coeff( 21)        *x32    *x51
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)            *x41*x51
c
      return
      end
      function x_r12p5_q3en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/  0.9170212E-02/
      data xmin/
     1 -0.99237E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99403E-01, 0.52704E-01, 0.82118E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21186988E-02, 0.25808421E+00, 0.13970692E+00,-0.19358210E+00,
     +  0.23390580E-01, 0.30704191E-01,-0.23155645E-01,-0.58935969E-02,
     +  0.37261830E-02, 0.53164805E-02,-0.36978663E-02,-0.24223919E-02,
     + -0.79322206E-02,-0.13586742E-01,-0.77275946E-02, 0.17169920E-02,
     +  0.44932350E-03, 0.15812689E-02, 0.16964213E-02, 0.21082766E-02,
     + -0.21772406E-02,-0.69296798E-02,-0.37255806E-02,-0.32130655E-03,
     +  0.35997829E-03,-0.84951613E-03, 0.34305584E-03,-0.10342188E-02,
     + -0.54210867E-03,-0.11867234E-02,-0.16845624E-02, 0.27795401E-03,
     + -0.75371313E-05, 0.42500824E-03, 0.48351050E-02, 0.18637025E-02,
     + -0.18944670E-02,-0.42681810E-02,-0.64324196E-02,-0.29998732E-03,
     + -0.33314866E-02, 0.23413585E-02,-0.14371771E-01,-0.28010413E-01,
     +  0.95463700E-04, 0.22373360E-02,-0.48365904E-03, 0.46634828E-03,
     +  0.33704311E-03,-0.27525160E-02, 0.31367649E-03,-0.27562315E-02,
     + -0.43333727E-02,-0.13024459E-02, 0.56471821E-03, 0.78333431E-03,
     + -0.87088934E-03, 0.12992510E-02, 0.63825701E-03,-0.26924657E-02,
     + -0.15656407E-02, 0.32805281E-02,-0.21404624E-02,-0.13217629E-01,
     +  0.16408814E-03, 0.29828423E-03, 0.20305740E-03, 0.35501219E-03,
     + -0.33035662E-03, 0.56331966E-03,-0.34623884E-03, 0.11072186E-03,
     + -0.19921472E-03, 0.46807906E-03, 0.38246650E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_r12p5_q3en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)                *x52
      x_r12p5_q3en=x_r12p5_q3en
     9  +coeff(  9)*x11            *x51
     1  +coeff( 10)*x12                
     2  +coeff( 11)            *x42    
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)*x12*x23            
     8  +coeff( 17)            *x41    
      x_r12p5_q3en=x_r12p5_q3en
     9  +coeff( 18)    *x21    *x41    
     1  +coeff( 19)    *x21*x31        
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)        *x32        
     4  +coeff( 22)*x11*x22            
     5  +coeff( 23)*x11    *x31*x41    
     6  +coeff( 24)    *x25            
     7  +coeff( 25)        *x31        
     8  +coeff( 26)            *x42*x51
      x_r12p5_q3en=x_r12p5_q3en
     9  +coeff( 27)                *x53
     1  +coeff( 28)*x11        *x42    
     2  +coeff( 29)*x12            *x51
     3  +coeff( 30)*x13                
     4  +coeff( 31)*x11    *x32        
     5  +coeff( 32)*x13        *x41    
     6  +coeff( 33)*x13    *x31        
     7  +coeff( 34)*x11    *x31        
     8  +coeff( 35)*x11*x23            
      x_r12p5_q3en=x_r12p5_q3en
     9  +coeff( 36)*x11*x22    *x41    
     1  +coeff( 37)*x11*x21    *x42    
     2  +coeff( 38)*x11*x22    *x42    
     3  +coeff( 39)*x11*x22*x31*x41    
     4  +coeff( 40)    *x24*x31*x41    
     5  +coeff( 41)*x11*x22*x32        
     6  +coeff( 42)*x11*x24*x31        
     7  +coeff( 43)*x12*x23    *x42    
     8  +coeff( 44)*x12*x23*x31*x41    
      x_r12p5_q3en=x_r12p5_q3en
     9  +coeff( 45)            *x41*x51
     1  +coeff( 46)    *x23            
     2  +coeff( 47)    *x21        *x52
     3  +coeff( 48)*x11*x21    *x41    
     4  +coeff( 49)*x11            *x52
     5  +coeff( 50)    *x24            
     6  +coeff( 51)        *x32    *x51
     7  +coeff( 52)    *x22*x31*x41    
     8  +coeff( 53)*x12*x22            
      x_r12p5_q3en=x_r12p5_q3en
     9  +coeff( 54)    *x22*x32        
     1  +coeff( 55)*x12        *x42    
     2  +coeff( 56)    *x21*x32    *x51
     3  +coeff( 57)*x11*x21    *x41*x52
     4  +coeff( 58)*x13*x21            
     5  +coeff( 59)    *x25        *x51
     6  +coeff( 60)*x11    *x32*x42    
     7  +coeff( 61)*x11    *x31*x45    
     8  +coeff( 62)    *x25*x31*x41*x51
      x_r12p5_q3en=x_r12p5_q3en
     9  +coeff( 63)*x13*x24            
     1  +coeff( 64)*x12*x23*x32        
     2  +coeff( 65)*x11        *x41    
     3  +coeff( 66)    *x22*x31        
     4  +coeff( 67)*x11*x21        *x51
     5  +coeff( 68)        *x31*x41*x51
     6  +coeff( 69)    *x22    *x41*x51
     7  +coeff( 70)*x12*x21            
     8  +coeff( 71)        *x31*x43    
      x_r12p5_q3en=x_r12p5_q3en
     9  +coeff( 72)*x12    *x31        
     1  +coeff( 73)    *x23    *x41*x51
     2  +coeff( 74)    *x22    *x42*x51
     3  +coeff( 75)    *x22    *x41*x52
c
      return
      end
      function y_r12p5_q3en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 59)
      data ncoeff/ 58/
      data avdat/  0.1030970E-02/
      data xmin/
     1 -0.99237E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99403E-01, 0.52704E-01, 0.82118E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.56848293E-02, 0.48803440E-02, 0.11102650E+00,-0.52868240E-02,
     + -0.15742473E+00,-0.71636267E-01, 0.56482099E-01,-0.32282203E-01,
     +  0.45961294E-01, 0.39601073E-01, 0.19128181E-01,-0.37938297E-01,
     +  0.31992879E-01,-0.28465018E-02,-0.22216143E-01,-0.96563826E-03,
     + -0.12024452E-01, 0.29304312E-02,-0.69693844E-02,-0.42648655E-02,
     + -0.14873259E-01, 0.13580149E-01,-0.70182416E-02,-0.24124510E-02,
     + -0.63191042E-02,-0.35862881E-02,-0.37585872E-02, 0.91196789E-03,
     + -0.13376424E-01,-0.45149289E-02, 0.28007885E-03, 0.17899830E-02,
     +  0.52635022E-03, 0.40329926E-03, 0.30672078E-02, 0.14398434E-02,
     + -0.96083881E-03,-0.21668584E-02, 0.19414718E-02,-0.10483962E-01,
     + -0.14186932E-01,-0.50267335E-02,-0.58469167E-02,-0.74541820E-02,
     +  0.10542391E-02,-0.52991421E-02, 0.16752128E-02, 0.38776184E-02,
     +  0.27744658E-02,-0.13457057E-02, 0.62922067E-02, 0.11250779E-02,
     +  0.40934468E-02, 0.36423103E-03,-0.21013147E-02,-0.41820714E-02,
     + -0.53733150E-02, 0.55883424E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x25 = x24*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_r12p5_q3en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)        *x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x21*x31        
      y_r12p5_q3en=y_r12p5_q3en
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)*x11    *x31        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)*x11*x21    *x41    
     5  +coeff( 14)*x11                
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)*x11        *x41*x51
     8  +coeff( 17)*x12        *x41    
      y_r12p5_q3en=y_r12p5_q3en
     9  +coeff( 18)    *x22            
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)        *x31*x42    
     4  +coeff( 22)*x11*x21*x31        
     5  +coeff( 23)*x12    *x31        
     6  +coeff( 24)*x11*x21            
     7  +coeff( 25)            *x43    
     8  +coeff( 26)    *x21*x31    *x51
      y_r12p5_q3en=y_r12p5_q3en
     9  +coeff( 27)        *x31    *x52
     1  +coeff( 28)*x12                
     2  +coeff( 29)        *x32*x41    
     3  +coeff( 30)        *x33        
     4  +coeff( 31)        *x33*x41    
     5  +coeff( 32)            *x42    
     6  +coeff( 33)    *x21        *x51
     7  +coeff( 34)                *x52
     8  +coeff( 35)        *x31*x41    
      y_r12p5_q3en=y_r12p5_q3en
     9  +coeff( 36)        *x32        
     1  +coeff( 37)*x11    *x31    *x51
     2  +coeff( 38)    *x22    *x41*x51
     3  +coeff( 39)*x12        *x41*x51
     4  +coeff( 40)*x12*x22    *x41    
     5  +coeff( 41)*x12*x22*x31        
     6  +coeff( 42)*x12    *x31*x42    
     7  +coeff( 43)*x13*x21    *x41    
     8  +coeff( 44)*x13*x21*x31        
      y_r12p5_q3en=y_r12p5_q3en
     9  +coeff( 45)    *x23*x31        
     1  +coeff( 46)*x11*x22    *x41    
     2  +coeff( 47)*x12*x22            
     3  +coeff( 48)*x12*x21    *x41    
     4  +coeff( 49)*x11*x21    *x43    
     5  +coeff( 50)*x11        *x43*x51
     6  +coeff( 51)*x11*x21*x31*x42    
     7  +coeff( 52)*x13*x21            
     8  +coeff( 53)    *x25    *x41    
      y_r12p5_q3en=y_r12p5_q3en
     9  +coeff( 54)*x13            *x51
     1  +coeff( 55)    *x24*x31    *x51
     2  +coeff( 56)*x11*x21    *x43*x51
     3  +coeff( 57)*x12    *x32*x41    
     4  +coeff( 58)    *x23    *x43*x52
c
      return
      end
      function t_r12p5_q3en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/ -0.7202234E-02/
      data xmin/
     1 -0.99237E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99403E-01, 0.52704E-01, 0.82118E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.33030999E-02, 0.45189064E-01,-0.86136334E-01, 0.22723956E-01,
     +  0.99971592E-02,-0.67343456E-02,-0.42432514E-02,-0.35670467E-02,
     + -0.33038240E-02,-0.14338724E-02,-0.12094959E-02, 0.39188094E-02,
     +  0.12239158E-02, 0.86482259E-03,-0.12399523E-02, 0.34162272E-02,
     +  0.81114878E-03, 0.55076013E-03, 0.86442759E-03, 0.17053910E-02,
     +  0.12276350E-02, 0.11018312E-02,-0.41466119E-03, 0.47683442E-03,
     +  0.48906083E-03,-0.35184729E-03,-0.85551714E-04,-0.13174127E-03,
     + -0.40436993E-03,-0.33071998E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_r12p5_q3en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11*x21            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)*x12                
      t_r12p5_q3en=t_r12p5_q3en
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)*x11        *x42    
     5  +coeff( 14)*x12            *x51
     6  +coeff( 15)*x11*x21        *x51
     7  +coeff( 16)*x13*x22            
     8  +coeff( 17)        *x32        
      t_r12p5_q3en=t_r12p5_q3en
     9  +coeff( 18)            *x42    
     1  +coeff( 19)*x11    *x32        
     2  +coeff( 20)    *x21*x32        
     3  +coeff( 21)*x11    *x31*x41    
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)        *x32    *x51
     6  +coeff( 24)            *x42*x51
     7  +coeff( 25)    *x21        *x52
     8  +coeff( 26)*x13        *x41    
      t_r12p5_q3en=t_r12p5_q3en
     9  +coeff( 27)        *x31        
     1  +coeff( 28)*x11    *x31        
     2  +coeff( 29)    *x21*x31        
     3  +coeff( 30)    *x21    *x41    
c
      return
      end
      function p_r12p5_q3en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 38)
      data ncoeff/ 37/
      data avdat/  0.5383748E-03/
      data xmin/
     1 -0.99237E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99403E-01, 0.52704E-01, 0.82118E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.18320525E-02, 0.11072148E-02,-0.39036389E-01, 0.10271245E-01,
     + -0.90891519E-03,-0.92908731E-02,-0.13082425E-01, 0.68952083E-02,
     +  0.11849507E-01, 0.63167224E-02, 0.66622687E-02,-0.96292133E-02,
     +  0.95116580E-02,-0.34639654E-02,-0.77105964E-04,-0.41620512E-02,
     + -0.67287061E-03, 0.62724808E-03,-0.27509402E-02,-0.12318155E-02,
     + -0.14980526E-02, 0.16284779E-02,-0.11598128E-02,-0.13952821E-02,
     +  0.80327044E-03,-0.44518220E-03, 0.11915395E-03, 0.23903222E-03,
     + -0.30790155E-02,-0.32908918E-03,-0.45263529E-03,-0.11270981E-02,
     + -0.17622560E-02, 0.43055337E-03, 0.25453794E-03,-0.44049384E-03,
     +  0.20954228E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x13 = x12*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_r12p5_q3en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31    *x51
      p_r12p5_q3en=p_r12p5_q3en
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11    *x31        
     2  +coeff( 11)*x11        *x41    
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)*x11*x21    *x41    
     5  +coeff( 14)*x12        *x41    
     6  +coeff( 15)*x13                
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)*x11                
      p_r12p5_q3en=p_r12p5_q3en
     9  +coeff( 18)    *x22            
     1  +coeff( 19)        *x31*x42    
     2  +coeff( 20)            *x43    
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)*x11*x21*x31        
     5  +coeff( 23)*x11    *x31    *x51
     6  +coeff( 24)*x12    *x31        
     7  +coeff( 25)        *x34*x41    
     8  +coeff( 26)*x11*x21            
      p_r12p5_q3en=p_r12p5_q3en
     9  +coeff( 27)*x11            *x51
     1  +coeff( 28)*x12                
     2  +coeff( 29)        *x32*x41    
     3  +coeff( 30)            *x41*x52
     4  +coeff( 31)*x11        *x41*x51
     5  +coeff( 32)    *x22    *x41*x51
     6  +coeff( 33)*x13*x21    *x41    
     7  +coeff( 34)        *x31*x41    
     8  +coeff( 35)            *x42    
      p_r12p5_q3en=p_r12p5_q3en
     9  +coeff( 36)        *x33        
     1  +coeff( 37)*x12*x21            
c
      return
      end
      function l_r12p5_q3en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/ -0.2522556E-01/
      data xmin/
     1 -0.99237E-01,-0.48813E-01,-0.94327E-01,-0.31909E-01,-0.49954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.99403E-01, 0.52704E-01, 0.82118E-01, 0.29584E-01, 0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11698964E-01,-0.35788873E+00, 0.10696790E-02, 0.86182554E-03,
     +  0.21468821E+00, 0.41546449E-01,-0.20449130E-01,-0.85159100E-03,
     +  0.11075090E-02, 0.67606213E-03, 0.10592469E-03, 0.84790104E-03,
     + -0.32524016E-01,-0.36181830E-01,-0.11794481E-01,-0.14392188E-01,
     + -0.57555330E-02,-0.69421586E-02,-0.45453347E-02, 0.87738130E-02,
     +  0.86710183E-02, 0.71963412E-02, 0.36703737E-02, 0.82960669E-02,
     + -0.15179124E-02, 0.62579727E-02, 0.10088680E-01, 0.35139823E-02,
     +  0.26686650E-02,-0.45947954E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_r12p5_q3en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)*x11*x21            
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)    *x22        *x51
      l_r12p5_q3en=l_r12p5_q3en
     9  +coeff(  9)        *x32    *x51
     1  +coeff( 10)            *x42*x51
     2  +coeff( 11)                *x53
     3  +coeff( 12)    *x24            
     4  +coeff( 13)                *x51
     5  +coeff( 14)    *x22            
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)*x12                
     8  +coeff( 17)        *x32        
      l_r12p5_q3en=l_r12p5_q3en
     9  +coeff( 18)        *x31*x41    
     1  +coeff( 19)            *x42    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)*x12            *x51
     6  +coeff( 24)    *x23*x31*x41    
     7  +coeff( 25)    *x21    *x41    
     8  +coeff( 26)    *x21*x32        
      l_r12p5_q3en=l_r12p5_q3en
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)*x11    *x32        
     2  +coeff( 29)*x11        *x42    
     3  +coeff( 30)*x11*x21        *x51
c
      return
      end
