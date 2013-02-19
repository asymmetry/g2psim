C Right HRS with septum forward functions. 
C Generated with larger x0 range for Vince (-1.6cm < x0 < 1.6 cm)
C Usage is the same as before.
C                                                   -JJL 2/5/07
      function x_sr6_largex0_ep3    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.1169751E+00/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.31198049E-02, 0.23674698E-09, 0.19829839E-01, 0.32464854E-01,
     +  0.48946647E-04, 0.74692857E-05,-0.46703308E-05, 0.98696737E-05,
     +  0.38373496E-05, 0.15878927E-05, 0.46072148E-06,-0.20652062E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      x_sr6_largex0_ep3    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22    *x41    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)            *x43    
      x_sr6_largex0_ep3    =x_sr6_largex0_ep3    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)                *x51
     2  +coeff( 11)        *x32        
     3  +coeff( 12)        *x31*x41    
c
      return
      end
      function t_sr6_largex0_ep3    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1074092E+00/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17243368E-02, 0.15852484E-04, 0.29980291E-01, 0.39095743E-03,
     +  0.11598411E-03, 0.47661997E-04, 0.85155101E-04, 0.53717788E-04,
     + -0.32217858E-05, 0.24799905E-04, 0.13399607E-04, 0.32448512E-04,
     +  0.20328736E-04,-0.28580724E-04, 0.10385470E-04, 0.84884632E-04,
     + -0.82379593E-04, 0.89814384E-05,-0.76801871E-05, 0.64626256E-05,
     + -0.35156005E-04,-0.88891338E-05, 0.21560681E-04,-0.13289758E-04,
     + -0.34087400E-05,-0.21954406E-05,-0.20143366E-05, 0.10336796E-04,
     +  0.11577779E-04,-0.56447598E-05,-0.49574428E-05,-0.25074331E-04,
     + -0.16139527E-04,-0.58609153E-05, 0.18508437E-04,-0.19148700E-04,
     + -0.13409256E-05,-0.27866752E-04,-0.72788857E-05, 0.44496865E-04,
     + -0.44150162E-04,-0.54428201E-05, 0.12519700E-04,-0.13098019E-04,
     +  0.22413442E-04, 0.21170106E-04, 0.45417670E-04,-0.44981825E-04,
     +  0.37471294E-04, 0.32917527E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sr6_largex0_ep3    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)*x11*x21            
     6  +coeff(  6)                *x51
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_sr6_largex0_ep3    =t_sr6_largex0_ep3    
     9  +coeff(  9)*x12    *x32        
     1  +coeff( 10)        *x32        
     2  +coeff( 11)*x12                
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)*x11*x21    *x41    
     5  +coeff( 14)        *x31*x43    
     6  +coeff( 15)    *x22*x33        
     7  +coeff( 16)    *x22*x31*x42    
     8  +coeff( 17)    *x22*x33*x42    
      t_sr6_largex0_ep3    =t_sr6_largex0_ep3    
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)*x11*x21*x31        
     3  +coeff( 21)        *x32*x42    
     4  +coeff( 22)*x11*x22*x31        
     5  +coeff( 23)*x11*x21*x31*x42    
     6  +coeff( 24)    *x23*x32*x41*x51
     7  +coeff( 25)    *x21*x31        
     8  +coeff( 26)            *x41*x51
      t_sr6_largex0_ep3    =t_sr6_largex0_ep3    
     9  +coeff( 27)                *x52
     1  +coeff( 28)    *x22*x31        
     2  +coeff( 29)        *x32*x41    
     3  +coeff( 30)*x11*x21        *x51
     4  +coeff( 31)    *x21*x32*x41    
     5  +coeff( 32)        *x33*x41    
     6  +coeff( 33)    *x23        *x51
     7  +coeff( 34)*x12    *x31*x41    
     8  +coeff( 35)    *x22    *x43    
      t_sr6_largex0_ep3    =t_sr6_largex0_ep3    
     9  +coeff( 36)        *x32*x43    
     1  +coeff( 37)*x12*x21        *x51
     2  +coeff( 38)*x11*x21*x31*x41*x51
     3  +coeff( 39)*x12*x23            
     4  +coeff( 40)    *x22*x33*x41    
     5  +coeff( 41)    *x22*x31*x43    
     6  +coeff( 42)    *x22    *x42*x52
     7  +coeff( 43)    *x23        *x53
     8  +coeff( 44)*x11*x22    *x41*x52
      t_sr6_largex0_ep3    =t_sr6_largex0_ep3    
     9  +coeff( 45)*x12*x23        *x51
     1  +coeff( 46)*x11*x22    *x43*x51
     2  +coeff( 47)*x11*x21*x31*x41*x53
     3  +coeff( 48)*x11*x23*x33*x41*x51
     4  +coeff( 49)*x11*x22*x31*x42*x53
     5  +coeff( 50)*x12*x23*x32*x42    
c
      return
      end
      function y_sr6_largex0_ep3    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.1542783E-02/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14869246E-02, 0.94804221E-07, 0.18628216E-06, 0.15387720E-06,
     + -0.59780043E-01,-0.15989149E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x51 = x5
c
c                  function
c
      y_sr6_largex0_ep3    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)*x11                
c
      return
      end
      function p_sr6_largex0_ep3    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.1035878E-02/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.98478433E-03,-0.54822687E-01, 0.20456669E-03, 0.11195314E-03,
     +  0.56313766E-04, 0.10221541E-04, 0.30648680E-04,-0.69342236E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x51 = x5
c
c                  function
c
      p_sr6_largex0_ep3    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)*x11        *x41    
     6  +coeff(  6)*x11                
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x21    *x41*x51
c
      return
      end
      function l_sr6_largex0_ep3    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 22)
      data ncoeff/ 21/
      data avdat/  0.1089585E+01/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.69440302E-03,-0.11859817E-03,-0.19412181E-03, 0.39063341E-07,
     +  0.16409582E-02, 0.25787122E-06, 0.49027707E-03, 0.32285400E-05,
     + -0.22742684E-05, 0.27927772E-07, 0.10957646E-05, 0.13551188E-05,
     + -0.27757376E-05, 0.14795171E-05, 0.13108171E-05,-0.22270613E-05,
     +  0.35086721E-05,-0.27007638E-05,-0.22550630E-05,-0.22263489E-05,
     +  0.21174064E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sr6_largex0_ep3    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21            
      l_sr6_largex0_ep3    =l_sr6_largex0_ep3    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)*x12            *x51
     3  +coeff( 12)        *x33    *x51
     4  +coeff( 13)    *x21*x31*x41*x51
     5  +coeff( 14)    *x21    *x41*x52
     6  +coeff( 15)*x11*x22    *x41    
     7  +coeff( 16)*x11*x21        *x52
     8  +coeff( 17)    *x22*x32*x41    
      l_sr6_largex0_ep3    =l_sr6_largex0_ep3    
     9  +coeff( 18)    *x22*x32    *x51
     1  +coeff( 19)    *x22*x33    *x51
     2  +coeff( 20)    *x23    *x42*x51
     3  +coeff( 21)*x11*x23*x32        
c
      return
      end
      function x_sr6_largex0_ep4    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1421731E+00/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.36980279E-02, 0.18263549E-07, 0.19766670E-01, 0.38876370E-01,
     +  0.86372012E-04, 0.10378096E-03, 0.55794306E-04, 0.10454930E-03,
     +  0.40355819E-04,-0.55745199E-05,-0.15500281E-03,-0.99214390E-06,
     +  0.95758602E-04, 0.33072109E-04,-0.12171529E-03, 0.11224384E-04,
     + -0.10734770E-04,-0.59866907E-04,-0.66062938E-04, 0.22419517E-04,
     + -0.36398958E-05, 0.33414868E-04,-0.36780701E-04,-0.14672910E-04,
     + -0.13023053E-03,-0.92016868E-04, 0.14856905E-04,-0.76106908E-05,
     +  0.92093387E-05, 0.19906838E-05,-0.32548933E-05, 0.11587416E-04,
     + -0.73594478E-04,-0.87692088E-05,-0.26784648E-04,-0.12871176E-04,
     + -0.62297049E-05,-0.73371608E-04,-0.15778383E-03,-0.13183481E-04,
     + -0.12790365E-03, 0.48355341E-04,-0.27758990E-04,-0.52124356E-04,
     + -0.48364280E-04,-0.17238026E-05, 0.35355094E-05, 0.31987149E-05,
     + -0.20010302E-05,-0.89574314E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x_sr6_largex0_ep4    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)            *x42    
      x_sr6_largex0_ep4    =x_sr6_largex0_ep4    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)        *x33*x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x22*x33*x41    
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x22*x31*x41    
     7  +coeff( 16)        *x31*x43    
     8  +coeff( 17)*x11*x21    *x41    
      x_sr6_largex0_ep4    =x_sr6_largex0_ep4    
     9  +coeff( 18)*x11*x21*x31*x41    
     1  +coeff( 19)*x11*x21    *x42    
     2  +coeff( 20)        *x32        
     3  +coeff( 21)                *x52
     4  +coeff( 22)        *x31*x42    
     5  +coeff( 23)    *x22*x32        
     6  +coeff( 24)*x11*x21*x32        
     7  +coeff( 25)    *x22*x31*x43    
     8  +coeff( 26)*x11*x21*x31*x43    
      x_sr6_largex0_ep4    =x_sr6_largex0_ep4    
     9  +coeff( 27)        *x32*x41    
     1  +coeff( 28)    *x22        *x51
     2  +coeff( 29)    *x21*x31*x42    
     3  +coeff( 30)    *x22*x31    *x51
     4  +coeff( 31)    *x22    *x41*x51
     5  +coeff( 32)    *x23*x31    *x51
     6  +coeff( 33)    *x22*x32*x42    
     7  +coeff( 34)*x11*x21*x31    *x51
     8  +coeff( 35)*x11*x22*x31*x41    
      x_sr6_largex0_ep4    =x_sr6_largex0_ep4    
     9  +coeff( 36)*x11*x22    *x42    
     1  +coeff( 37)*x11*x23    *x42    
     2  +coeff( 38)*x11*x21*x32*x42    
     3  +coeff( 39)    *x22*x33*x42*x51
     4  +coeff( 40)*x11*x21    *x43*x51
     5  +coeff( 41)    *x22*x32*x43*x51
     6  +coeff( 42)*x11*x23    *x43    
     7  +coeff( 43)    *x23*x33    *x53
     8  +coeff( 44)*x11*x23*x32*x43    
      x_sr6_largex0_ep4    =x_sr6_largex0_ep4    
     9  +coeff( 45)*x11*x22*x32*x42*x52
     1  +coeff( 46)    *x21*x31        
     2  +coeff( 47)    *x21    *x41    
     3  +coeff( 48)        *x33        
     4  +coeff( 49)    *x21        *x52
     5  +coeff( 50)    *x21*x32*x41    
c
      return
      end
      function t_sr6_largex0_ep4    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1273352E+00/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21574972E-02, 0.18624864E-03, 0.30307142E-01, 0.10232345E-02,
     + -0.21192634E-02,-0.11409142E-02,-0.35251720E-04,-0.18794143E-02,
     +  0.74871161E-04,-0.95008966E-03,-0.11832546E-02,-0.53703152E-04,
     + -0.32894398E-03,-0.54996635E-03,-0.36235413E-03,-0.63322892E-04,
     + -0.75578061E-03, 0.78752183E-03, 0.19100196E-03, 0.25693010E-03,
     +  0.18594328E-03, 0.50976141E-05,-0.38932871E-04,-0.20076247E-03,
     + -0.16608468E-04, 0.52662599E-04,-0.18043521E-03,-0.36755049E-04,
     + -0.76878496E-04, 0.21533910E-03,-0.61195151E-05,-0.43471166E-03,
     + -0.28193573E-03, 0.23894091E-03,-0.61420182E-03,-0.91039961E-04,
     +  0.14296493E-03, 0.66787080E-03,-0.46260390E-03, 0.41608313E-04,
     +  0.21320036E-04,-0.33439268E-04,-0.61461666E-04, 0.27483185E-04,
     + -0.55552719E-04,-0.48379490E-04,-0.59037324E-04,-0.10052563E-03,
     + -0.78780540E-04, 0.43638942E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sr6_largex0_ep4    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22    *x42    
     6  +coeff(  6)    *x22*x31*x43    
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x22*x31*x41    
      t_sr6_largex0_ep4    =t_sr6_largex0_ep4    
     9  +coeff(  9)        *x32        
     1  +coeff( 10)*x11*x21*x31*x41    
     2  +coeff( 11)*x11*x21    *x42    
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x22*x32        
     6  +coeff( 15)*x11*x21*x32        
     7  +coeff( 16)    *x22*x31    *x52
     8  +coeff( 17)    *x22*x32*x42    
      t_sr6_largex0_ep4    =t_sr6_largex0_ep4    
     9  +coeff( 18)*x12*x22    *x43    
     1  +coeff( 19)    *x22            
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)            *x42    
     4  +coeff( 22)            *x41*x51
     5  +coeff( 23)*x12                
     6  +coeff( 24)    *x22*x31        
     7  +coeff( 25)    *x21    *x42    
     8  +coeff( 26)        *x31*x42    
      t_sr6_largex0_ep4    =t_sr6_largex0_ep4    
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x22*x31    *x51
     2  +coeff( 29)*x11*x21*x31    *x51
     3  +coeff( 30)    *x22*x31*x42    
     4  +coeff( 31)*x11*x23*x31        
     5  +coeff( 32)*x11*x22*x31*x41    
     6  +coeff( 33)*x12*x22    *x41    
     7  +coeff( 34)    *x23*x31*x42    
     8  +coeff( 35)*x11*x21*x31*x43    
      t_sr6_largex0_ep4    =t_sr6_largex0_ep4    
     9  +coeff( 36)*x11*x21    *x43*x51
     1  +coeff( 37)*x12*x23    *x41    
     2  +coeff( 38)*x11*x22*x31*x43    
     3  +coeff( 39)*x12*x22*x31*x43    
     4  +coeff( 40)        *x32*x41    
     5  +coeff( 41)    *x21*x31    *x51
     6  +coeff( 42)            *x41*x52
     7  +coeff( 43)*x11*x21*x31        
     8  +coeff( 44)*x12    *x31        
      t_sr6_largex0_ep4    =t_sr6_largex0_ep4    
     9  +coeff( 45)    *x23*x31        
     1  +coeff( 46)    *x22    *x41*x51
     2  +coeff( 47)    *x21*x31*x41*x51
     3  +coeff( 48)        *x31*x42*x51
     4  +coeff( 49)            *x43*x51
     5  +coeff( 50)    *x21*x31    *x52
c
      return
      end
      function y_sr6_largex0_ep4    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.1557288E-02/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14906939E-02, 0.33981115E-06,-0.83851552E-06,-0.49880862E-06,
     + -0.71710788E-01,-0.15964152E-01, 0.20921729E-03, 0.10729938E-03,
     +  0.57566940E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x51 = x5
c
c                  function
c
      y_sr6_largex0_ep4    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21*x31        
      y_sr6_largex0_ep4    =y_sr6_largex0_ep4    
     9  +coeff(  9)*x11        *x41    
c
      return
      end
      function p_sr6_largex0_ep4    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.8546625E-03/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.80439553E-03,-0.53985097E-01, 0.10516364E-02, 0.47748338E-03,
     +  0.21268705E-03, 0.28129047E-03, 0.14523270E-03, 0.63457264E-03,
     +  0.44449025E-05,-0.42945463E-04, 0.15946261E-03, 0.56922121E-03,
     +  0.12279488E-02, 0.78378065E-03, 0.73043506E-04,-0.18213115E-04,
     + -0.32616023E-04, 0.15474590E-03,-0.24460498E-05, 0.57257671E-03,
     +  0.23282127E-03, 0.13244838E-04,-0.24369876E-04, 0.18018311E-04,
     + -0.25610001E-04, 0.45183300E-04, 0.16214492E-03,-0.13384552E-04,
     + -0.51182310E-05, 0.98862962E-04, 0.24994844E-03, 0.18955069E-03,
     + -0.10263336E-04,-0.24177967E-04, 0.77916302E-05, 0.56138575E-04,
     +  0.55743835E-05,-0.42341471E-05,-0.66879111E-05, 0.13707860E-03,
     + -0.85604825E-05, 0.34154651E-04, 0.13791518E-03, 0.61565870E-05,
     + -0.31854470E-04,-0.22255603E-04, 0.22967748E-04,-0.24519086E-04,
     + -0.77472352E-04, 0.62020677E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      p_sr6_largex0_ep4    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)*x11                
     6  +coeff(  6)*x11        *x41    
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x21    *x42    
      p_sr6_largex0_ep4    =p_sr6_largex0_ep4    
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21*x31*x42    
     5  +coeff( 14)    *x21    *x43    
     6  +coeff( 15)    *x23*x32*x41    
     7  +coeff( 16)    *x23*x31*x42    
     8  +coeff( 17)    *x21    *x41*x51
      p_sr6_largex0_ep4    =p_sr6_largex0_ep4    
     9  +coeff( 18)*x11        *x42    
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x21*x32*x41    
     3  +coeff( 21)    *x21*x31*x43    
     4  +coeff( 22)*x11*x22*x31*x41    
     5  +coeff( 23)*x11    *x33*x41    
     6  +coeff( 24)    *x23*x33        
     7  +coeff( 25)    *x21*x31    *x51
     8  +coeff( 26)*x11    *x32        
      p_sr6_largex0_ep4    =p_sr6_largex0_ep4    
     9  +coeff( 27)*x11    *x31*x41    
     1  +coeff( 28)*x11        *x41*x51
     2  +coeff( 29)    *x23*x31        
     3  +coeff( 30)    *x21*x33        
     4  +coeff( 31)*x11    *x31*x42    
     5  +coeff( 32)*x11        *x43    
     6  +coeff( 33)*x11    *x32    *x51
     7  +coeff( 34)*x11    *x32*x42    
     8  +coeff( 35)*x11    *x31*x43    
      p_sr6_largex0_ep4    =p_sr6_largex0_ep4    
     9  +coeff( 36)*x11    *x33*x42    
     1  +coeff( 37)    *x22    *x42    
     2  +coeff( 38)    *x22*x31    *x51
     3  +coeff( 39)*x11*x23            
     4  +coeff( 40)*x11    *x32*x41    
     5  +coeff( 41)*x11            *x53
     6  +coeff( 42)    *x21*x33*x41    
     7  +coeff( 43)    *x21*x32*x42    
     8  +coeff( 44)        *x32*x43    
      p_sr6_largex0_ep4    =p_sr6_largex0_ep4    
     9  +coeff( 45)    *x21    *x41*x53
     1  +coeff( 46)*x11*x23    *x41    
     2  +coeff( 47)*x11*x22    *x42    
     3  +coeff( 48)*x11*x21*x31*x42    
     4  +coeff( 49)    *x23    *x43    
     5  +coeff( 50)    *x21*x32*x43    
c
      return
      end
      function l_sr6_largex0_ep4    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 19)
      data ncoeff/ 18/
      data avdat/  0.1308936E+01/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.69060444E-03, 0.64195675E-03, 0.12116693E-02, 0.24868788E-07,
     +  0.19676243E-02, 0.26332209E-05, 0.59048209E-03, 0.35573594E-05,
     +  0.21969627E-05,-0.10094931E-04, 0.23855105E-06, 0.24370295E-05,
     + -0.77639197E-05,-0.33109861E-05,-0.44125509E-05,-0.15901387E-05,
     + -0.23021396E-05,-0.31955653E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_sr6_largex0_ep4    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21            
      l_sr6_largex0_ep4    =l_sr6_largex0_ep4    
     9  +coeff(  9)                *x51
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x22*x33        
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)*x11*x21    *x41    
     6  +coeff( 15)    *x22*x31*x41    
     7  +coeff( 16)*x11*x23            
     8  +coeff( 17)*x11*x23*x31        
      l_sr6_largex0_ep4    =l_sr6_largex0_ep4    
     9  +coeff( 18)*x12        *x42*x51
c
      return
      end
      function x_sr6_largex0_ep5    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1728655E+00/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.43138121E-02,-0.18752314E-05, 0.19646754E-01, 0.44832237E-01,
     +  0.48295068E-03, 0.78113058E-04,-0.66829694E-03, 0.14809449E-03,
     +  0.18046242E-03,-0.79157302E-03,-0.11847159E-02,-0.42665645E-03,
     +  0.29114068E-04,-0.19498038E-03,-0.28699037E-03,-0.73431797E-05,
     + -0.14306422E-03,-0.91750029E-03,-0.21646103E-04,-0.12734623E-03,
     + -0.12231441E-03,-0.19034503E-03, 0.18425166E-04,-0.95458767E-04,
     + -0.73984714E-03, 0.70548867E-03, 0.52555936E-03, 0.22917096E-04,
     + -0.47814406E-05,-0.20286836E-04,-0.64854816E-04,-0.32191645E-04,
     + -0.33618730E-04, 0.22773320E-04,-0.82445498E-04,-0.75574950E-04,
     + -0.92186223E-04,-0.12183288E-03,-0.10087029E-03, 0.18309476E-03,
     +  0.35294468E-03, 0.15099630E-03,-0.12773564E-03, 0.24394132E-03,
     +  0.77046807E-05, 0.46295627E-05,-0.67823489E-05, 0.42740439E-05,
     + -0.17549512E-04, 0.16210952E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x_sr6_largex0_ep5    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21    *x42    
     8  +coeff(  8)        *x31*x41    
      x_sr6_largex0_ep5    =x_sr6_largex0_ep5    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x22*x31*x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)*x11*x21*x31*x41    
     4  +coeff( 13)        *x32        
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22*x32        
     7  +coeff( 16)    *x22*x33        
     8  +coeff( 17)*x11*x21*x32        
      x_sr6_largex0_ep5    =x_sr6_largex0_ep5    
     9  +coeff( 18)    *x22*x31*x43    
     1  +coeff( 19)                *x52
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x22*x33*x41    
     5  +coeff( 23)*x11*x23*x31        
     6  +coeff( 24)*x11*x21    *x43    
     7  +coeff( 25)*x11*x21*x31*x43    
     8  +coeff( 26)    *x22*x33*x43    
      x_sr6_largex0_ep5    =x_sr6_largex0_ep5    
     9  +coeff( 27)*x11*x21*x33*x43    
     1  +coeff( 28)*x11*x21            
     2  +coeff( 29)    *x21*x32*x41    
     3  +coeff( 30)    *x22*x31    *x51
     4  +coeff( 31)*x11*x21*x31        
     5  +coeff( 32)    *x23    *x42    
     6  +coeff( 33)*x11*x21*x31    *x51
     7  +coeff( 34)*x11*x21    *x41*x51
     8  +coeff( 35)*x11*x23    *x41    
      x_sr6_largex0_ep5    =x_sr6_largex0_ep5    
     9  +coeff( 36)*x11*x22*x31*x41    
     1  +coeff( 37)*x11*x22    *x42    
     2  +coeff( 38)    *x23*x32*x43    
     3  +coeff( 39)*x11*x21    *x43*x51
     4  +coeff( 40)*x11*x21    *x42*x52
     5  +coeff( 41)*x11*x23    *x43    
     6  +coeff( 42)*x11*x22    *x42*x52
     7  +coeff( 43)*x11*x22*x32*x43    
     8  +coeff( 44)*x11*x21*x31*x43*x52
      x_sr6_largex0_ep5    =x_sr6_largex0_ep5    
     9  +coeff( 45)    *x21    *x41    
     1  +coeff( 46)    *x21        *x51
     2  +coeff( 47)            *x41*x51
     3  +coeff( 48)    *x21*x32        
     4  +coeff( 49)    *x21*x31*x41    
     5  +coeff( 50)            *x43    
c
      return
      end
      function t_sr6_largex0_ep5    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1640479E+00/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27406984E-02, 0.48823902E-04, 0.30073244E-01, 0.28409555E-02,
     + -0.49844170E-04, 0.56890240E-04,-0.11883391E-03,-0.56839380E-02,
     + -0.77154087E-02,-0.26428022E-02,-0.15748701E-02,-0.36159924E-02,
     +  0.69372043E-04,-0.28899848E-03,-0.10936316E-02,-0.12527721E-02,
     + -0.19736437E-02,-0.32456589E-02,-0.17249747E-03,-0.49314246E-03,
     + -0.82479342E-03,-0.38034569E-02, 0.15641177E-02,-0.98830918E-04,
     + -0.56044635E-04,-0.11174138E-03,-0.28991919E-04, 0.30726561E-03,
     +  0.56885416E-03,-0.78658119E-03,-0.41311834E-03, 0.20178474E-03,
     +  0.54123165E-03,-0.13752208E-02,-0.10801475E-02, 0.42987624E-03,
     + -0.84452331E-03,-0.19323435E-02, 0.73359674E-03,-0.42275689E-03,
     +  0.36790303E-02, 0.84821862E-03,-0.18184170E-02, 0.95613126E-03,
     +  0.95677991E-04, 0.33537723E-04, 0.52847030E-04,-0.17338867E-03,
     +  0.85519408E-04, 0.82506951E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sr6_largex0_ep5    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x22*x31*x41    
      t_sr6_largex0_ep5    =t_sr6_largex0_ep5    
     9  +coeff(  9)    *x22    *x42    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x22*x32        
     3  +coeff( 12)*x11*x21    *x42    
     4  +coeff( 13)    *x22*x33        
     5  +coeff( 14)*x11*x23*x31*x41    
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)*x11*x21    *x41    
     8  +coeff( 17)*x11*x21*x31*x41    
      t_sr6_largex0_ep5    =t_sr6_largex0_ep5    
     9  +coeff( 18)    *x22*x31*x43    
     1  +coeff( 19)                *x52
     2  +coeff( 20)*x11*x21*x31        
     3  +coeff( 21)*x11*x21*x32        
     4  +coeff( 22)*x11*x21*x31*x43    
     5  +coeff( 23)    *x22*x33*x43    
     6  +coeff( 24)*x12                
     7  +coeff( 25)    *x23            
     8  +coeff( 26)    *x22*x31    *x51
      t_sr6_largex0_ep5    =t_sr6_largex0_ep5    
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)    *x22*x32*x41    
     2  +coeff( 29)    *x22    *x42*x52
     3  +coeff( 30)*x11*x21*x33*x41    
     4  +coeff( 31)*x11*x23*x31    *x51
     5  +coeff( 32)*x11*x23    *x41*x51
     6  +coeff( 33)*x11*x22    *x42*x51
     7  +coeff( 34)*x11*x21*x31*x42*x51
     8  +coeff( 35)*x11*x21    *x43*x51
      t_sr6_largex0_ep5    =t_sr6_largex0_ep5    
     9  +coeff( 36)*x11*x21*x31    *x53
     1  +coeff( 37)*x12        *x43*x51
     2  +coeff( 38)    *x22*x31*x42*x52
     3  +coeff( 39)*x11*x23    *x43    
     4  +coeff( 40)*x12        *x42*x53
     5  +coeff( 41)*x11*x21*x33*x43    
     6  +coeff( 42)*x12*x22*x32*x42    
     7  +coeff( 43)    *x22*x32*x43*x52
     8  +coeff( 44)*x12        *x43*x53
      t_sr6_largex0_ep5    =t_sr6_largex0_ep5    
     9  +coeff( 45)        *x31*x41    
     1  +coeff( 46)*x11        *x41    
     2  +coeff( 47)    *x21*x32        
     3  +coeff( 48)            *x43    
     4  +coeff( 49)    *x21    *x41*x51
     5  +coeff( 50)            *x42*x51
c
      return
      end
      function y_sr6_largex0_ep5    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.1711239E-02/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16336213E-02,-0.20988482E-05,-0.70192492E-07,-0.28552747E-05,
     + -0.83045535E-01,-0.15902454E-01, 0.33239511E-03, 0.16329273E-03,
     +  0.13792668E-03, 0.69145870E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x51 = x5
c
c                  function
c
      y_sr6_largex0_ep5    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21*x31        
      y_sr6_largex0_ep5    =y_sr6_largex0_ep5    
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)*x11    *x31        
c
      return
      end
      function p_sr6_largex0_ep5    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.8225234E-03/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.77178882E-03,-0.53907607E-01, 0.97783946E-03, 0.21979515E-03,
     +  0.42473627E-03, 0.27662065E-03, 0.27730791E-02, 0.13682196E-05,
     +  0.11902922E-03, 0.22982240E-02, 0.17606262E-02, 0.38712798E-02,
     +  0.27822058E-02, 0.37890284E-04,-0.30682670E-04, 0.54204790E-03,
     +  0.58455288E-03, 0.23522873E-03, 0.66788279E-03, 0.26222807E-03,
     + -0.46426267E-04, 0.47738475E-03, 0.40570874E-03, 0.88610314E-03,
     +  0.56594244E-03,-0.19357793E-03, 0.14625962E-04,-0.52041578E-04,
     +  0.13627544E-03,-0.52850999E-04, 0.54177395E-04, 0.82687555E-04,
     +  0.17806405E-03, 0.14908312E-03, 0.18692743E-03, 0.24842017E-03,
     + -0.75492069E-04,-0.33940224E-03, 0.22111174E-03,-0.28936384E-03,
     +  0.10003397E-04, 0.22948920E-04,-0.12785451E-04,-0.14185021E-04,
     + -0.14827869E-04, 0.15930360E-04, 0.17416764E-04,-0.12674653E-04,
     + -0.14266379E-04,-0.59154991E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sr6_largex0_ep5    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)*x11        *x41    
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31*x41    
      p_sr6_largex0_ep5    =p_sr6_largex0_ep5    
     9  +coeff(  9)*x11    *x31        
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)    *x21*x32*x41    
     3  +coeff( 12)    *x21*x31*x42    
     4  +coeff( 13)    *x21    *x43    
     5  +coeff( 14)    *x23*x32        
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)*x11        *x42    
      p_sr6_largex0_ep5    =p_sr6_largex0_ep5    
     9  +coeff( 18)    *x21*x33        
     1  +coeff( 19)    *x21*x31*x43    
     2  +coeff( 20)*x11*x22*x31*x41    
     3  +coeff( 21)    *x21        *x51
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)*x11    *x32*x41    
     6  +coeff( 24)*x11    *x31*x42    
     7  +coeff( 25)*x11        *x43    
     8  +coeff( 26)    *x21*x33*x41    
      p_sr6_largex0_ep5    =p_sr6_largex0_ep5    
     9  +coeff( 27)*x11*x22*x32        
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)*x11    *x32        
     3  +coeff( 30)    *x23    *x41    
     4  +coeff( 31)*x11    *x33        
     5  +coeff( 32)*x11*x22    *x42    
     6  +coeff( 33)*x11    *x31*x43    
     7  +coeff( 34)*x12*x21    *x41*x51
     8  +coeff( 35)*x12*x21*x31*x42    
      p_sr6_largex0_ep5    =p_sr6_largex0_ep5    
     9  +coeff( 36)*x12*x21    *x42*x51
     1  +coeff( 37)    *x21    *x43*x53
     2  +coeff( 38)*x11*x22*x33*x41    
     3  +coeff( 39)*x12*x21*x31*x43    
     4  +coeff( 40)*x12*x21*x32*x42*x51
     5  +coeff( 41)    *x22*x31        
     6  +coeff( 42)            *x42*x51
     7  +coeff( 43)    *x21        *x52
     8  +coeff( 44)*x11*x21    *x41    
      p_sr6_largex0_ep5    =p_sr6_largex0_ep5    
     9  +coeff( 45)*x11*x21        *x51
     1  +coeff( 46)*x11        *x41*x51
     2  +coeff( 47)*x12*x21            
     3  +coeff( 48)*x12            *x51
     4  +coeff( 49)    *x22*x31    *x51
     5  +coeff( 50)    *x21    *x42*x51
c
      return
      end
      function l_sr6_largex0_ep5    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/  0.1519711E+01/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.58199844E-03, 0.39660167E-05, 0.17055228E-02, 0.35962318E-02,
     +  0.22764623E-02, 0.86766258E-05, 0.69278764E-03, 0.24296234E-04,
     +  0.14546649E-04,-0.58257383E-05,-0.30885993E-04,-0.19748928E-04,
     + -0.15764434E-04,-0.55155066E-04,-0.56571553E-04, 0.19658983E-05,
     + -0.73104202E-05,-0.13405148E-05,-0.72786888E-05,-0.13912942E-04,
     + -0.19148338E-04,-0.18200111E-04,-0.22097887E-04,-0.19140127E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      l_sr6_largex0_ep5    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)                *x51
      l_sr6_largex0_ep5    =l_sr6_largex0_ep5    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)*x11*x21    *x41    
     5  +coeff( 14)    *x22*x31*x41    
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)        *x32        
     8  +coeff( 17)        *x31*x42    
      l_sr6_largex0_ep5    =l_sr6_largex0_ep5    
     9  +coeff( 18)            *x43    
     1  +coeff( 19)*x11*x21*x31        
     2  +coeff( 20)    *x22*x32        
     3  +coeff( 21)*x11*x21*x31*x41    
     4  +coeff( 22)*x11*x21    *x42    
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)    *x22*x33*x42    
c
      return
      end
      function x_sr6_largex0_ep6    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2118092E+00/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.52541564E-02, 0.86472783E-05, 0.19494250E-01, 0.50669599E-01,
     +  0.12547545E-02,-0.34420940E-02, 0.19563068E-03,-0.45038437E-05,
     + -0.25353120E-02, 0.23525577E-04,-0.83416572E-03,-0.72461978E-03,
     +  0.11613978E-03,-0.14293465E-02,-0.23248084E-03, 0.23708175E-03,
     + -0.48350033E-03,-0.50674065E-03,-0.95037278E-03,-0.14202778E-02,
     + -0.78475816E-04,-0.19302359E-03,-0.34273398E-03,-0.11798382E-02,
     +  0.52536203E-03, 0.44953843E-03, 0.76759112E-03, 0.39594557E-04,
     +  0.21261294E-04, 0.90391994E-04,-0.47797069E-03,-0.10096333E-03,
     +  0.82240505E-04,-0.12401282E-02, 0.14674338E-03,-0.86102012E-03,
     +  0.70170303E-04,-0.79485675E-03,-0.33871908E-03,-0.62551594E-03,
     +  0.23195181E-03, 0.45434933E-03,-0.18024495E-04,-0.15982978E-04,
     + -0.23009192E-04, 0.18918159E-04,-0.29710591E-04, 0.11338739E-04,
     + -0.54492593E-04, 0.37866292E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x_sr6_largex0_ep6    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22    *x42    
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)*x11*x21            
      x_sr6_largex0_ep6    =x_sr6_largex0_ep6    
     9  +coeff(  9)    *x22*x31*x41    
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22*x32        
     4  +coeff( 13)    *x22*x33        
     5  +coeff( 14)*x11*x21    *x42    
     6  +coeff( 15)*x11*x23*x31*x41    
     7  +coeff( 16)            *x42    
     8  +coeff( 17)    *x22*x31        
      x_sr6_largex0_ep6    =x_sr6_largex0_ep6    
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)*x11*x21*x31*x41    
     2  +coeff( 20)    *x22*x31*x43    
     3  +coeff( 21)                *x52
     4  +coeff( 22)*x11*x21*x31        
     5  +coeff( 23)*x11*x21*x32        
     6  +coeff( 24)*x11*x21*x31*x43    
     7  +coeff( 25)    *x22*x33*x43    
     8  +coeff( 26)    *x22*x33*x42*x52
      x_sr6_largex0_ep6    =x_sr6_largex0_ep6    
     9  +coeff( 27)*x11*x21*x33*x43    
     1  +coeff( 28)        *x32        
     2  +coeff( 29)            *x43    
     3  +coeff( 30)            *x42*x52
     4  +coeff( 31)    *x22    *x43    
     5  +coeff( 32)*x11*x21*x31    *x51
     6  +coeff( 33)*x11*x21    *x41*x51
     7  +coeff( 34)    *x22*x33*x42    
     8  +coeff( 35)*x11*x21    *x43    
      x_sr6_largex0_ep6    =x_sr6_largex0_ep6    
     9  +coeff( 36)    *x22*x32*x43    
     1  +coeff( 37)*x11*x22        *x52
     2  +coeff( 38)    *x22*x31*x42*x52
     3  +coeff( 39)*x11*x23    *x42    
     4  +coeff( 40)*x11*x21    *x43*x51
     5  +coeff( 41)*x11*x21    *x42*x52
     6  +coeff( 42)*x11*x21    *x43*x53
     7  +coeff( 43)        *x31    *x51
     8  +coeff( 44)            *x41*x51
      x_sr6_largex0_ep6    =x_sr6_largex0_ep6    
     9  +coeff( 45)    *x23            
     1  +coeff( 46)    *x21*x31    *x51
     2  +coeff( 47)            *x41*x52
     3  +coeff( 48)*x11    *x31        
     4  +coeff( 49)    *x23*x31        
     5  +coeff( 50)    *x21*x33        
c
      return
      end
      function t_sr6_largex0_ep6    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2011861E+00/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.34422940E-02,-0.11541248E-03, 0.29518889E-01, 0.46518384E-02,
     + -0.31979199E-03,-0.15393973E-03,-0.47593794E-03,-0.13619431E-04,
     + -0.10298392E-01,-0.14927451E-01,-0.66686678E-02, 0.49190176E-05,
     + -0.18107727E-02,-0.30922648E-02,-0.26330720E-02,-0.52875187E-02,
     + -0.13142278E-03,-0.17000058E-02,-0.85029536E-03,-0.19483654E-02,
     +  0.67989575E-04,-0.15691273E-03, 0.34743169E-03,-0.28680356E-02,
     + -0.14402667E-02, 0.34995282E-02,-0.17824729E-02,-0.72518662E-02,
     + -0.90540657E-02,-0.13249064E-01,-0.17688031E-02, 0.98952362E-02,
     +  0.93615223E-02,-0.83446459E-04,-0.21298244E-03,-0.68665424E-04,
     +  0.17938730E-04,-0.38924412E-03, 0.17774156E-03,-0.42692485E-03,
     +  0.23602550E-02,-0.20942723E-02,-0.13824294E-02,-0.14005214E-02,
     + -0.17931929E-02,-0.95586281E-03, 0.19545767E-02,-0.11261000E-02,
     +  0.84100710E-03,-0.78848971E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sr6_largex0_ep6    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)        *x31*x42    
      t_sr6_largex0_ep6    =t_sr6_largex0_ep6    
     9  +coeff(  9)    *x22*x31*x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x22    *x43    
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22*x32        
     7  +coeff( 16)*x11*x21    *x42    
     8  +coeff( 17)*x11*x21*x31*x42    
      t_sr6_largex0_ep6    =t_sr6_largex0_ep6    
     9  +coeff( 18)*x11*x23*x31*x41    
     1  +coeff( 19)*x11*x21*x31        
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)        *x32*x42    
     4  +coeff( 22)        *x32    *x52
     5  +coeff( 23)*x11*x23            
     6  +coeff( 24)*x11*x21*x31*x41    
     7  +coeff( 25)    *x22*x32*x41    
     8  +coeff( 26)    *x22*x32*x42    
      t_sr6_largex0_ep6    =t_sr6_largex0_ep6    
     9  +coeff( 27)*x11*x23*x32        
     1  +coeff( 28)    *x22*x33*x42    
     2  +coeff( 29)    *x22*x31*x42*x52
     3  +coeff( 30)    *x22*x31*x43*x52
     4  +coeff( 31)*x12*x22    *x42*x52
     5  +coeff( 32)    *x22*x33*x42*x52
     6  +coeff( 33)    *x22*x33*x43*x52
     7  +coeff( 34)            *x41*x51
     8  +coeff( 35)                *x52
      t_sr6_largex0_ep6    =t_sr6_largex0_ep6    
     9  +coeff( 36)    *x23*x31        
     1  +coeff( 37)    *x22    *x41*x51
     2  +coeff( 38)        *x31*x42*x51
     3  +coeff( 39)*x12        *x43    
     4  +coeff( 40)*x12        *x41*x52
     5  +coeff( 41)    *x22*x31*x41*x52
     6  +coeff( 42)*x11*x21*x31*x43    
     7  +coeff( 43)*x12*x22*x31*x41    
     8  +coeff( 44)*x12*x22    *x41*x51
      t_sr6_largex0_ep6    =t_sr6_largex0_ep6    
     9  +coeff( 45)*x12*x23*x31*x41    
     1  +coeff( 46)*x12*x23    *x42    
     2  +coeff( 47)*x12*x22    *x43    
     3  +coeff( 48)*x12*x22    *x42*x51
     4  +coeff( 49)*x12*x22*x31    *x52
     5  +coeff( 50)*x12*x23*x33        
c
      return
      end
      function y_sr6_largex0_ep6    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 19)
      data ncoeff/ 18/
      data avdat/ -0.1978577E-02/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.18891032E-02,-0.63034929E-06,-0.73942812E-06,-0.94514951E-01,
     + -0.15829599E-01, 0.43476952E-03, 0.23212565E-03, 0.21937818E-03,
     +  0.11944734E-02,-0.35508009E-03, 0.10576416E-03, 0.24548711E-03,
     +  0.10765922E-02,-0.10586034E-03, 0.13613647E-02, 0.11035991E-02,
     + -0.53695554E-03, 0.48724774E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
c
c                  function
c
      y_sr6_largex0_ep6    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)*x11        *x41    
      y_sr6_largex0_ep6    =y_sr6_largex0_ep6    
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x21*x33*x41    
     2  +coeff( 11)*x11    *x31        
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x21*x31*x42    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x21*x32*x42    
      y_sr6_largex0_ep6    =y_sr6_largex0_ep6    
     9  +coeff( 18)    *x21*x34*x41    
c
      return
      end
      function p_sr6_largex0_ep6    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.8625893E-03/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.80510409E-03,-0.54280747E-01, 0.42119171E-03, 0.18704680E-03,
     +  0.27353156E-03, 0.55305832E-02, 0.11987788E-03, 0.30570454E-03,
     +  0.40925331E-02, 0.19825615E-02, 0.56006759E-02, 0.27218986E-04,
     +  0.13701305E-03, 0.10146720E-02,-0.41519328E-04,-0.28171414E-03,
     +  0.62427945E-02,-0.13259731E-02,-0.12688863E-02, 0.81616262E-03,
     +  0.27163915E-03, 0.28050202E-03,-0.16368293E-03, 0.10435533E-02,
     + -0.14309058E-03,-0.35242090E-03, 0.22447604E-03,-0.55146727E-04,
     +  0.13519789E-03, 0.16151197E-03, 0.76328107E-03, 0.96426555E-03,
     +  0.57943439E-03,-0.89233901E-04, 0.24213683E-03, 0.85447056E-04,
     + -0.38294902E-03,-0.21518751E-03,-0.14331803E-03,-0.30768212E-03,
     + -0.85091720E-04, 0.58276614E-03, 0.25871865E-03, 0.15074050E-03,
     + -0.36050726E-03, 0.21784424E-03, 0.50320348E-03, 0.23833620E-03,
     +  0.22415184E-03,-0.52631687E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sr6_largex0_ep6    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x23*x31*x41    
     8  +coeff(  8)*x11        *x41    
      p_sr6_largex0_ep6    =p_sr6_largex0_ep6    
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)    *x21*x32*x41    
     2  +coeff( 11)    *x21    *x43    
     3  +coeff( 12)    *x23*x31*x42    
     4  +coeff( 13)*x11    *x31        
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21*x31*x42    
      p_sr6_largex0_ep6    =p_sr6_largex0_ep6    
     9  +coeff( 18)    *x21*x32*x42    
     1  +coeff( 19)    *x23*x33*x41    
     2  +coeff( 20)*x11        *x42    
     3  +coeff( 21)    *x21*x33        
     4  +coeff( 22)*x11*x22*x31*x41    
     5  +coeff( 23)*x11    *x33*x41    
     6  +coeff( 24)    *x23*x31*x43    
     7  +coeff( 25)*x11*x22*x33*x41    
     8  +coeff( 26)*x11*x22*x31*x43    
      p_sr6_largex0_ep6    =p_sr6_largex0_ep6    
     9  +coeff( 27)    *x23*x31*x42*x52
     1  +coeff( 28)    *x21        *x51
     2  +coeff( 29)    *x21    *x41*x51
     3  +coeff( 30)*x11    *x32        
     4  +coeff( 31)*x11    *x31*x41    
     5  +coeff( 32)*x11    *x31*x42    
     6  +coeff( 33)*x11        *x43    
     7  +coeff( 34)    *x23*x32        
     8  +coeff( 35)*x12*x21    *x41    
      p_sr6_largex0_ep6    =p_sr6_largex0_ep6    
     9  +coeff( 36)    *x22*x32    *x51
     1  +coeff( 37)    *x21    *x42*x52
     2  +coeff( 38)    *x21    *x41*x53
     3  +coeff( 39)*x11    *x32*x42    
     4  +coeff( 40)*x11*x21*x31*x41*x51
     5  +coeff( 41)*x11*x21        *x53
     6  +coeff( 42)    *x23*x32*x41    
     7  +coeff( 43)*x12*x21    *x41*x51
     8  +coeff( 44)    *x23    *x42*x51
      p_sr6_largex0_ep6    =p_sr6_largex0_ep6    
     9  +coeff( 45)    *x21    *x43*x52
     1  +coeff( 46)*x11*x22*x32*x41    
     2  +coeff( 47)*x11    *x32*x43    
     3  +coeff( 48)*x12*x21*x31*x42    
     4  +coeff( 49)    *x22*x33*x42    
     5  +coeff( 50)*x12*x21    *x43    
c
      return
      end
      function l_sr6_largex0_ep6    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.1732157E+01/
      data xmin/
     1 -0.15988E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29982E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.44798368E-03, 0.42518072E-05, 0.25596386E-02, 0.58842786E-02,
     +  0.25892728E-02, 0.18996876E-04, 0.79997390E-03, 0.85877342E-04,
     +  0.38260125E-04,-0.10732886E-04,-0.67316840E-04,-0.30549717E-03,
     +  0.40754412E-05,-0.41113861E-04,-0.51966086E-04,-0.24000340E-03,
     + -0.18511681E-03,-0.38438342E-04,-0.15444205E-05,-0.64093772E-04,
     + -0.87932312E-04,-0.91181202E-04,-0.17587129E-05,-0.83344741E-04,
     + -0.18920854E-03,-0.14491061E-05,-0.69242007E-04,-0.27173508E-04,
     + -0.20864729E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      l_sr6_largex0_ep6    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)                *x51
      l_sr6_largex0_ep6    =l_sr6_largex0_ep6    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)        *x32        
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)*x11*x21    *x41    
     7  +coeff( 16)    *x22*x31*x41    
     8  +coeff( 17)    *x22    *x43    
      l_sr6_largex0_ep6    =l_sr6_largex0_ep6    
     9  +coeff( 18)    *x22*x33*x42    
     1  +coeff( 19)        *x31*x42    
     2  +coeff( 20)    *x22*x32        
     3  +coeff( 21)*x11*x21*x31*x41    
     4  +coeff( 22)*x11*x21    *x42    
     5  +coeff( 23)    *x22*x33        
     6  +coeff( 24)    *x22*x32*x41    
     7  +coeff( 25)    *x22*x31*x42    
     8  +coeff( 26)*x11*x23*x31        
      l_sr6_largex0_ep6    =l_sr6_largex0_ep6    
     9  +coeff( 27)    *x22*x31*x43    
     1  +coeff( 28)*x11*x23*x32        
     2  +coeff( 29)*x11*x21*x31        
c
      return
      end
      function x_sr6_largex0_ep7    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2614051E+00/
      data xmin/
     1 -0.15981E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.64687831E-02, 0.19020488E-01, 0.56062292E-01, 0.23962345E-02,
     + -0.30199985E-03,-0.20664898E-03,-0.74015586E-02, 0.96887222E-03,
     +  0.21098026E-03,-0.15613512E-02,-0.53297840E-02, 0.36346511E-04,
     + -0.86436060E-03,-0.13096940E-02,-0.20494501E-02, 0.51592295E-04,
     + -0.69801725E-03, 0.34793810E-03,-0.16489098E-03, 0.65294524E-04,
     +  0.31411113E-04,-0.42517943E-03,-0.95493527E-03,-0.43533271E-03,
     + -0.23012368E-02, 0.23251447E-03,-0.16655215E-02, 0.56796364E-03,
     +  0.87766489E-03, 0.62942825E-04,-0.17372742E-02, 0.63172229E-04,
     +  0.13468529E-04,-0.83422077E-04,-0.23422200E-03,-0.10391907E-02,
     + -0.25146428E-04,-0.67443587E-03,-0.19751866E-02, 0.39892067E-04,
     +  0.18459118E-03,-0.10590086E-02, 0.41267951E-03,-0.97182509E-03,
     + -0.94798172E-03, 0.60208078E-03,-0.15414404E-02, 0.55647036E-03,
     + -0.13099643E-02, 0.36238565E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x_sr6_largex0_ep7    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)*x11*x21            
     7  +coeff(  7)    *x22    *x42    
     8  +coeff(  8)    *x22*x33*x41    
      x_sr6_largex0_ep7    =x_sr6_largex0_ep7    
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x22*x31*x41    
     3  +coeff( 12)        *x32        
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x22*x32        
     6  +coeff( 15)*x11*x21    *x42    
     7  +coeff( 16)*x11*x21*x31*x42    
     8  +coeff( 17)*x11*x23*x31*x41    
      x_sr6_largex0_ep7    =x_sr6_largex0_ep7    
     9  +coeff( 18)            *x42    
     1  +coeff( 19)                *x52
     2  +coeff( 20)            *x43    
     3  +coeff( 21)        *x32*x42    
     4  +coeff( 22)*x11*x21*x31        
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)    *x22*x32*x41    
     7  +coeff( 25)    *x22    *x43    
     8  +coeff( 26)*x11*x23            
      x_sr6_largex0_ep7    =x_sr6_largex0_ep7    
     9  +coeff( 27)*x11*x21*x31*x41    
     1  +coeff( 28)    *x22*x32*x42    
     2  +coeff( 29)    *x22    *x42*x52
     3  +coeff( 30)*x11*x23*x32        
     4  +coeff( 31)*x11*x23    *x42    
     5  +coeff( 32)        *x33        
     6  +coeff( 33)    *x21*x31    *x51
     7  +coeff( 34)    *x22*x31    *x51
     8  +coeff( 35)    *x21    *x42*x51
      x_sr6_largex0_ep7    =x_sr6_largex0_ep7    
     9  +coeff( 36)    *x22*x31*x42    
     1  +coeff( 37)        *x31*x42*x52
     2  +coeff( 38)*x11*x21*x32        
     3  +coeff( 39)    *x22*x31*x43    
     4  +coeff( 40)*x11*x21*x31    *x51
     5  +coeff( 41)*x11*x21    *x41*x51
     6  +coeff( 42)    *x22*x33*x42    
     7  +coeff( 43)    *x23*x31*x42*x51
     8  +coeff( 44)    *x22*x31*x42*x52
      x_sr6_largex0_ep7    =x_sr6_largex0_ep7    
     9  +coeff( 45)*x11*x21*x31*x43    
     1  +coeff( 46)*x11*x22    *x42*x51
     2  +coeff( 47)*x11*x21*x31*x42*x51
     3  +coeff( 48)    *x23*x32*x42*x51
     4  +coeff( 49)*x11*x21    *x43*x51
     5  +coeff( 50)*x11    *x31*x43*x51
c
      return
      end
      function t_sr6_largex0_ep7    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2224681E+00/
      data xmin/
     1 -0.15981E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.42409548E-02, 0.16103491E-03, 0.30288603E-01, 0.56949896E-02,
     + -0.11040537E-02,-0.84830844E-03,-0.23274792E-03,-0.60458463E-02,
     + -0.88234316E-03,-0.11018180E-01,-0.16217349E-01, 0.21701494E-03,
     + -0.12490287E-03,-0.29298742E-02,-0.25631667E-02,-0.45687794E-02,
     + -0.21381194E-02,-0.35136384E-02,-0.54072128E-02,-0.36862455E-03,
     + -0.41321828E-03,-0.17005455E-02, 0.14054247E-02,-0.39106116E-03,
     + -0.89177786E-03,-0.27025132E-02, 0.21414843E-02,-0.37027120E-02,
     + -0.14758473E-02, 0.25712200E-02,-0.11651902E-02, 0.46515912E-02,
     + -0.16392814E-04,-0.23254668E-03,-0.94962015E-03,-0.10945712E-02,
     +  0.58853446E-03,-0.85992750E-03,-0.29375430E-02,-0.42297135E-03,
     + -0.37971709E-03, 0.42226576E-03,-0.21734505E-02, 0.10370667E-02,
     + -0.47431316E-03, 0.82712102E-03,-0.50188694E-03, 0.69264113E-03,
     +  0.34973159E-03, 0.83816034E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sr6_largex0_ep7    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x22    *x41    
      t_sr6_largex0_ep7    =t_sr6_largex0_ep7    
     9  +coeff(  9)        *x31*x42    
     1  +coeff( 10)    *x22*x31*x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x22*x33        
     4  +coeff( 13)        *x32        
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x22*x32        
     7  +coeff( 16)*x11*x21    *x42    
     8  +coeff( 17)*x11*x23    *x41    
      t_sr6_largex0_ep7    =t_sr6_largex0_ep7    
     9  +coeff( 18)*x11*x23*x31*x41    
     1  +coeff( 19)*x11*x23    *x42    
     2  +coeff( 20)*x11*x23*x33        
     3  +coeff( 21)                *x52
     4  +coeff( 22)            *x43    
     5  +coeff( 23)        *x32*x42    
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)    *x21    *x42*x51
     8  +coeff( 26)*x11*x21    *x43    
      t_sr6_largex0_ep7    =t_sr6_largex0_ep7    
     9  +coeff( 27)*x12        *x43    
     1  +coeff( 28)    *x22*x31*x43    
     2  +coeff( 29)    *x22*x31*x42*x51
     3  +coeff( 30)    *x22    *x42*x52
     4  +coeff( 31)*x11*x23    *x41*x51
     5  +coeff( 32)*x11*x23    *x42*x52
     6  +coeff( 33)            *x41*x51
     7  +coeff( 34)*x12                
     8  +coeff( 35)*x11*x21*x31        
      t_sr6_largex0_ep7    =t_sr6_largex0_ep7    
     9  +coeff( 36)*x12        *x41    
     1  +coeff( 37)        *x33*x41    
     2  +coeff( 38)*x11*x21*x32        
     3  +coeff( 39)*x11*x21*x31*x41    
     4  +coeff( 40)*x11*x21*x31    *x51
     5  +coeff( 41)*x11*x21        *x52
     6  +coeff( 42)*x12        *x42    
     7  +coeff( 43)    *x22    *x43    
     8  +coeff( 44)        *x32*x43    
      t_sr6_largex0_ep7    =t_sr6_largex0_ep7    
     9  +coeff( 45)*x12        *x41*x51
     1  +coeff( 46)    *x22*x31*x41*x51
     2  +coeff( 47)    *x21*x32*x41*x51
     3  +coeff( 48)    *x22    *x41*x52
     4  +coeff( 49)    *x21    *x42*x52
     5  +coeff( 50)*x12*x22    *x41    
c
      return
      end
      function y_sr6_largex0_ep7    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.2711254E-02/
      data xmin/
     1 -0.15981E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26049246E-02, 0.16097142E-05, 0.36552360E-05,-0.10703532E+00,
     + -0.15781429E-01,-0.25256304E-03, 0.17739930E-02, 0.20420768E-03,
     + -0.14924679E-02, 0.12318304E-04, 0.31294444E-03, 0.39934614E-03,
     +  0.24446514E-02, 0.38808798E-02, 0.30857129E-02, 0.26441659E-04,
     +  0.63142454E-03,-0.96508506E-04,-0.58111633E-04, 0.14352720E-03,
     + -0.94486924E-04, 0.19904840E-03, 0.16468607E-02, 0.12830319E-02,
     +  0.25168603E-03, 0.55333097E-02, 0.30521522E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      y_sr6_largex0_ep7    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x21*x33        
      y_sr6_largex0_ep7    =y_sr6_largex0_ep7    
     9  +coeff(  9)    *x21*x33*x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)*x11        *x41    
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21*x31*x42    
     6  +coeff( 15)    *x21    *x43    
     7  +coeff( 16)*x11    *x33        
     8  +coeff( 17)    *x21*x32*x42    
      y_sr6_largex0_ep7    =y_sr6_largex0_ep7    
     9  +coeff( 18)    *x21*x34*x41    
     1  +coeff( 19)    *x21        *x51
     2  +coeff( 20)*x11    *x31        
     3  +coeff( 21)    *x23            
     4  +coeff( 22)*x11        *x42    
     5  +coeff( 23)    *x21*x32*x41    
     6  +coeff( 24)    *x21    *x44    
     7  +coeff( 25)*x11    *x33*x41    
     8  +coeff( 26)    *x21*x33*x43    
      y_sr6_largex0_ep7    =y_sr6_largex0_ep7    
     9  +coeff( 27)    *x21*x34*x44    
c
      return
      end
      function p_sr6_largex0_ep7    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1043073E-02/
      data xmin/
     1 -0.15981E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.98960439E-03,-0.53847034E-01, 0.51975838E-03, 0.67449301E-02,
     +  0.24244700E-03, 0.37726178E-03, 0.49217776E-02, 0.64234952E-02,
     + -0.31701487E-03, 0.10316682E-02, 0.58855359E-04, 0.70969472E-02,
     + -0.13781792E-02, 0.60701882E-03, 0.13577382E-03, 0.10058040E-02,
     + -0.43642175E-03, 0.25409306E-02,-0.22063950E-04,-0.82347955E-03,
     +  0.57264344E-04,-0.12479657E-03, 0.17596697E-03,-0.48419653E-04,
     +  0.92120678E-03, 0.27595711E-03, 0.12081411E-02, 0.10309847E-02,
     +  0.71358401E-03, 0.19176953E-03, 0.37636526E-03, 0.58736384E-03,
     + -0.50908554E-03, 0.25794865E-04, 0.21198121E-03,-0.61284642E-04,
     +  0.61751591E-04, 0.54251112E-03, 0.12250523E-03,-0.70204842E-04,
     +  0.23745447E-03,-0.41289791E-03, 0.34299603E-03, 0.26280928E-03,
     + -0.69834793E-03, 0.19546234E-03, 0.16511510E-03,-0.50140353E-03,
     + -0.28098089E-03,-0.16092388E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sr6_largex0_ep7    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21    *x42    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)    *x21    *x43    
      p_sr6_largex0_ep7    =p_sr6_largex0_ep7    
     9  +coeff(  9)    *x23*x31*x42    
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x23*x31        
     3  +coeff( 12)    *x21*x31*x42    
     4  +coeff( 13)    *x21*x32*x42    
     5  +coeff( 14)    *x23*x32*x41    
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)*x11        *x42    
     8  +coeff( 17)    *x23    *x41    
      p_sr6_largex0_ep7    =p_sr6_largex0_ep7    
     9  +coeff( 18)    *x21*x32*x41    
     1  +coeff( 19)*x11*x22*x31        
     2  +coeff( 20)    *x21*x33*x41    
     3  +coeff( 21)*x11*x22*x31*x41    
     4  +coeff( 22)    *x21        *x51
     5  +coeff( 23)*x11    *x31        
     6  +coeff( 24)*x11*x22            
     7  +coeff( 25)*x11    *x31*x41    
     8  +coeff( 26)    *x21*x33        
      p_sr6_largex0_ep7    =p_sr6_largex0_ep7    
     9  +coeff( 27)*x11    *x31*x42    
     1  +coeff( 28)*x11        *x43    
     2  +coeff( 29)    *x21*x31*x43    
     3  +coeff( 30)    *x22    *x42*x51
     4  +coeff( 31)*x12*x21    *x41*x51
     5  +coeff( 32)*x12*x21*x31*x42    
     6  +coeff( 33)*x11*x21*x31*x41*x53
     7  +coeff( 34)    *x21*x31    *x51
     8  +coeff( 35)*x11    *x32        
      p_sr6_largex0_ep7    =p_sr6_largex0_ep7    
     9  +coeff( 36)*x11*x21        *x51
     1  +coeff( 37)            *x43*x51
     2  +coeff( 38)*x11    *x32*x41    
     3  +coeff( 39)    *x21*x31*x42*x51
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)    *x21*x32    *x52
     6  +coeff( 42)    *x21    *x42*x52
     7  +coeff( 43)*x11*x22    *x42    
     8  +coeff( 44)*x12*x21*x31*x41    
      p_sr6_largex0_ep7    =p_sr6_largex0_ep7    
     9  +coeff( 45)    *x21*x32*x43    
     1  +coeff( 46)    *x23    *x42*x51
     2  +coeff( 47)    *x23    *x41*x52
     3  +coeff( 48)    *x21    *x43*x52
     4  +coeff( 49)*x11*x23*x31*x41    
     5  +coeff( 50)*x11*x23    *x42    
c
      return
      end
      function l_sr6_largex0_ep7    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.1962502E+01/
      data xmin/
     1 -0.15981E-01,-0.54809E-01,-0.19951E-01,-0.29937E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15990E-01, 0.54912E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16191983E-03, 0.44523813E-02, 0.11603788E-01, 0.34843475E-03,
     +  0.28905498E-02, 0.26302758E-04, 0.93149714E-03,-0.55248707E-04,
     +  0.14519703E-05, 0.65224602E-04,-0.28696604E-03,-0.79655950E-03,
     + -0.11722678E-02, 0.10014342E-04,-0.51846140E-03,-0.14893457E-03,
     + -0.32705848E-04,-0.17126999E-03,-0.19192709E-03, 0.19489896E-05,
     + -0.30912721E-03,-0.37764202E-03,-0.41724456E-03, 0.12861273E-04,
     + -0.33869522E-03,-0.11061913E-04,-0.38499642E-04,-0.75676660E-04,
     + -0.36724141E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      l_sr6_largex0_ep7    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)*x11*x21            
      l_sr6_largex0_ep7    =l_sr6_largex0_ep7    
     9  +coeff(  9)        *x32        
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22*x31*x41    
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)    *x22*x33        
     6  +coeff( 15)    *x22    *x43    
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)        *x31*x42    
      l_sr6_largex0_ep7    =l_sr6_largex0_ep7    
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)    *x22*x32        
     2  +coeff( 20)        *x31*x43    
     3  +coeff( 21)*x11*x21*x31*x41    
     4  +coeff( 22)*x11*x21    *x42    
     5  +coeff( 23)    *x22*x31*x42    
     6  +coeff( 24)*x11*x23*x31        
     7  +coeff( 25)    *x22*x31*x43    
     8  +coeff( 26)                *x52
      l_sr6_largex0_ep7    =l_sr6_largex0_ep7    
     9  +coeff( 27)        *x32*x41    
     1  +coeff( 28)*x11*x21*x31        
     2  +coeff( 29)        *x31*x42*x51
c
      return
      end
      function x_sr6_largex0_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.3932923E-02/
      data xmin/
     1 -0.15981E-01,-0.54082E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.54577E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.32847964E-02, 0.14496994E+00, 0.83540985E-02, 0.48856222E-03,
     +  0.34609528E-02,-0.56180754E-03,-0.12724224E-01, 0.19765331E-02,
     + -0.94046583E-02,-0.14561698E-01, 0.44398813E-03,-0.20040369E-03,
     + -0.10658561E-02,-0.42720819E-02,-0.33492121E-03,-0.28138129E-05,
     + -0.53132302E-03,-0.73545375E-02,-0.15701946E-01, 0.33056541E-03,
     + -0.58079459E-03,-0.19298798E-02,-0.13716724E-04,-0.74686017E-03,
     + -0.22558990E-03,-0.16249549E-02,-0.39937007E-02,-0.61790115E-03,
     +  0.10705999E-04,-0.28180759E-04,-0.61445782E-03, 0.11425205E-04,
     + -0.31731302E-04,-0.60893787E-03,-0.23034422E-03,-0.13786580E-02,
     +  0.21277394E-02,-0.17203381E-02,-0.55512687E-03,-0.34243512E-03,
     +  0.23356683E-03,-0.15796241E-02,-0.16493285E-02, 0.39546159E-02,
     + -0.23654629E-02,-0.11417370E-02, 0.33296648E-03,-0.51365205E-03,
     +  0.11450103E-02, 0.26048250E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x_sr6_largex0_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31*x41    
      x_sr6_largex0_q1ex   =x_sr6_largex0_q1ex   
     9  +coeff(  9)    *x21*x31*x41    
     1  +coeff( 10)    *x21    *x43    
     2  +coeff( 11)    *x23*x32        
     3  +coeff( 12)*x11*x22    *x41    
     4  +coeff( 13)    *x23*x31*x42    
     5  +coeff( 14)    *x23    *x43    
     6  +coeff( 15)*x11    *x31        
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)*x11        *x41    
      x_sr6_largex0_q1ex   =x_sr6_largex0_q1ex   
     9  +coeff( 18)    *x21*x32*x41    
     1  +coeff( 19)    *x21*x31*x42    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x21*x32        
     5  +coeff( 23)    *x21    *x41*x51
     6  +coeff( 24)    *x21*x33        
     7  +coeff( 25)    *x21    *x42*x51
     8  +coeff( 26)*x11        *x42    
      x_sr6_largex0_q1ex   =x_sr6_largex0_q1ex   
     9  +coeff( 27)    *x21*x31*x43    
     1  +coeff( 28)    *x23    *x41*x51
     2  +coeff( 29)*x11*x22*x31*x41    
     3  +coeff( 30)*x11    *x33*x41    
     4  +coeff( 31)*x11*x22    *x42    
     5  +coeff( 32)            *x41*x51
     6  +coeff( 33)    *x23        *x51
     7  +coeff( 34)    *x21*x31*x41*x51
     8  +coeff( 35)*x11    *x32        
      x_sr6_largex0_q1ex   =x_sr6_largex0_q1ex   
     9  +coeff( 36)*x11    *x31*x41    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x21*x32*x42    
     3  +coeff( 39)    *x21*x32    *x52
     4  +coeff( 40)    *x21    *x42*x52
     5  +coeff( 41)*x11*x21    *x42    
     6  +coeff( 42)*x11    *x31*x42    
     7  +coeff( 43)*x11        *x43    
     8  +coeff( 44)    *x21*x32*x43    
      x_sr6_largex0_q1ex   =x_sr6_largex0_q1ex   
     9  +coeff( 45)    *x23    *x42*x51
     1  +coeff( 46)    *x22    *x43*x51
     2  +coeff( 47)*x11*x23    *x41    
     3  +coeff( 48)*x11*x22    *x41*x51
     4  +coeff( 49)    *x22*x31*x43*x51
     5  +coeff( 50)    *x23    *x42*x52
c
      return
      end
      function t_sr6_largex0_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.5139456E-03/
      data xmin/
     1 -0.15981E-01,-0.54082E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.54577E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.46849152E-03,-0.90568792E-02,-0.36869626E-05,-0.15646020E-06,
     + -0.65287594E-02, 0.95530973E-04, 0.29828921E-02,-0.36118752E-02,
     +  0.26805923E-03,-0.12931066E-03, 0.29221858E-03,-0.30823250E-03,
     + -0.24611130E-02,-0.47596052E-03,-0.57822581E-06,-0.13408359E-03,
     + -0.39472156E-02,-0.37820856E-02,-0.12264645E-03, 0.49706391E-03,
     +  0.19298615E-04,-0.42939070E-03,-0.51496661E-03, 0.78711477E-04,
     + -0.18262060E-02,-0.59887357E-05,-0.63513708E-03,-0.54389914E-03,
     + -0.78539386E-04,-0.91344838E-04,-0.74894386E-04,-0.16066343E-03,
     + -0.13209329E-03,-0.23539802E-03,-0.27860681E-03,-0.80786459E-03,
     + -0.49793348E-03, 0.86817157E-03, 0.10695399E-02,-0.43428282E-03,
     +  0.46454478E-03,-0.11940965E-02, 0.47942394E-03, 0.11927736E-02,
     +  0.11624560E-04,-0.33676803E-04,-0.59640086E-04,-0.23472405E-04,
     + -0.22397848E-04,-0.96120741E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sr6_largex0_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x21    *x42    
      t_sr6_largex0_q1ex   =t_sr6_largex0_q1ex   
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)*x11            *x51
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)    *x21*x31*x42    
      t_sr6_largex0_q1ex   =t_sr6_largex0_q1ex   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)*x11*x22    *x42    
     2  +coeff( 20)    *x23*x32*x41    
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)*x11    *x31*x41    
     5  +coeff( 23)*x11        *x42    
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x21*x32*x41    
     8  +coeff( 26)*x11*x22*x31        
      t_sr6_largex0_q1ex   =t_sr6_largex0_q1ex   
     9  +coeff( 27)*x11    *x31*x42    
     1  +coeff( 28)*x11        *x43    
     2  +coeff( 29)    *x21*x31    *x51
     3  +coeff( 30)*x11    *x32        
     4  +coeff( 31)*x11        *x41*x51
     5  +coeff( 32)    *x21*x33        
     6  +coeff( 33)    *x21    *x42*x51
     7  +coeff( 34)*x11    *x32*x41    
     8  +coeff( 35)    *x21*x32*x42    
      t_sr6_largex0_q1ex   =t_sr6_largex0_q1ex   
     9  +coeff( 36)    *x21*x31*x43    
     1  +coeff( 37)    *x23    *x41*x51
     2  +coeff( 38)    *x23*x31*x42    
     3  +coeff( 39)    *x21*x32*x43    
     4  +coeff( 40)    *x22    *x43*x51
     5  +coeff( 41)*x12*x21    *x43    
     6  +coeff( 42)*x12*x21    *x43*x51
     7  +coeff( 43)*x12*x21*x33*x42    
     8  +coeff( 44)*x12*x21*x32*x43*x51
      t_sr6_largex0_q1ex   =t_sr6_largex0_q1ex   
     9  +coeff( 45)                *x51
     1  +coeff( 46)*x11    *x31        
     2  +coeff( 47)*x11*x22            
     3  +coeff( 48)*x11    *x31    *x51
     4  +coeff( 49)*x11            *x52
     5  +coeff( 50)    *x22    *x42    
c
      return
      end
      function y_sr6_largex0_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1092882E-01/
      data xmin/
     1 -0.15981E-01,-0.54082E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.54577E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10404312E-01, 0.29603664E-01, 0.15751262E+00, 0.17450266E-01,
     +  0.52005990E-03, 0.77187124E-03,-0.62214425E-02,-0.19095142E-02,
     + -0.10759024E-01,-0.62344153E-02,-0.90507148E-02, 0.12279334E-02,
     + -0.10035389E-01,-0.35149094E-02,-0.13108538E-01,-0.80238190E-03,
     +  0.34523834E-02,-0.78292296E-03, 0.85866930E-04, 0.51328208E-03,
     + -0.38228079E-02,-0.27948846E-02,-0.47109458E-02,-0.20100672E-01,
     + -0.13532110E-01,-0.37827142E-01,-0.62926196E-01,-0.92164643E-01,
     +  0.12130870E-02,-0.96664904E-02, 0.10809304E+00,-0.20538296E-02,
     + -0.13342556E-02, 0.33879635E-03,-0.81010442E-03,-0.27341326E-02,
     +  0.10784095E-02,-0.32657731E-01,-0.26395700E-02,-0.63372619E-01,
     +  0.13323108E-01,-0.16285831E-01, 0.95453188E-02,-0.39994910E-01,
     + -0.28275061E-01,-0.43397844E-01,-0.43754116E-01,-0.53490434E-01,
     +  0.37467252E-01, 0.45982100E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sr6_largex0_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)*x11*x21            
      y_sr6_largex0_q1ex   =y_sr6_largex0_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)*x11*x21    *x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x22    *x41*x51
     4  +coeff( 13)*x11*x21    *x42    
     5  +coeff( 14)    *x22*x33        
     6  +coeff( 15)    *x22*x33*x41    
     7  +coeff( 16)    *x22*x32    *x52
     8  +coeff( 17)    *x22*x31*x41*x52
      y_sr6_largex0_q1ex   =y_sr6_largex0_q1ex   
     9  +coeff( 18)    *x22*x33*x41*x52
     1  +coeff( 19)        *x32        
     2  +coeff( 20)        *x33        
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)*x11*x21*x31        
     5  +coeff( 23)    *x22*x32        
     6  +coeff( 24)    *x22*x31*x41    
     7  +coeff( 25)*x11*x21*x31*x41    
     8  +coeff( 26)    *x22    *x43    
      y_sr6_largex0_q1ex   =y_sr6_largex0_q1ex   
     9  +coeff( 27)    *x22*x32*x42    
     1  +coeff( 28)    *x22    *x44    
     2  +coeff( 29)        *x34*x43    
     3  +coeff( 30)    *x22*x33*x42    
     4  +coeff( 31)    *x22*x32*x44    
     5  +coeff( 32)            *x41*x51
     6  +coeff( 33)                *x52
     7  +coeff( 34)        *x31*x42    
     8  +coeff( 35)        *x33    *x51
      y_sr6_largex0_q1ex   =y_sr6_largex0_q1ex   
     9  +coeff( 36)*x11*x21*x32        
     1  +coeff( 37)*x11*x21    *x41*x51
     2  +coeff( 38)    *x22*x31*x42    
     3  +coeff( 39)    *x22    *x42*x51
     4  +coeff( 40)    *x22*x31*x43    
     5  +coeff( 41)    *x22    *x42*x52
     6  +coeff( 42)    *x22*x34*x41    
     7  +coeff( 43)    *x23*x31*x43*x51
     8  +coeff( 44)*x11*x23    *x44    
      y_sr6_largex0_q1ex   =y_sr6_largex0_q1ex   
     9  +coeff( 45)*x11*x23    *x43*x51
     1  +coeff( 46)    *x22*x34*x44    
     2  +coeff( 47)*x12*x22    *x43*x51
     3  +coeff( 48)*x12*x22    *x42*x52
     4  +coeff( 49)*x12*x22*x32*x43*x51
     5  +coeff( 50)*x12*x22    *x42*x54
c
      return
      end
      function p_sr6_largex0_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4784687E-02/
      data xmin/
     1 -0.15981E-01,-0.54082E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.54577E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.48297592E-02,-0.12391109E-03, 0.92776474E-02, 0.68579242E-01,
     +  0.92637949E-02,-0.26809722E-02,-0.13008781E-03,-0.75419794E-03,
     + -0.34614506E-02,-0.18323755E-01,-0.24754955E-01,-0.75084256E-03,
     + -0.19607200E-03,-0.19785068E-02,-0.30449023E-02,-0.14875103E-02,
     + -0.33891920E-02,-0.74686231E-02,-0.30115476E-01,-0.21466958E-03,
     +  0.66162148E-02,-0.48930635E-03,-0.84377313E-03,-0.13637524E-02,
     + -0.65778759E-02,-0.15136166E-01,-0.69841370E-02, 0.51847182E-03,
     +  0.96279086E-03, 0.62091240E-04,-0.16840829E-02, 0.61849790E-03,
     + -0.84441938E-02,-0.50429575E-03,-0.16064500E-02,-0.38572822E-02,
     + -0.23684942E-02, 0.67617069E-02,-0.25917299E-02,-0.43298802E-02,
     +  0.27468011E-01,-0.65595419E-02,-0.13584182E-01,-0.72246878E-02,
     +  0.15737626E-02,-0.20690330E-02, 0.35756691E-02, 0.50740871E-02,
     + -0.13087937E-02,-0.16723511E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sr6_largex0_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)*x11*x21            
      p_sr6_largex0_q1ex   =p_sr6_largex0_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x22*x31*x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x22*x33        
     4  +coeff( 13)        *x32        
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)*x11*x21    *x41    
     8  +coeff( 17)    *x22*x32        
      p_sr6_largex0_q1ex   =p_sr6_largex0_q1ex   
     9  +coeff( 18)*x11*x21    *x42    
     1  +coeff( 19)    *x22    *x43    
     2  +coeff( 20)*x11*x23*x31*x41    
     3  +coeff( 21)    *x22*x33*x42    
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)                *x52
     6  +coeff( 24)*x11*x21*x31        
     7  +coeff( 25)*x11*x21*x31*x41    
     8  +coeff( 26)    *x22*x31*x42    
      p_sr6_largex0_q1ex   =p_sr6_largex0_q1ex   
     9  +coeff( 27)*x12*x22    *x42*x51
     1  +coeff( 28)            *x43    
     2  +coeff( 29)        *x32*x42    
     3  +coeff( 30)    *x21*x32    *x51
     4  +coeff( 31)*x11*x21*x32        
     5  +coeff( 32)*x11*x21    *x41*x51
     6  +coeff( 33)    *x22*x32*x41    
     7  +coeff( 34)        *x32*x43    
     8  +coeff( 35)    *x21    *x43*x51
      p_sr6_largex0_q1ex   =p_sr6_largex0_q1ex   
     9  +coeff( 36)*x11*x21    *x43    
     1  +coeff( 37)    *x22    *x43*x51
     2  +coeff( 38)    *x22    *x42*x52
     3  +coeff( 39)*x11*x23    *x42    
     4  +coeff( 40)*x11*x21    *x43*x51
     5  +coeff( 41)    *x22*x32*x43    
     6  +coeff( 42)*x11*x21*x33*x42    
     7  +coeff( 43)*x11*x23    *x43    
     8  +coeff( 44)*x11*x23    *x42*x51
      p_sr6_largex0_q1ex   =p_sr6_largex0_q1ex   
     9  +coeff( 45)*x11*x21*x32    *x53
     1  +coeff( 46)*x12*x23*x31*x41    
     2  +coeff( 47)    *x23*x31*x43*x51
     3  +coeff( 48)    *x22*x33*x41*x52
     4  +coeff( 49)*x12*x23*x33        
     5  +coeff( 50)*x12                
c
      return
      end
      function l_sr6_largex0_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.3943235E+01/
      data xmin/
     1 -0.15981E-01,-0.54082E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.54577E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.14136885E-02, 0.22541219E-02, 0.51036687E-02, 0.38845211E-02,
     +  0.58684224E-03, 0.35516478E-02, 0.61566636E-04, 0.84308366E-03,
     + -0.15540379E-03,-0.13510532E-02,-0.21510303E-02,-0.48571139E-04,
     + -0.71620423E-03,-0.21746985E-03,-0.27490156E-02,-0.49903593E-03,
     + -0.16017426E-02, 0.36434267E-04, 0.56268225E-04, 0.65703491E-04,
     + -0.31845600E-05,-0.13415824E-03,-0.74233176E-04, 0.93352974E-04,
     + -0.10004133E-03,-0.66921971E-03,-0.20528836E-02,-0.41002862E-03,
     + -0.18049017E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      l_sr6_largex0_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      l_sr6_largex0_q1ex   =l_sr6_largex0_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)*x11*x23    *x43    
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)    *x22*x31*x41    
     5  +coeff( 14)*x11*x21    *x42    
     6  +coeff( 15)    *x22    *x43    
     7  +coeff( 16)    *x22*x33*x42    
     8  +coeff( 17)*x11*x23*x33*x42    
      l_sr6_largex0_q1ex   =l_sr6_largex0_q1ex   
     9  +coeff( 18)    *x21            
     1  +coeff( 19)                *x51
     2  +coeff( 20)                *x52
     3  +coeff( 21)        *x33        
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)*x11*x21        *x51
     6  +coeff( 24)        *x32*x42    
     7  +coeff( 25)    *x22    *x41*x51
     8  +coeff( 26)    *x22*x32*x41    
      l_sr6_largex0_q1ex   =l_sr6_largex0_q1ex   
     9  +coeff( 27)    *x22*x31*x42    
     1  +coeff( 28)    *x22    *x42*x51
     2  +coeff( 29)        *x31*x41*x53
c
      return
      end
      function x_sr6_largex0_dent   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.5104202E+01/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28561212E-02,-0.12341393E+00, 0.55539236E-04,-0.70516057E-02,
     +  0.11212689E-01, 0.34373500E-02, 0.10429644E-02, 0.25268884E-01,
     + -0.23482782E-02, 0.16167004E-02, 0.18242406E-01, 0.24336595E-01,
     + -0.12842368E-02,-0.81233674E-03, 0.49659554E-02, 0.98832846E-02,
     +  0.21701463E-03,-0.84190274E-03,-0.64932735E-03, 0.80799771E-03,
     +  0.12479778E-01, 0.26771542E-01,-0.17649117E-03, 0.12528636E-02,
     +  0.39147171E-02, 0.44937228E-03, 0.48201668E-03, 0.55050524E-03,
     +  0.18015096E-02, 0.24206103E-02, 0.28863456E-02, 0.48496444E-02,
     +  0.25015036E-02,-0.13281637E-01, 0.42819274E-02,-0.12768912E-02,
     + -0.10784095E-03,-0.23877835E-03, 0.29992015E-03,-0.20447578E-02,
     + -0.77675359E-03, 0.47937073E-03,-0.18880692E-02,-0.24292960E-02,
     +  0.24231414E-03, 0.21727782E-03, 0.24718111E-02,-0.74410047E-02,
     +  0.25740720E-02, 0.24198019E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x_sr6_largex0_dent   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x42    
      x_sr6_largex0_dent   =x_sr6_largex0_dent   
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x23*x32        
     6  +coeff( 15)    *x23*x31*x42    
     7  +coeff( 16)    *x23    *x43    
     8  +coeff( 17)    *x21    *x41*x51
      x_sr6_largex0_dent   =x_sr6_largex0_dent   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x21*x32*x41    
     4  +coeff( 22)    *x21*x31*x42    
     5  +coeff( 23)*x11*x22*x31        
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)    *x21*x32        
     8  +coeff( 26)    *x22        *x51
      x_sr6_largex0_dent   =x_sr6_largex0_dent   
     9  +coeff( 27)    *x21        *x52
     1  +coeff( 28)*x11    *x31        
     2  +coeff( 29)    *x21*x33        
     3  +coeff( 30)*x11    *x31*x41    
     4  +coeff( 31)*x11        *x42    
     5  +coeff( 32)    *x21*x31*x43    
     6  +coeff( 33)    *x23    *x41*x51
     7  +coeff( 34)    *x21*x32*x43    
     8  +coeff( 35)    *x22    *x43*x51
      x_sr6_largex0_dent   =x_sr6_largex0_dent   
     9  +coeff( 36)    *x23    *x41*x52
     1  +coeff( 37)                *x51
     2  +coeff( 38)    *x22    *x41    
     3  +coeff( 39)    *x21*x31    *x51
     4  +coeff( 40)    *x23    *x41    
     5  +coeff( 41)    *x22    *x41*x51
     6  +coeff( 42)*x11    *x32        
     7  +coeff( 43)    *x21*x33*x41    
     8  +coeff( 44)    *x23    *x42    
      x_sr6_largex0_dent   =x_sr6_largex0_dent   
     9  +coeff( 45)*x11        *x41*x51
     1  +coeff( 46)*x11            *x52
     2  +coeff( 47)*x11    *x31*x42    
     3  +coeff( 48)    *x21*x33*x42    
     4  +coeff( 49)*x11        *x43    
     5  +coeff( 50)    *x23    *x42*x51
c
      return
      end
      function t_sr6_largex0_dent   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1301393E+01/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.26530172E-02, 0.55289593E-01,-0.72109336E-02, 0.37534856E-02,
     + -0.24553059E-03, 0.30108413E-02,-0.13727231E-01, 0.60797657E-03,
     +  0.67529239E-03,-0.54563629E-03,-0.89352988E-02,-0.82538446E-03,
     + -0.66932559E-03, 0.69257471E-03,-0.18371932E-02,-0.33335251E-03,
     + -0.12840098E-01,-0.11196554E-01,-0.57135080E-03, 0.14775981E-03,
     + -0.65639906E-03,-0.81616169E-03,-0.64280228E-03, 0.14222104E-02,
     + -0.55890195E-02, 0.38845508E-04,-0.98413066E-03, 0.29912342E-02,
     +  0.48030060E-05,-0.14557813E-03,-0.69924390E-04,-0.37963872E-03,
     + -0.57503948E-03, 0.32243074E-03,-0.50218677E-03,-0.82749076E-03,
     + -0.91973587E-03,-0.48443468E-03,-0.11408956E-02,-0.24222911E-02,
     + -0.10326464E-02,-0.12790263E-03,-0.39606728E-02,-0.61992917E-03,
     +  0.31176843E-02,-0.19130831E-02, 0.12601844E-03,-0.15235620E-03,
     + -0.39108694E-03,-0.16334638E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sr6_largex0_dent   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31*x41    
      t_sr6_largex0_dent   =t_sr6_largex0_dent   
     9  +coeff(  9)                *x51
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x21*x31*x42    
      t_sr6_largex0_dent   =t_sr6_largex0_dent   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x23*x32*x41    
     2  +coeff( 20)            *x41    
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)    *x23            
     5  +coeff( 23)    *x22    *x41    
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x21*x32*x41    
     8  +coeff( 26)*x11*x22*x31        
      t_sr6_largex0_dent   =t_sr6_largex0_dent   
     9  +coeff( 27)*x11*x22    *x42    
     1  +coeff( 28)    *x21*x32*x43    
     2  +coeff( 29)        *x31        
     3  +coeff( 30)        *x32        
     4  +coeff( 31)        *x31    *x51
     5  +coeff( 32)*x11    *x31        
     6  +coeff( 33)            *x43    
     7  +coeff( 34)            *x42*x51
     8  +coeff( 35)    *x21        *x52
      t_sr6_largex0_dent   =t_sr6_largex0_dent   
     9  +coeff( 36)*x11    *x31*x41    
     1  +coeff( 37)*x11        *x42    
     2  +coeff( 38)    *x21*x33        
     3  +coeff( 39)    *x22    *x42    
     4  +coeff( 40)    *x21*x31*x43    
     5  +coeff( 41)    *x23    *x41*x51
     6  +coeff( 42)        *x31*x42*x52
     7  +coeff( 43)    *x23    *x43    
     8  +coeff( 44)*x12*x21    *x41*x51
      t_sr6_largex0_dent   =t_sr6_largex0_dent   
     9  +coeff( 45)    *x23    *x42*x52
     1  +coeff( 46)    *x22*x31*x42*x52
     2  +coeff( 47)            *x41*x51
     3  +coeff( 48)    *x22*x31        
     4  +coeff( 49)        *x31*x42    
     5  +coeff( 50)    *x21*x31    *x51
c
      return
      end
      function y_sr6_largex0_dent   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.7575093E-02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.77847666E-02, 0.60088001E-02, 0.10994225E+00, 0.18325070E-01,
     + -0.97202610E-04, 0.36697980E-03, 0.26029334E-02, 0.17417682E-01,
     +  0.50875526E-02,-0.56712185E-02,-0.15827282E-02,-0.11066395E-01,
     + -0.98223693E-03, 0.14256722E-03,-0.28721677E-03,-0.66757011E-02,
     + -0.32844141E-01,-0.29142179E-01,-0.20726905E-02,-0.11899984E-01,
     + -0.24008586E-02,-0.54859560E-01, 0.26178476E-03,-0.50056824E-02,
     +  0.25672574E-01,-0.57153875E-03, 0.13986557E-02, 0.36560951E-03,
     + -0.57779965E-02,-0.71873078E-02,-0.16978725E-02,-0.34719252E-02,
     + -0.10221455E-01,-0.17566493E-01,-0.38093720E-01,-0.10571231E-01,
     + -0.37986025E-01,-0.10730065E-01,-0.92522260E-02, 0.56689005E-01,
     + -0.93369856E-02,-0.23328403E-01,-0.26871009E-01, 0.12191226E-02,
     + -0.85300731E-03,-0.92128367E-03,-0.48509709E-03, 0.27966381E-02,
     + -0.34006936E-02,-0.30468125E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sr6_largex0_dent   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      y_sr6_largex0_dent   =y_sr6_largex0_dent   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x22            
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)        *x34        
     6  +coeff( 15)    *x21*x33        
     7  +coeff( 16)*x11*x21    *x41    
     8  +coeff( 17)    *x22*x31*x41    
      y_sr6_largex0_dent   =y_sr6_largex0_dent   
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x22    *x41*x51
     2  +coeff( 20)*x11*x21    *x42    
     3  +coeff( 21)    *x22*x33        
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)*x11*x23*x31        
     6  +coeff( 24)*x11*x21*x33*x41    
     7  +coeff( 25)    *x22*x33*x42    
     8  +coeff( 26)    *x21            
      y_sr6_largex0_dent   =y_sr6_largex0_dent   
     9  +coeff( 27)                *x52
     1  +coeff( 28)        *x33        
     2  +coeff( 29)    *x22*x31        
     3  +coeff( 30)    *x22*x32        
     4  +coeff( 31)    *x22*x31    *x51
     5  +coeff( 32)*x11*x21*x32        
     6  +coeff( 33)*x11*x21*x31*x41    
     7  +coeff( 34)    *x22*x32*x41    
     8  +coeff( 35)    *x22*x31*x42    
      y_sr6_largex0_dent   =y_sr6_largex0_dent   
     9  +coeff( 36)    *x22    *x42*x51
     1  +coeff( 37)    *x22    *x44    
     2  +coeff( 38)    *x22    *x43*x51
     3  +coeff( 39)*x11*x21    *x43*x51
     4  +coeff( 40)    *x22*x32*x43    
     5  +coeff( 41)    *x22*x33*x41*x51
     6  +coeff( 42)    *x22*x33*x43    
     7  +coeff( 43)*x11*x23    *x44    
     8  +coeff( 44)    *x21*x31        
      y_sr6_largex0_dent   =y_sr6_largex0_dent   
     9  +coeff( 45)            *x41*x52
     1  +coeff( 46)    *x21    *x41*x51
     2  +coeff( 47)*x12                
     3  +coeff( 48)        *x33*x41    
     4  +coeff( 49)        *x31*x43    
     5  +coeff( 50)*x11*x21*x31        
c
      return
      end
      function p_sr6_largex0_dent   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1770968E-02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16170467E-02,-0.45161862E-04,-0.23653872E-01,-0.17326403E-02,
     +  0.18800652E-04, 0.32794684E-04,-0.35949789E-02, 0.36334965E-04,
     + -0.80100239E-04,-0.43939779E-03, 0.45521343E-02, 0.59795263E-03,
     +  0.33641863E-03, 0.14062065E-04,-0.13937195E-02,-0.25727708E-03,
     + -0.90434813E-04,-0.40754289E-03, 0.36867084E-04, 0.76237228E-03,
     +  0.50168494E-02, 0.34221255E-04,-0.90648056E-04, 0.36402384E-02,
     +  0.54441375E-03, 0.32296961E-04,-0.33807250E-04,-0.10521655E-03,
     +  0.17377299E-02, 0.16143911E-03,-0.66928039E-02, 0.44019122E-03,
     +  0.70134073E-03, 0.30065794E-03, 0.15506448E-02, 0.36903084E-02,
     + -0.56678947E-03, 0.84460847E-03, 0.66103181E-03,-0.16017152E-03,
     +  0.49781157E-02, 0.37410526E-03, 0.28202008E-03, 0.58530003E-03,
     +  0.28214819E-03, 0.74702676E-03,-0.48631930E-03,-0.19890038E-03,
     +  0.25699148E-03,-0.55348920E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sr6_largex0_dent   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)        *x32        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sr6_largex0_dent   =p_sr6_largex0_dent   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)        *x33        
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)        *x31*x42    
     8  +coeff( 17)        *x31    *x52
      p_sr6_largex0_dent   =p_sr6_largex0_dent   
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)*x12    *x31        
     2  +coeff( 20)    *x22*x32        
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)        *x33    *x51
     5  +coeff( 23)    *x22*x33        
     6  +coeff( 24)    *x22*x31*x42    
     7  +coeff( 25)        *x33*x42    
     8  +coeff( 26)    *x22*x31    *x52
      p_sr6_largex0_dent   =p_sr6_largex0_dent   
     9  +coeff( 27)        *x33    *x52
     1  +coeff( 28)*x12    *x33        
     2  +coeff( 29)*x11*x23    *x42    
     3  +coeff( 30)    *x22*x33    *x52
     4  +coeff( 31)        *x31        
     5  +coeff( 32)    *x22            
     6  +coeff( 33)        *x31    *x51
     7  +coeff( 34)*x11        *x41    
     8  +coeff( 35)    *x23    *x41    
      p_sr6_largex0_dent   =p_sr6_largex0_dent   
     9  +coeff( 36)    *x22*x31*x41    
     1  +coeff( 37)    *x22        *x51
     2  +coeff( 38)*x11*x21    *x41    
     3  +coeff( 39)    *x21    *x43    
     4  +coeff( 40)*x11*x23            
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)*x11*x21            
     7  +coeff( 43)    *x23            
     8  +coeff( 44)    *x21*x31*x41    
      p_sr6_largex0_dent   =p_sr6_largex0_dent   
     9  +coeff( 45)        *x32*x41    
     1  +coeff( 46)    *x21    *x42    
     2  +coeff( 47)            *x43    
     3  +coeff( 48)    *x21    *x41*x51
     4  +coeff( 49)*x11*x21*x31        
     5  +coeff( 50)    *x22    *x41*x51
c
      return
      end
      function l_sr6_largex0_dent   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.1074794E+02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21030321E-02,-0.23504019E+00, 0.22353357E-01, 0.16268065E-01,
     +  0.72698150E-03,-0.13091178E-01, 0.40749079E-02, 0.27042453E-02,
     +  0.51274365E-02, 0.44445191E-01,-0.74044988E-03, 0.18893302E-02,
     + -0.37158146E-02,-0.25899180E-02,-0.17621707E-02, 0.33163026E-01,
     +  0.41634619E-01, 0.10562832E-01, 0.68504023E-02, 0.10006211E-01,
     +  0.17712738E-02, 0.20179527E-02, 0.12452649E-02, 0.30769436E-02,
     +  0.21298698E-02, 0.39546858E-03, 0.15940947E-01, 0.48614029E-01,
     +  0.24747355E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      l_sr6_largex0_dent   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)            *x41    
     8  +coeff(  8)    *x21*x31        
      l_sr6_largex0_dent   =l_sr6_largex0_dent   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x23*x31*x41    
     3  +coeff( 12)        *x31        
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x21    *x43    
      l_sr6_largex0_dent   =l_sr6_largex0_dent   
     9  +coeff( 18)    *x23*x32        
     1  +coeff( 19)    *x23*x31*x42    
     2  +coeff( 20)    *x23    *x43    
     3  +coeff( 21)            *x41*x51
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)*x11        *x42    
     8  +coeff( 26)    *x23*x31        
      l_sr6_largex0_dent   =l_sr6_largex0_dent   
     9  +coeff( 27)    *x21*x32*x41    
     1  +coeff( 28)    *x21*x31*x42    
     2  +coeff( 29)*x11*x22*x31        
c
      return
      end
      function x_sr6_largex0_dext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1539440E-01/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.11278605E-01, 0.33763430E+00, 0.10485929E-03, 0.13341391E+00,
     + -0.43062653E-01,-0.89220371E-03, 0.44704121E-01,-0.79002641E-01,
     +  0.12507457E-01,-0.37748928E-02,-0.53626645E-01,-0.33874370E-02,
     + -0.97651398E-02,-0.58529396E-02,-0.68941981E-01,-0.24443969E-01,
     + -0.38941663E-01,-0.46742051E-02,-0.23832233E-02, 0.15863281E-02,
     + -0.14394124E-02, 0.24494226E-02,-0.10050169E-01,-0.38595866E-01,
     + -0.75514607E-01, 0.10573974E-03, 0.14844115E-02,-0.24871359E-03,
     + -0.15895215E-02,-0.15408089E-02,-0.13088167E-02,-0.92099253E-02,
     + -0.50090812E-02, 0.50776191E-02, 0.61394875E-02,-0.93489895E-02,
     + -0.97725829E-02,-0.20627582E-01, 0.40322464E-01,-0.19398848E-01,
     + -0.17999855E-03,-0.27894820E-02,-0.31524880E-02,-0.49002911E-03,
     +  0.81228521E-02,-0.13321325E-02,-0.92779408E-03,-0.20838929E-02,
     +  0.11358637E-02,-0.36167025E-02,-0.12116619E-01,-0.94675934E-02,
     + -0.86964938E-04,-0.59596426E-03,-0.31096938E-02, 0.44424352E-03,
     + -0.56790388E-02, 0.17811893E-01, 0.27180070E-02,-0.23384231E-04,
     +  0.13699408E-01, 0.57962728E-02, 0.18287802E-01,-0.12272787E-02,
     + -0.45086746E-02,-0.26453456E-01, 0.22305626E-01, 0.64225378E-02,
     +  0.12240128E-01,-0.28704936E-01, 0.14965049E-01,-0.11689475E-01,
     +  0.22379730E-01,-0.45566775E-01,-0.50804850E-01, 0.22138219E-01,
     + -0.19728901E-01,-0.67051046E-03,-0.14293352E-02,-0.64817717E-03,
     + -0.16101240E-02,-0.23678467E-02,-0.13634162E-02, 0.17138920E-03,
     + -0.17482173E-02, 0.37970633E-03, 0.57899149E-03, 0.67986234E-03,
     + -0.81743562E-03,-0.23440872E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x_sr6_largex0_dext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x21    *x42    
      x_sr6_largex0_dext   =x_sr6_largex0_dext   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21    *x43    
     7  +coeff( 16)    *x23*x31*x42    
     8  +coeff( 17)    *x23    *x43    
      x_sr6_largex0_dext   =x_sr6_largex0_dext   
     9  +coeff( 18)            *x42    
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)    *x21*x31    *x51
     4  +coeff( 22)    *x23*x31        
     5  +coeff( 23)*x11        *x42    
     6  +coeff( 24)    *x21*x32*x41    
     7  +coeff( 25)    *x21*x31*x42    
     8  +coeff( 26)*x11*x22*x31        
      x_sr6_largex0_dext   =x_sr6_largex0_dext   
     9  +coeff( 27)*x11    *x33*x41    
     1  +coeff( 28)            *x41*x51
     2  +coeff( 29)*x11*x21            
     3  +coeff( 30)*x11    *x31        
     4  +coeff( 31)*x11*x22            
     5  +coeff( 32)*x11    *x31*x41    
     6  +coeff( 33)    *x21*x33        
     7  +coeff( 34)    *x21    *x42*x51
     8  +coeff( 35)    *x23*x31*x41    
      x_sr6_largex0_dext   =x_sr6_largex0_dext   
     9  +coeff( 36)    *x21*x31*x43    
     1  +coeff( 37)    *x23*x31*x41*x51
     2  +coeff( 38)    *x23    *x42*x51
     3  +coeff( 39)    *x21*x32*x43    
     4  +coeff( 40)*x12*x21    *x42*x51
     5  +coeff( 41)        *x31        
     6  +coeff( 42)        *x31*x41    
     7  +coeff( 43)    *x24            
     8  +coeff( 44)        *x31*x41*x51
      x_sr6_largex0_dext   =x_sr6_largex0_dext   
     9  +coeff( 45)    *x23    *x41    
     1  +coeff( 46)*x11    *x32        
     2  +coeff( 47)*x11        *x41*x51
     3  +coeff( 48)    *x21*x32    *x51
     4  +coeff( 49)        *x31*x43    
     5  +coeff( 50)            *x43*x51
     6  +coeff( 51)*x11    *x31*x42    
     7  +coeff( 52)*x11        *x43    
     8  +coeff( 53)    *x22    *x42*x51
      x_sr6_largex0_dext   =x_sr6_largex0_dext   
     9  +coeff( 54)    *x21*x32*x42    
     1  +coeff( 55)    *x21*x32    *x52
     2  +coeff( 56)    *x21    *x42*x52
     3  +coeff( 57)    *x22    *x43*x51
     4  +coeff( 58)    *x21*x33*x42    
     5  +coeff( 59)*x12*x21    *x41*x51
     6  +coeff( 60)*x12*x23        *x51
     7  +coeff( 61)    *x23    *x42*x52
     8  +coeff( 62)*x11    *x33*x42    
      x_sr6_largex0_dext   =x_sr6_largex0_dext   
     9  +coeff( 63)    *x21    *x43*x53
     1  +coeff( 64)*x11*x22*x31    *x53
     2  +coeff( 65)*x12*x21*x31*x43    
     3  +coeff( 66)*x12*x21    *x43*x51
     4  +coeff( 67)*x11*x23    *x43*x51
     5  +coeff( 68)*x12*x23        *x53
     6  +coeff( 69)*x11*x21*x31*x42*x53
     7  +coeff( 70)    *x23    *x43*x53
     8  +coeff( 71)*x11*x24*x31*x43    
      x_sr6_largex0_dext   =x_sr6_largex0_dext   
     9  +coeff( 72)*x12*x22    *x42*x52
     1  +coeff( 73)*x12*x21*x32*x42*x51
     2  +coeff( 74)*x12*x23*x32*x42    
     3  +coeff( 75)*x12*x23*x31*x43    
     4  +coeff( 76)*x11*x24    *x43*x52
     5  +coeff( 77)*x12*x21*x33*x42*x51
     6  +coeff( 78)        *x32        
     7  +coeff( 79)    *x23            
     8  +coeff( 80)    *x21        *x52
      x_sr6_largex0_dext   =x_sr6_largex0_dext   
     9  +coeff( 81)        *x32*x41    
     1  +coeff( 82)        *x31*x42    
     2  +coeff( 83)            *x43    
     3  +coeff( 84)        *x32    *x51
     4  +coeff( 85)            *x42*x51
     5  +coeff( 86)                *x53
     6  +coeff( 87)*x11*x21        *x51
     7  +coeff( 88)    *x22*x32        
     8  +coeff( 89)*x11*x22        *x51
      x_sr6_largex0_dext   =x_sr6_largex0_dext   
     9  +coeff( 90)        *x31*x42*x51
c
      return
      end
      function t_sr6_largex0_dext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.5302829E+00/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.34623628E-02,-0.77066280E-01, 0.15264535E-04, 0.22816109E-01,
     +  0.68826457E-02,-0.32565545E-02, 0.20718708E-03,-0.39863400E-02,
     +  0.83329272E-03,-0.12763551E-02, 0.15358668E-01,-0.61361108E-03,
     +  0.93690958E-03,-0.83010638E-03, 0.10173690E-01, 0.14062334E-01,
     + -0.26379398E-03, 0.14164359E-02, 0.15839681E-02, 0.56242724E-02,
     +  0.86431316E-03, 0.15145477E-02,-0.18795539E-03, 0.76985830E-02,
     +  0.15659805E-01, 0.67653687E-05, 0.88050961E-03,-0.17845219E-03,
     +  0.10864295E-03,-0.21281824E-03, 0.46454495E-03,-0.17681439E-04,
     +  0.21158089E-02, 0.85447729E-03, 0.53785474E-03, 0.13456357E-02,
     +  0.65759744E-03, 0.95023256E-03, 0.11938622E-02, 0.29081574E-02,
     + -0.42914087E-02,-0.31948560E-02, 0.29240882E-02, 0.16904715E-03,
     + -0.24450914E-03, 0.24152825E-03, 0.22961527E-03, 0.14644224E-03,
     + -0.14793624E-02, 0.44294610E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sr6_largex0_dext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_sr6_largex0_dext   =t_sr6_largex0_dext   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23*x31*x41    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x23*x32        
      t_sr6_largex0_dext   =t_sr6_largex0_dext   
     9  +coeff( 18)    *x23    *x41*x51
     1  +coeff( 19)    *x23*x31*x42    
     2  +coeff( 20)    *x23    *x43    
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)*x11        *x42    
     5  +coeff( 23)    *x23*x31        
     6  +coeff( 24)    *x21*x32*x41    
     7  +coeff( 25)    *x21*x31*x42    
     8  +coeff( 26)*x11    *x33        
      t_sr6_largex0_dext   =t_sr6_largex0_dext   
     9  +coeff( 27)    *x21*x32    *x52
     1  +coeff( 28)*x11    *x33*x41    
     2  +coeff( 29)        *x31        
     3  +coeff( 30)        *x32        
     4  +coeff( 31)*x11    *x31        
     5  +coeff( 32)    *x23            
     6  +coeff( 33)    *x21*x32        
     7  +coeff( 34)        *x31*x41*x51
     8  +coeff( 35)*x11*x22            
      t_sr6_largex0_dext   =t_sr6_largex0_dext   
     9  +coeff( 36)*x11    *x31*x41    
     1  +coeff( 37)    *x21*x33        
     2  +coeff( 38)    *x21    *x42*x51
     3  +coeff( 39)    *x22    *x43    
     4  +coeff( 40)    *x21*x31*x43    
     5  +coeff( 41)    *x21*x32*x43    
     6  +coeff( 42)    *x22*x31*x43*x51
     7  +coeff( 43)*x12*x21    *x43*x51
     8  +coeff( 44)        *x31    *x51
      t_sr6_largex0_dext   =t_sr6_largex0_dext   
     9  +coeff( 45)    *x22        *x51
     1  +coeff( 46)            *x42*x51
     2  +coeff( 47)*x11    *x32        
     3  +coeff( 48)*x11            *x52
     4  +coeff( 49)    *x23    *x41    
     5  +coeff( 50)    *x22    *x41*x51
c
      return
      end
      function y_sr6_largex0_dext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3901048E-03/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.87675918E-03,-0.37318029E-01, 0.11918847E-01, 0.18258056E-01,
     + -0.73099777E-03, 0.24226453E-03, 0.10361054E-02, 0.69704669E-03,
     +  0.86100399E-02, 0.55772576E-01, 0.26355029E-03, 0.60180812E-02,
     + -0.35630218E-02,-0.47004569E-01,-0.71629691E-02,-0.67112236E-02,
     + -0.85134367E-02,-0.17286576E-02,-0.26625998E-02,-0.55930922E-02,
     +  0.84134890E-03, 0.52731214E-02,-0.39999839E-02, 0.71618212E-02,
     +  0.65227658E-02, 0.89134456E-03,-0.48987399E-03,-0.50335797E-03,
     + -0.56206207E-02,-0.61161956E-02,-0.20533723E-02,-0.17374177E-02,
     + -0.23142113E-02,-0.39956582E-03,-0.42557459E-01,-0.58348747E-02,
     +  0.51952165E-03, 0.32436217E-02,-0.22090815E-02,-0.28985261E-02,
     + -0.36958058E-02,-0.55284025E-02,-0.24061987E-01,-0.11186434E-01,
     + -0.15460384E-02,-0.76725716E-02, 0.49594097E-03, 0.25044580E-02,
     +  0.47939625E-02,-0.14972850E-01,-0.56806016E-02, 0.12649963E-01,
     +  0.11593542E-01,-0.48445850E-02,-0.17494192E-01,-0.12676207E-01,
     +  0.37435394E-01,-0.14290005E-02, 0.44198185E-02, 0.73128757E-02,
     +  0.39029267E-01, 0.96082382E-01,-0.46204060E-01,-0.99499838E-03,
     +  0.12181818E-02,-0.10751558E-01,-0.34035284E-02,-0.86315079E-02,
     + -0.56688031E-02,-0.90416614E-03,-0.50489337E-03,-0.28318281E-02,
     +  0.19712641E-02,-0.24536401E-02,-0.21394088E-02,-0.10346181E-01,
     + -0.89478698E-02, 0.17731424E-02, 0.13387456E-01,-0.60160834E-01,
     +  0.37293362E-02,-0.12570473E-02, 0.10221185E-01,-0.76343492E-01,
     + -0.75722545E-01, 0.13693120E-01, 0.53008962E-02,-0.34999605E-02,
     +  0.89178160E-02,-0.61997954E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x35 = x34*x3
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
      y_sr6_largex0_dext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sr6_largex0_dext   =y_sr6_largex0_dext   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)        *x31*x42    
     8  +coeff( 17)            *x43    
      y_sr6_largex0_dext   =y_sr6_largex0_dext   
     9  +coeff( 18)        *x31*x41*x51
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)    *x22            
     3  +coeff( 21)*x11    *x31        
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)            *x41*x52
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x21    *x42    
     8  +coeff( 26)*x11            *x51
      y_sr6_largex0_dext   =y_sr6_largex0_dext   
     9  +coeff( 27)                *x53
     1  +coeff( 28)    *x21*x31    *x51
     2  +coeff( 29)    *x21    *x41*x51
     3  +coeff( 30)        *x32*x42    
     4  +coeff( 31)            *x44    
     5  +coeff( 32)*x11*x21            
     6  +coeff( 33)        *x31*x42*x51
     7  +coeff( 34)            *x43*x51
     8  +coeff( 35)    *x22    *x41    
      y_sr6_largex0_dext   =y_sr6_largex0_dext   
     9  +coeff( 36)    *x22        *x51
     1  +coeff( 37)    *x21*x32*x41    
     2  +coeff( 38)    *x23            
     3  +coeff( 39)        *x34*x41    
     4  +coeff( 40)*x11*x21*x31        
     5  +coeff( 41)*x11*x21    *x41    
     6  +coeff( 42)    *x22*x32        
     7  +coeff( 43)    *x22*x31*x41    
     8  +coeff( 44)    *x22    *x42    
      y_sr6_largex0_dext   =y_sr6_largex0_dext   
     9  +coeff( 45)*x11*x21        *x51
     1  +coeff( 46)    *x22    *x41*x51
     2  +coeff( 47)    *x21*x33    *x51
     3  +coeff( 48)        *x32*x44    
     4  +coeff( 49)        *x31*x45    
     5  +coeff( 50)*x11*x21    *x42    
     6  +coeff( 51)    *x22*x33        
     7  +coeff( 52)    *x22*x31*x42    
     8  +coeff( 53)    *x22    *x43    
      y_sr6_largex0_dext   =y_sr6_largex0_dext   
     9  +coeff( 54)*x11*x21    *x41*x51
     1  +coeff( 55)    *x22    *x42*x51
     2  +coeff( 56)    *x22*x31*x42*x51
     3  +coeff( 57)    *x23    *x43    
     4  +coeff( 58)*x11*x21*x33*x41    
     5  +coeff( 59)    *x22*x35        
     6  +coeff( 60)*x11*x22    *x41*x52
     7  +coeff( 61)    *x23*x32*x43    
     8  +coeff( 62)    *x23*x31*x44    
      y_sr6_largex0_dext   =y_sr6_largex0_dext   
     9  +coeff( 63)*x12*x22    *x43*x51
     1  +coeff( 64)        *x32    *x51
     2  +coeff( 65)    *x21*x32        
     3  +coeff( 66)        *x31*x43    
     4  +coeff( 67)        *x32*x41*x51
     5  +coeff( 68)    *x22*x31        
     6  +coeff( 69)    *x21    *x43    
     7  +coeff( 70)    *x21    *x42*x51
     8  +coeff( 71)*x12                
      y_sr6_largex0_dext   =y_sr6_largex0_dext   
     9  +coeff( 72)        *x32*x43    
     1  +coeff( 73)            *x45    
     2  +coeff( 74)    *x22*x31    *x51
     3  +coeff( 75)*x11*x21*x32        
     4  +coeff( 76)*x11*x21*x31*x41    
     5  +coeff( 77)    *x22*x31*x41*x51
     6  +coeff( 78)    *x21*x33*x42    
     7  +coeff( 79)    *x21    *x45    
     8  +coeff( 80)    *x22    *x44    
      y_sr6_largex0_dext   =y_sr6_largex0_dext   
     9  +coeff( 81)    *x21*x31*x41*x53
     1  +coeff( 82)        *x34*x44    
     2  +coeff( 83)    *x23    *x42*x51
     3  +coeff( 84)    *x22*x31*x44    
     4  +coeff( 85)    *x22    *x45    
     5  +coeff( 86)    *x23    *x44    
     6  +coeff( 87)*x12    *x32*x41*x51
     7  +coeff( 88)*x12*x22    *x41    
     8  +coeff( 89)    *x22*x33*x43    
      y_sr6_largex0_dext   =y_sr6_largex0_dext   
     9  +coeff( 90)    *x22*x31*x45    
c
      return
      end
      function p_sr6_largex0_dext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1183870E-02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.88239467E-03,-0.88209235E-04,-0.91338726E-02,-0.13035971E-01,
     +  0.13565923E-02,-0.10202427E-02,-0.69971116E-04,-0.99147316E-02,
     + -0.40757706E-03,-0.13771558E-02, 0.13260097E-02, 0.10056433E-01,
     +  0.12394397E-02, 0.11001937E-02,-0.65729078E-02,-0.14504023E-02,
     +  0.13423369E-02,-0.52428717E-03,-0.33625055E-03,-0.24881616E-03,
     +  0.42401534E-03,-0.60061540E-03, 0.48741110E-03,-0.15470891E-02,
     + -0.13702147E-02,-0.72551961E-03, 0.34592484E-03,-0.16163583E-02,
     +  0.27431638E-02, 0.27139743E-02,-0.21314602E-02, 0.20578101E-02,
     +  0.19046958E-03, 0.16277419E-03,-0.80586068E-03, 0.20079191E-03,
     + -0.58821111E-03,-0.42026633E-03,-0.25757562E-03, 0.48448140E-03,
     + -0.65423339E-03, 0.53314312E-03, 0.84785226E-03,-0.16228579E-02,
     + -0.13351366E-02, 0.33367733E-02, 0.45303893E-02, 0.14706706E-02,
     + -0.32917082E-02, 0.37100715E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sr6_largex0_dext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      p_sr6_largex0_dext   =p_sr6_largex0_dext   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)    *x23    *x41    
      p_sr6_largex0_dext   =p_sr6_largex0_dext   
     9  +coeff( 18)    *x22            
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x22*x31        
     5  +coeff( 23)    *x21    *x42    
     6  +coeff( 24)            *x43    
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)            *x42*x51
      p_sr6_largex0_dext   =p_sr6_largex0_dext   
     9  +coeff( 27)    *x23*x31        
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)    *x23*x31*x41    
     3  +coeff( 30)    *x23    *x42    
     4  +coeff( 31)*x11*x21    *x42*x51
     5  +coeff( 32)*x11*x22    *x42*x52
     6  +coeff( 33)*x11    *x31        
     7  +coeff( 34)*x11            *x51
     8  +coeff( 35)        *x31*x42    
      p_sr6_largex0_dext   =p_sr6_largex0_dext   
     9  +coeff( 36)    *x21*x31    *x51
     1  +coeff( 37)        *x31*x41*x51
     2  +coeff( 38)            *x41*x52
     3  +coeff( 39)*x11*x21*x31        
     4  +coeff( 40)    *x21    *x43    
     5  +coeff( 41)    *x22*x31    *x51
     6  +coeff( 42)    *x23*x32        
     7  +coeff( 43)    *x23    *x41*x51
     8  +coeff( 44)    *x22    *x42*x51
      p_sr6_largex0_dext   =p_sr6_largex0_dext   
     9  +coeff( 45)*x11*x21*x31*x41*x51
     1  +coeff( 46)    *x23*x31*x42    
     2  +coeff( 47)    *x23    *x43    
     3  +coeff( 48)    *x23    *x42*x51
     4  +coeff( 49)*x12*x22    *x43*x51
     5  +coeff( 50)*x11                
c
      return
      end
      function l_sr6_largex0_dext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1737264E+02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.22410078E-01, 0.48208794E+00, 0.99035643E-01,-0.54861072E-01,
     +  0.31376056E-01, 0.13955897E-02, 0.43195229E-01,-0.10696811E+00,
     +  0.10592190E-01,-0.41373502E-02,-0.74228756E-01, 0.48506171E-02,
     + -0.56583621E-02, 0.37373276E-02,-0.14080929E-01, 0.26084455E-02,
     + -0.10785713E+00,-0.10071434E+00, 0.12758988E-01, 0.20279195E-02,
     + -0.74644515E-03, 0.70345500E-02,-0.59539907E-01, 0.16465285E-05,
     + -0.24585698E-02,-0.16164441E-02, 0.10615343E-02, 0.61404482E-02,
     + -0.18591591E-02,-0.17560989E-02,-0.48758774E-02,-0.12186220E-01,
     + -0.14740883E-01,-0.67108078E-02,-0.14531055E-01, 0.77140080E-02,
     + -0.20158757E-01,-0.51029772E-02,-0.49803745E-01,-0.22991451E-02,
     +  0.42494781E-01,-0.16649161E-02, 0.16047517E-02,-0.31414088E-02,
     + -0.18837932E-02, 0.75859844E-03,-0.10611714E-02,-0.42795169E-03,
     + -0.86867111E-02, 0.21789980E-02,-0.19225092E-02,-0.18442465E-01,
     + -0.16211584E-01,-0.11353324E-01,-0.20156484E-01,-0.37774447E-03,
     +  0.83972085E-02, 0.47307666E-01,-0.16264049E-01,-0.18024370E-01,
     + -0.13062424E-01, 0.20751392E-02,-0.69664686E-03,-0.29826360E-01,
     + -0.21844676E-01,-0.21597365E-01,-0.15904643E-01,-0.16979113E-01,
     + -0.42648078E-02, 0.33161761E-02, 0.64277854E-02,-0.77992245E-02,
     +  0.51933054E-01,-0.24963228E-01, 0.23612117E-01, 0.36050279E-01,
     + -0.18516714E-01,-0.48112452E-01,-0.33169528E-03, 0.16269932E-02,
     + -0.90821221E-03,-0.15325675E-02,-0.16729945E-02, 0.13295324E-02,
     +  0.18969174E-02,-0.76538877E-03,-0.22672785E-02,-0.30382208E-02,
     + -0.76334164E-02, 0.10039265E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sr6_largex0_dext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x21    *x42    
      l_sr6_largex0_dext   =l_sr6_largex0_dext   
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)            *x41    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x21*x31*x42    
      l_sr6_largex0_dext   =l_sr6_largex0_dext   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x23*x32*x41    
     2  +coeff( 20)        *x31        
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x32*x41    
     6  +coeff( 24)*x11*x22*x31        
     7  +coeff( 25)*x11*x22    *x42    
     8  +coeff( 26)*x11*x22*x33*x41    
      l_sr6_largex0_dext   =l_sr6_largex0_dext   
     9  +coeff( 27)        *x31*x41    
     1  +coeff( 28)            *x42    
     2  +coeff( 29)                *x52
     3  +coeff( 30)*x11    *x31        
     4  +coeff( 31)    *x22    *x41    
     5  +coeff( 32)*x11    *x31*x41    
     6  +coeff( 33)*x11        *x42    
     7  +coeff( 34)    *x21*x33        
     8  +coeff( 35)    *x22    *x42    
      l_sr6_largex0_dext   =l_sr6_largex0_dext   
     9  +coeff( 36)    *x21    *x42*x51
     1  +coeff( 37)    *x21*x31*x43    
     2  +coeff( 38)*x11*x22*x32        
     3  +coeff( 39)    *x23    *x43    
     4  +coeff( 40)    *x21*x33    *x53
     5  +coeff( 41)*x12*x23    *x41*x51
     6  +coeff( 42)    *x22*x31        
     7  +coeff( 43)    *x22        *x51
     8  +coeff( 44)            *x42*x51
      l_sr6_largex0_dext   =l_sr6_largex0_dext   
     9  +coeff( 45)    *x21        *x52
     1  +coeff( 46)*x11*x21        *x51
     2  +coeff( 47)*x11        *x41*x51
     3  +coeff( 48)*x11            *x52
     4  +coeff( 49)    *x22*x31*x41    
     5  +coeff( 50)        *x31*x43    
     6  +coeff( 51)    *x21*x31*x41*x51
     7  +coeff( 52)*x11    *x31*x42    
     8  +coeff( 53)*x11        *x43    
      l_sr6_largex0_dext   =l_sr6_largex0_dext   
     9  +coeff( 54)    *x21*x32*x42    
     1  +coeff( 55)    *x23    *x41*x51
     2  +coeff( 56)    *x21*x32    *x52
     3  +coeff( 57)    *x21    *x42*x52
     4  +coeff( 58)    *x21*x32*x43    
     5  +coeff( 59)*x12*x21    *x41*x51
     6  +coeff( 60)    *x23    *x42*x51
     7  +coeff( 61)    *x22    *x43*x51
     8  +coeff( 62)        *x32*x41*x53
      l_sr6_largex0_dext   =l_sr6_largex0_dext   
     9  +coeff( 63)*x11    *x33*x42    
     1  +coeff( 64)*x12*x21*x31*x42    
     2  +coeff( 65)    *x23*x31*x43    
     3  +coeff( 66)*x12*x21    *x42*x51
     4  +coeff( 67)    *x22*x31*x42*x52
     5  +coeff( 68)    *x21*x32*x42*x52
     6  +coeff( 69)        *x32*x43*x52
     7  +coeff( 70)*x12*x21        *x53
     8  +coeff( 71)        *x32*x42*x53
      l_sr6_largex0_dext   =l_sr6_largex0_dext   
     9  +coeff( 72)*x12*x21*x31*x43    
     1  +coeff( 73)    *x23    *x43*x52
     2  +coeff( 74)    *x21*x32*x43*x52
     3  +coeff( 75)*x11*x23    *x43*x51
     4  +coeff( 76)*x12*x21*x33*x42    
     5  +coeff( 77)*x12*x23*x31*x41*x51
     6  +coeff( 78)*x12*x23    *x43*x51
     7  +coeff( 79)        *x32        
     8  +coeff( 80)            *x41*x51
      l_sr6_largex0_dext   =l_sr6_largex0_dext   
     9  +coeff( 81)    *x23            
     1  +coeff( 82)    *x21*x31    *x51
     2  +coeff( 83)*x11    *x32        
     3  +coeff( 84)    *x23*x31        
     4  +coeff( 85)    *x23        *x51
     5  +coeff( 86)    *x21*x32    *x51
     6  +coeff( 87)            *x43*x51
     7  +coeff( 88)    *x21    *x41*x52
     8  +coeff( 89)*x11    *x32*x41    
      l_sr6_largex0_dext   =l_sr6_largex0_dext   
     9  +coeff( 90)*x11*x21*x31    *x51
c
      return
      end
      function x_sr6_largex0_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.7152329E-02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.32037555E-02, 0.21399041E+00,-0.22212848E-03, 0.13911122E+00,
     + -0.30067669E-01, 0.18834535E-01,-0.77014638E-03, 0.34679890E-01,
     + -0.53513803E-01, 0.87997451E-03,-0.25298316E-02,-0.59280633E-02,
     + -0.35193540E-01,-0.32967334E-02,-0.30404485E-02,-0.65180990E-02,
     + -0.20268965E-02,-0.44931952E-01,-0.91526061E-02,-0.25754236E-01,
     + -0.12878266E-02,-0.62220730E-03,-0.29555525E-03, 0.11109583E-02,
     + -0.12449010E-02, 0.87136560E-03, 0.48919101E-02,-0.23641706E-01,
     + -0.50605826E-01,-0.25225250E-03, 0.14642256E-02,-0.62666899E-04,
     +  0.65317383E-03,-0.66571170E-02,-0.73833750E-02,-0.23837145E-02,
     +  0.28870425E-02,-0.95902765E-02,-0.80086133E-02,-0.13201639E-01,
     + -0.15754739E-03,-0.16542245E-02, 0.20042380E-01,-0.22898220E-01,
     + -0.37723646E-03,-0.67448226E-03,-0.74074074E-03,-0.43797810E-03,
     + -0.16468442E-02,-0.73253858E-03,-0.41120667E-02,-0.29294782E-02,
     + -0.84081729E-03,-0.51961206E-02,-0.39537963E-02, 0.17203884E-02,
     + -0.56993105E-02,-0.36271384E-02,-0.21451726E-02,-0.56850831E-02,
     + -0.13560791E-01,-0.59638922E-02, 0.13953323E-01,-0.59757028E-02,
     + -0.12389302E-01, 0.28127786E-02, 0.10901931E-01,-0.11722708E-01,
     + -0.10234840E-01, 0.10963762E-01,-0.93124479E-04,-0.14342860E-02,
     +  0.24927477E-03, 0.12318705E-02,-0.12944845E-02,-0.15610989E-03,
     + -0.91408135E-03,-0.21442014E-02,-0.75004355E-03, 0.98161458E-03,
     + -0.64473943E-03,-0.84925565E-03,-0.31396050E-02,-0.25792527E-02,
     +  0.62562461E-03,-0.48036473E-02, 0.52699964E-02, 0.12530525E-02,
     + -0.17014405E-02, 0.21924409E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x_sr6_largex0_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sr6_largex0_q3en   =x_sr6_largex0_q3en   
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x23*x31*x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)    *x21    *x41*x51
      x_sr6_largex0_q3en   =x_sr6_largex0_q3en   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x23*x31*x42    
     2  +coeff( 20)    *x23    *x43    
     3  +coeff( 21)        *x31*x41    
     4  +coeff( 22)            *x41*x51
     5  +coeff( 23)*x11        *x41    
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)    *x21*x31    *x51
     8  +coeff( 26)    *x23*x31        
      x_sr6_largex0_q3en   =x_sr6_largex0_q3en   
     9  +coeff( 27)    *x23    *x41    
     1  +coeff( 28)    *x21*x32*x41    
     2  +coeff( 29)    *x21*x31*x42    
     3  +coeff( 30)*x11    *x33        
     4  +coeff( 31)*x11*x22    *x42    
     5  +coeff( 32)        *x31        
     6  +coeff( 33)*x11            *x51
     7  +coeff( 34)*x11    *x31*x41    
     8  +coeff( 35)*x11        *x42    
      x_sr6_largex0_q3en   =x_sr6_largex0_q3en   
     9  +coeff( 36)    *x21*x33        
     1  +coeff( 37)    *x21    *x42*x51
     2  +coeff( 38)*x11    *x31*x42    
     3  +coeff( 39)*x11        *x43    
     4  +coeff( 40)    *x21*x31*x43    
     5  +coeff( 41)*x11*x22*x32        
     6  +coeff( 42)    *x21*x31*x41*x53
     7  +coeff( 43)    *x23    *x43*x52
     8  +coeff( 44)*x12*x21    *x43*x51
      x_sr6_largex0_q3en   =x_sr6_largex0_q3en   
     9  +coeff( 45)        *x32        
     1  +coeff( 46)*x11    *x31        
     2  +coeff( 47)    *x22*x31        
     3  +coeff( 48)    *x22    *x41    
     4  +coeff( 49)    *x24            
     5  +coeff( 50)            *x42*x51
     6  +coeff( 51)    *x22*x31*x41    
     7  +coeff( 52)    *x22    *x42    
     8  +coeff( 53)*x11        *x41*x51
      x_sr6_largex0_q3en   =x_sr6_largex0_q3en   
     9  +coeff( 54)    *x23    *x41*x51
     1  +coeff( 55)*x11    *x32*x41    
     2  +coeff( 56)    *x22        *x53
     3  +coeff( 57)    *x21*x32*x42    
     4  +coeff( 58)    *x21*x31*x42*x51
     5  +coeff( 59)    *x21*x32    *x52
     6  +coeff( 60)    *x24    *x42    
     7  +coeff( 61)    *x23    *x42*x51
     8  +coeff( 62)    *x22    *x43*x51
      x_sr6_largex0_q3en   =x_sr6_largex0_q3en   
     9  +coeff( 63)    *x21*x32*x43    
     1  +coeff( 64)    *x21*x32*x41*x52
     2  +coeff( 65)*x12*x21    *x42*x51
     3  +coeff( 66)*x12*x21        *x53
     4  +coeff( 67)*x12*x21*x32*x41*x51
     5  +coeff( 68)*x12*x23*x31*x41*x51
     6  +coeff( 69)*x12*x22*x31*x42*x51
     7  +coeff( 70)*x12*x23    *x41*x53
     8  +coeff( 71)        *x31    *x51
      x_sr6_largex0_q3en   =x_sr6_largex0_q3en   
     9  +coeff( 72)*x11*x22            
     1  +coeff( 73)*x11*x21        *x51
     2  +coeff( 74)    *x23        *x51
     3  +coeff( 75)*x11    *x32        
     4  +coeff( 76)*x11    *x31    *x51
     5  +coeff( 77)    *x21*x32    *x51
     6  +coeff( 78)    *x21*x31*x41*x51
     7  +coeff( 79)    *x21        *x53
     8  +coeff( 80)        *x31*x43    
      x_sr6_largex0_q3en   =x_sr6_largex0_q3en   
     9  +coeff( 81)*x11*x22        *x51
     1  +coeff( 82)            *x43*x51
     2  +coeff( 83)    *x22*x31*x42    
     3  +coeff( 84)    *x22    *x43    
     4  +coeff( 85)    *x22*x32    *x51
     5  +coeff( 86)    *x21*x32*x41*x51
     6  +coeff( 87)    *x21    *x43*x51
     7  +coeff( 88)    *x21    *x42*x52
     8  +coeff( 89)    *x21    *x41*x53
      x_sr6_largex0_q3en   =x_sr6_largex0_q3en   
     9  +coeff( 90)*x11*x22*x31*x41    
c
      return
      end
      function t_sr6_largex0_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4344000E-02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28511882E-02,-0.77838935E-01, 0.62561601E-04, 0.22797789E-01,
     +  0.70216111E-02,-0.47328542E-02, 0.57811692E-03,-0.51133465E-02,
     +  0.53352944E-03,-0.14012600E-02, 0.13059382E-02, 0.15939178E-01,
     + -0.95408701E-03,-0.62825589E-03, 0.11083577E-01, 0.13504189E-01,
     + -0.38246240E-03,-0.18645298E-05, 0.61397878E-02, 0.87220798E-03,
     +  0.15447945E-02,-0.46709854E-04,-0.14406512E-02, 0.65555591E-02,
     +  0.17782703E-01, 0.45025135E-04, 0.23551739E-02, 0.45868580E-03,
     + -0.12686144E-03, 0.21939282E-04,-0.60629864E-04, 0.37838906E-03,
     +  0.43934368E-03, 0.24254075E-02, 0.59727230E-03, 0.37535085E-03,
     +  0.12159490E-02, 0.82126912E-03, 0.23467161E-02, 0.36094899E-02,
     +  0.58968540E-03,-0.41576524E-03, 0.19921949E-02, 0.20066353E-02,
     + -0.50095161E-02, 0.14121996E-03, 0.46283644E-03,-0.39293812E-03,
     +  0.26036685E-03, 0.23309690E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      t_sr6_largex0_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_sr6_largex0_q3en   =t_sr6_largex0_q3en   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x23*x31*x41    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x23*x32        
      t_sr6_largex0_q3en   =t_sr6_largex0_q3en   
     9  +coeff( 18)    *x23*x31*x42    
     1  +coeff( 19)    *x23    *x43    
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)*x11        *x42    
     4  +coeff( 22)    *x23*x31        
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x21*x32*x41    
     7  +coeff( 25)    *x21*x31*x42    
     8  +coeff( 26)*x11    *x33        
      t_sr6_largex0_q3en   =t_sr6_largex0_q3en   
     9  +coeff( 27)    *x23    *x41*x51
     1  +coeff( 28)    *x21*x32    *x52
     2  +coeff( 29)*x11    *x33*x41    
     3  +coeff( 30)        *x31        
     4  +coeff( 31)            *x42    
     5  +coeff( 32)    *x23            
     6  +coeff( 33)    *x22*x31        
     7  +coeff( 34)    *x21*x32        
     8  +coeff( 35)    *x22    *x41    
      t_sr6_largex0_q3en   =t_sr6_largex0_q3en   
     9  +coeff( 36)*x11*x22            
     1  +coeff( 37)*x11    *x31*x41    
     2  +coeff( 38)    *x21*x33        
     3  +coeff( 39)    *x22    *x42    
     4  +coeff( 40)    *x21*x31*x43    
     5  +coeff( 41)    *x23*x31    *x51
     6  +coeff( 42)        *x31*x43*x51
     7  +coeff( 43)    *x22*x33*x41    
     8  +coeff( 44)    *x22    *x43*x51
      t_sr6_largex0_q3en   =t_sr6_largex0_q3en   
     9  +coeff( 45)    *x23    *x43*x52
     1  +coeff( 46)        *x32        
     2  +coeff( 47)*x11    *x31        
     3  +coeff( 48)    *x22        *x51
     4  +coeff( 49)    *x21        *x52
     5  +coeff( 50)*x11    *x32        
c
      return
      end
      function y_sr6_largex0_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1702326E-02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24162434E-03,-0.48349775E-01,-0.15530429E-02, 0.20204861E-01,
     + -0.13477546E-02,-0.62501221E-03, 0.41883150E-02, 0.10932704E-01,
     +  0.69099389E-01, 0.30279200E-03, 0.73792241E-02,-0.51905788E-02,
     + -0.62856339E-01,-0.88381944E-02,-0.25257664E-02,-0.15585209E-02,
     + -0.32913878E-02,-0.54694656E-02, 0.66621671E-02,-0.48006307E-02,
     +  0.61125942E-02, 0.69809328E-02, 0.10534674E-02,-0.20880776E-03,
     +  0.20856797E-02,-0.13472025E-01,-0.13318643E-01,-0.17939780E-02,
     + -0.14365755E-02,-0.18811983E-02,-0.44274807E-01, 0.12981900E-02,
     + -0.86223921E-02, 0.20415092E-01, 0.23187576E-01, 0.20979729E-02,
     +  0.43152389E-02,-0.44327378E-02,-0.20553246E-01,-0.34614597E-02,
     + -0.37537257E-02,-0.27091989E-01,-0.35341837E-01,-0.82213432E-03,
     + -0.24706260E-02,-0.12839074E-02,-0.76766638E-02, 0.22304475E-02,
     +  0.13341233E-01,-0.15684763E-01,-0.57831379E-02,-0.54206266E-02,
     + -0.16203254E-01,-0.44197123E-01,-0.70564882E-02, 0.12490462E-02,
     +  0.14127421E-02,-0.16429600E-02, 0.50321463E-02,-0.55940274E-01,
     + -0.16274080E-01,-0.17495865E-01, 0.10729288E-03, 0.27257735E-02,
     + -0.12074362E-02, 0.86726749E-03, 0.23720644E-02,-0.77617904E-02,
     + -0.90553248E-02, 0.47456664E-02, 0.37498481E-02,-0.20591594E-01,
     + -0.61593000E-02,-0.41573821E-02, 0.17675130E-01, 0.14777996E-01,
     +  0.33491314E-02,-0.13632791E-02,-0.32194424E-02,-0.11762292E-01,
     + -0.28086011E-02,-0.54732938E-02,-0.26178164E-02, 0.17858846E-02,
     + -0.12771319E-01,-0.46404939E-01, 0.32729451E-01, 0.94587179E-02,
     + -0.20141836E-01,-0.23384126E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sr6_largex0_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x31    *x51
      y_sr6_largex0_q3en   =y_sr6_largex0_q3en   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11                
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x21        *x51
     6  +coeff( 15)        *x32*x41    
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)            *x42*x51
      y_sr6_largex0_q3en   =y_sr6_largex0_q3en   
     9  +coeff( 18)    *x22            
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)    *x21*x31    *x51
     7  +coeff( 25)        *x33*x41    
     8  +coeff( 26)        *x31*x43    
      y_sr6_largex0_q3en   =y_sr6_largex0_q3en   
     9  +coeff( 27)            *x44    
     1  +coeff( 28)*x11*x21            
     2  +coeff( 29)        *x32*x41*x51
     3  +coeff( 30)            *x43*x51
     4  +coeff( 31)    *x22    *x41    
     5  +coeff( 32)        *x31*x41*x52
     6  +coeff( 33)    *x22        *x51
     7  +coeff( 34)    *x21*x31*x42    
     8  +coeff( 35)    *x21    *x43    
      y_sr6_largex0_q3en   =y_sr6_largex0_q3en   
     9  +coeff( 36)    *x21*x31*x41*x51
     1  +coeff( 37)    *x23            
     2  +coeff( 38)        *x32*x43    
     3  +coeff( 39)            *x45    
     4  +coeff( 40)*x11*x21*x31        
     5  +coeff( 41)*x11*x21    *x41    
     6  +coeff( 42)    *x22*x31*x41    
     7  +coeff( 43)    *x22    *x42    
     8  +coeff( 44)*x11*x21        *x51
      y_sr6_largex0_q3en   =y_sr6_largex0_q3en   
     9  +coeff( 45)        *x31*x42*x52
     1  +coeff( 46)            *x43*x52
     2  +coeff( 47)    *x22    *x41*x51
     3  +coeff( 48)    *x21*x31*x42*x51
     4  +coeff( 49)    *x23    *x41    
     5  +coeff( 50)*x11*x21    *x42    
     6  +coeff( 51)        *x31*x44*x51
     7  +coeff( 52)    *x22*x33        
     8  +coeff( 53)    *x22*x31*x42    
      y_sr6_largex0_q3en   =y_sr6_largex0_q3en   
     9  +coeff( 54)    *x22    *x43    
     1  +coeff( 55)*x11*x21    *x41*x51
     2  +coeff( 56)    *x21    *x41*x53
     3  +coeff( 57)    *x22*x34        
     4  +coeff( 58)*x11*x21*x33*x41    
     5  +coeff( 59)    *x22*x35        
     6  +coeff( 60)    *x22    *x44*x51
     7  +coeff( 61)*x11*x21*x31*x41*x53
     8  +coeff( 62)*x12*x22    *x42*x51
      y_sr6_largex0_q3en   =y_sr6_largex0_q3en   
     9  +coeff( 63)        *x33        
     1  +coeff( 64)            *x43    
     2  +coeff( 65)        *x32    *x51
     3  +coeff( 66)*x11    *x31        
     4  +coeff( 67)    *x21*x32        
     5  +coeff( 68)    *x21    *x41*x51
     6  +coeff( 69)    *x22*x31        
     7  +coeff( 70)    *x21*x32*x41    
     8  +coeff( 71)    *x21    *x42*x51
      y_sr6_largex0_q3en   =y_sr6_largex0_q3en   
     9  +coeff( 72)        *x31*x44    
     1  +coeff( 73)    *x22*x32        
     2  +coeff( 74)    *x22*x31    *x51
     3  +coeff( 75)    *x21*x31*x43    
     4  +coeff( 76)    *x21    *x44    
     5  +coeff( 77)    *x23*x31        
     6  +coeff( 78)*x12        *x41    
     7  +coeff( 79)*x11*x21*x32        
     8  +coeff( 80)*x11*x21*x31*x41    
      y_sr6_largex0_q3en   =y_sr6_largex0_q3en   
     9  +coeff( 81)*x11*x21*x31    *x51
     1  +coeff( 82)    *x22*x31*x41*x51
     2  +coeff( 83)*x12        *x42    
     3  +coeff( 84)*x11        *x42*x52
     4  +coeff( 85)*x11*x21    *x42*x51
     5  +coeff( 86)    *x22*x31*x43*x51
     6  +coeff( 87)    *x22*x32*x45    
     7  +coeff( 88)    *x22    *x42*x54
     8  +coeff( 89)*x11*x23    *x45    
      y_sr6_largex0_q3en   =y_sr6_largex0_q3en   
     9  +coeff( 90)*x12*x22    *x43*x51
c
      return
      end
      function p_sr6_largex0_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1273859E-02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.90547773E-03,-0.69291746E-04,-0.11016385E-01,-0.12451623E-01,
     +  0.22725975E-02,-0.71607542E-03,-0.12785061E-02,-0.16785448E-03,
     + -0.12908679E-01,-0.33967631E-03,-0.15260035E-02, 0.15901611E-02,
     +  0.13004278E-01, 0.16453606E-02, 0.13767505E-02, 0.39886902E-03,
     + -0.93148528E-02, 0.10008960E-02,-0.18125726E-02,-0.32659783E-02,
     +  0.22044824E-02,-0.34157097E-03,-0.33290437E-03, 0.29631596E-03,
     + -0.13673033E-02, 0.66052208E-03,-0.21053883E-02,-0.11798490E-02,
     + -0.38829402E-03, 0.28449206E-02,-0.12312217E-02, 0.31409105E-02,
     + -0.39228406E-02,-0.19355001E-02, 0.22495997E-02, 0.23908784E-04,
     +  0.15978343E-03,-0.31908028E-03,-0.13878644E-02,-0.53976156E-03,
     + -0.31256955E-03,-0.30264512E-03,-0.19181634E-02,-0.25987842E-02,
     +  0.24731371E-02,-0.10982837E-02, 0.64577546E-03,-0.45846710E-02,
     +  0.12648742E-02, 0.39367704E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      p_sr6_largex0_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sr6_largex0_q3en   =p_sr6_largex0_q3en   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22    *x41    
      p_sr6_largex0_q3en   =p_sr6_largex0_q3en   
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)    *x22    *x41*x51
     3  +coeff( 21)    *x23*x31*x41    
     4  +coeff( 22)        *x31*x41    
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)*x11    *x31        
     7  +coeff( 25)    *x22*x31        
     8  +coeff( 26)    *x21*x31*x41    
      p_sr6_largex0_q3en   =p_sr6_largex0_q3en   
     9  +coeff( 27)            *x43    
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)*x11*x21        *x51
     3  +coeff( 30)    *x23    *x41    
     4  +coeff( 31)    *x22*x31    *x51
     5  +coeff( 32)    *x23    *x42    
     6  +coeff( 33)    *x22    *x42*x51
     7  +coeff( 34)*x11*x23    *x41*x51
     8  +coeff( 35)    *x22*x33*x41*x51
      p_sr6_largex0_q3en   =p_sr6_largex0_q3en   
     9  +coeff( 36)*x11                
     1  +coeff( 37)*x11            *x51
     2  +coeff( 38)        *x32*x41    
     3  +coeff( 39)        *x31*x42    
     4  +coeff( 40)            *x42*x51
     5  +coeff( 41)            *x41*x52
     6  +coeff( 42)*x11*x21*x31        
     7  +coeff( 43)    *x22*x31*x41    
     8  +coeff( 44)    *x22    *x42    
      p_sr6_largex0_q3en   =p_sr6_largex0_q3en   
     9  +coeff( 45)    *x21    *x43    
     1  +coeff( 46)*x11*x21    *x42    
     2  +coeff( 47)    *x23*x32        
     3  +coeff( 48)    *x22*x31*x41*x51
     4  +coeff( 49)*x11*x22    *x42    
     5  +coeff( 50)    *x23*x31*x42    
c
      return
      end
      function l_sr6_largex0_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1837640E+02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.49938E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49950E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12366776E-01, 0.31348497E+00, 0.32125436E-01,-0.33347465E-01,
     +  0.29279580E-01,-0.12822441E-02, 0.19554211E-01,-0.67278281E-01,
     +  0.32087734E-02,-0.21350035E-02,-0.54087881E-02,-0.44928227E-01,
     +  0.52516791E-02, 0.75232810E-02, 0.31361324E-02,-0.13751284E-02,
     + -0.83743911E-02,-0.77976860E-01,-0.59913326E-01,-0.10235506E-01,
     +  0.22314345E-02,-0.23085461E-02,-0.18039632E-02, 0.26119237E-02,
     +  0.90772435E-02,-0.30530868E-01,-0.31306973E-03,-0.34063573E-02,
     +  0.19256216E-02, 0.17678426E-02,-0.12036143E-02, 0.22477801E-02,
     + -0.73154648E-02,-0.88541377E-02,-0.46537230E-02,-0.12938913E-01,
     +  0.59062871E-02,-0.18047415E-01,-0.27571134E-02,-0.27591013E-02,
     +  0.14034587E-03,-0.44021931E-01, 0.11168109E-01,-0.15326124E-01,
     +  0.23082249E-01,-0.82401233E-03,-0.15054840E-02, 0.45806140E-03,
     + -0.18049411E-03,-0.50139963E-03,-0.82328841E-02, 0.85052643E-05,
     + -0.22130362E-02, 0.27826289E-02,-0.87111499E-02,-0.81672473E-02,
     + -0.76932991E-02,-0.11880982E-01,-0.12425012E-01,-0.19982248E-02,
     + -0.16707641E-02, 0.16467929E-01,-0.13046876E-01,-0.14325153E-01,
     + -0.11399435E-01,-0.71114218E-02,-0.83127264E-02, 0.34002634E-02,
     + -0.18420106E-01, 0.17929578E-01,-0.14821894E-01, 0.24573328E-02,
     +  0.22518118E-02,-0.66702180E-02, 0.58622409E-01, 0.14095540E-01,
     + -0.88500129E-02, 0.13354729E-01, 0.28467720E-03, 0.26278550E-03,
     +  0.22780601E-03,-0.10235077E-02,-0.98051474E-04,-0.17908163E-02,
     + -0.19757808E-02,-0.13603548E-02,-0.83198521E-03,-0.13459023E-02,
     + -0.37028757E-03, 0.11850456E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sr6_largex0_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x21    *x42    
      l_sr6_largex0_q3en   =l_sr6_largex0_q3en   
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)            *x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)    *x21*x32        
      l_sr6_largex0_q3en   =l_sr6_largex0_q3en   
     9  +coeff( 18)    *x21*x31*x42    
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)    *x23*x32*x41    
     3  +coeff( 21)        *x31        
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x22    *x41    
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)    *x23    *x41    
     8  +coeff( 26)    *x21*x32*x41    
      l_sr6_largex0_q3en   =l_sr6_largex0_q3en   
     9  +coeff( 27)*x11*x22*x31        
     1  +coeff( 28)*x11*x22    *x42    
     2  +coeff( 29)        *x31*x41    
     3  +coeff( 30)            *x41*x51
     4  +coeff( 31)*x11    *x31        
     5  +coeff( 32)    *x22        *x51
     6  +coeff( 33)*x11    *x31*x41    
     7  +coeff( 34)*x11        *x42    
     8  +coeff( 35)    *x21*x33        
      l_sr6_largex0_q3en   =l_sr6_largex0_q3en   
     9  +coeff( 36)    *x22    *x42    
     1  +coeff( 37)    *x21    *x42*x51
     2  +coeff( 38)    *x21*x31*x43    
     3  +coeff( 39)    *x22    *x42*x51
     4  +coeff( 40)*x11*x22*x32        
     5  +coeff( 41)    *x22*x33*x41    
     6  +coeff( 42)    *x23    *x43    
     7  +coeff( 43)    *x21*x32*x43    
     8  +coeff( 44)    *x22    *x43*x51
      l_sr6_largex0_q3en   =l_sr6_largex0_q3en   
     9  +coeff( 45)*x12*x23    *x41*x51
     1  +coeff( 46)    *x21*x31    *x51
     2  +coeff( 47)            *x42*x51
     3  +coeff( 48)*x11*x21        *x51
     4  +coeff( 49)*x11        *x41*x51
     5  +coeff( 50)*x11            *x52
     6  +coeff( 51)    *x22*x31*x41    
     7  +coeff( 52)        *x31*x43    
     8  +coeff( 53)    *x21*x31*x41*x51
      l_sr6_largex0_q3en   =l_sr6_largex0_q3en   
     9  +coeff( 54)*x11*x22    *x41    
     1  +coeff( 55)*x11    *x31*x42    
     2  +coeff( 56)*x11        *x43    
     3  +coeff( 57)    *x21*x32*x42    
     4  +coeff( 58)    *x22    *x43    
     5  +coeff( 59)    *x23    *x41*x51
     6  +coeff( 60)    *x21*x32    *x52
     7  +coeff( 61)        *x31*x41*x53
     8  +coeff( 62)    *x21*x33*x42    
      l_sr6_largex0_q3en   =l_sr6_largex0_q3en   
     9  +coeff( 63)*x12*x21    *x41*x51
     1  +coeff( 64)    *x23    *x42*x51
     2  +coeff( 65)    *x22*x31*x42*x51
     3  +coeff( 66)*x11*x22*x32*x41    
     4  +coeff( 67)    *x22*x33*x42    
     5  +coeff( 68)*x12*x21    *x43    
     6  +coeff( 69)*x12*x21    *x42*x51
     7  +coeff( 70)    *x23    *x42*x52
     8  +coeff( 71)    *x22*x31*x42*x52
      l_sr6_largex0_q3en   =l_sr6_largex0_q3en   
     9  +coeff( 72)        *x33*x42*x52
     1  +coeff( 73)*x12*x21        *x53
     2  +coeff( 74)*x11*x22    *x43*x51
     3  +coeff( 75)    *x23*x32*x43    
     4  +coeff( 76)*x11*x23    *x43*x51
     5  +coeff( 77)*x11*x21*x32*x43*x51
     6  +coeff( 78)*x12*x21*x32*x42*x51
     7  +coeff( 79)        *x31    *x51
     8  +coeff( 80)                *x52
      l_sr6_largex0_q3en   =l_sr6_largex0_q3en   
     9  +coeff( 81)*x12                
     1  +coeff( 82)    *x23            
     2  +coeff( 83)        *x33        
     3  +coeff( 84)        *x32*x41    
     4  +coeff( 85)        *x31*x42    
     5  +coeff( 86)    *x21        *x52
     6  +coeff( 87)*x11    *x32        
     7  +coeff( 88)    *x22*x32        
     8  +coeff( 89)*x12        *x41    
      l_sr6_largex0_q3en   =l_sr6_largex0_q3en   
     9  +coeff( 90)    *x23        *x51
c
      return
      end
      function x_sr6_largex0_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.9637365E-03/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.45980E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49875E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12937909E-01, 0.71704797E-01,-0.59818882E-04,-0.18494602E-03,
     +  0.31686243E+00,-0.25687426E-01, 0.34983002E-01,-0.20733370E-01,
     +  0.15529077E-01,-0.34979027E-01,-0.33881710E-03,-0.60353028E-02,
     + -0.19462965E-01, 0.44230785E-03,-0.14281810E-02,-0.43797041E-02,
     + -0.29234649E-02,-0.91602997E-03,-0.17344502E-02,-0.11364687E-02,
     +  0.20108044E-02,-0.68559748E-03,-0.44016093E-02, 0.19437429E-02,
     +  0.20689489E-02,-0.28301012E-01,-0.25326142E-01,-0.42795613E-02,
     +  0.66257623E-03,-0.28855240E-03,-0.97833364E-03,-0.57249225E-03,
     + -0.26071693E-02,-0.13380541E-01,-0.19283625E-02, 0.48000962E-03,
     +  0.23559135E-03,-0.68614264E-02,-0.45620068E-03,-0.91873697E-03,
     + -0.39745215E-02,-0.49316003E-02, 0.17919638E-02,-0.93169993E-03,
     +  0.53051306E-03,-0.10855785E-04, 0.42503182E-03,-0.22986964E-02,
     + -0.18043621E-02, 0.25260684E-04,-0.95588865E-03,-0.70723257E-03,
     + -0.17127151E-02,-0.28779791E-03,-0.18552060E-02,-0.61359545E-02,
     +  0.78472253E-02,-0.50472431E-02, 0.97963205E-02,-0.38144980E-02,
     + -0.62986570E-02, 0.14384357E-02, 0.64006648E-02, 0.54525239E-02,
     + -0.73763398E-02,-0.22165298E-01,-0.13976428E-01, 0.77317436E-02,
     + -0.19462207E-01,-0.37768279E-03, 0.42970787E-03,-0.17654037E-03,
     + -0.86481974E-03,-0.23045790E-03,-0.85430342E-03,-0.23812929E-02,
     +  0.31437675E-03, 0.42968697E-03,-0.17395512E-02,-0.69906493E-03,
     + -0.15549511E-02,-0.50558615E-03,-0.12119086E-02,-0.72402885E-03,
     + -0.59216213E-03,-0.13005923E-02, 0.10078871E-02,-0.24298739E-02,
     + -0.30162702E-02,-0.38788983E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x_sr6_largex0_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_sr6_largex0_q3ex   =x_sr6_largex0_q3ex   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)    *x21    *x41*x51
      x_sr6_largex0_q3ex   =x_sr6_largex0_q3ex   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)            *x41*x51
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)    *x21*x31    *x51
     5  +coeff( 23)    *x24            
     6  +coeff( 24)                *x53
     7  +coeff( 25)    *x23    *x41    
     8  +coeff( 26)    *x21*x31*x42    
      x_sr6_largex0_q3ex   =x_sr6_largex0_q3ex   
     9  +coeff( 27)    *x21    *x43    
     1  +coeff( 28)*x11*x22    *x42    
     2  +coeff( 29)    *x23*x32*x41    
     3  +coeff( 30)        *x31    *x51
     4  +coeff( 31)*x11        *x41    
     5  +coeff( 32)*x11            *x51
     6  +coeff( 33)*x11    *x31*x41    
     7  +coeff( 34)    *x21*x32*x41    
     8  +coeff( 35)    *x21    *x42*x51
      x_sr6_largex0_q3ex   =x_sr6_largex0_q3ex   
     9  +coeff( 36)*x11*x22*x31        
     1  +coeff( 37)*x11    *x33        
     2  +coeff( 38)    *x21*x31*x43    
     3  +coeff( 39)*x11*x22*x31*x41    
     4  +coeff( 40)    *x23*x33        
     5  +coeff( 41)    *x23*x31*x41*x51
     6  +coeff( 42)    *x23    *x42*x51
     7  +coeff( 43)*x11    *x31*x43    
     8  +coeff( 44)*x11    *x31        
      x_sr6_largex0_q3ex   =x_sr6_largex0_q3ex   
     9  +coeff( 45)    *x22    *x41    
     1  +coeff( 46)*x12                
     2  +coeff( 47)*x11*x21    *x41    
     3  +coeff( 48)    *x22*x31*x41    
     4  +coeff( 49)*x11        *x42    
     5  +coeff( 50)*x11        *x41*x51
     6  +coeff( 51)    *x22        *x52
     7  +coeff( 52)    *x21*x32    *x51
     8  +coeff( 53)            *x43*x51
      x_sr6_largex0_q3ex   =x_sr6_largex0_q3ex   
     9  +coeff( 54)*x11*x24            
     1  +coeff( 55)*x11*x22*x32        
     2  +coeff( 56)    *x23    *x43    
     3  +coeff( 57)    *x21*x32*x43    
     4  +coeff( 58)*x11*x22*x31*x42    
     5  +coeff( 59)    *x23    *x42*x52
     6  +coeff( 60)    *x22*x31*x41*x53
     7  +coeff( 61)*x12*x21    *x42*x51
     8  +coeff( 62)*x12*x21        *x53
      x_sr6_largex0_q3ex   =x_sr6_largex0_q3ex   
     9  +coeff( 63)*x11*x22*x32*x42    
     1  +coeff( 64)    *x24*x31*x43    
     2  +coeff( 65)*x12*x21    *x43*x51
     3  +coeff( 66)    *x23*x32*x42*x52
     4  +coeff( 67)*x11*x22    *x43*x53
     5  +coeff( 68)    *x23*x33*x41*x53
     6  +coeff( 69)*x11*x23*x31*x43*x53
     7  +coeff( 70)            *x42*x51
     8  +coeff( 71)            *x41*x52
      x_sr6_largex0_q3ex   =x_sr6_largex0_q3ex   
     9  +coeff( 72)*x11*x21        *x51
     1  +coeff( 73)    *x23        *x51
     2  +coeff( 74)*x11    *x31    *x51
     3  +coeff( 75)    *x21*x33        
     4  +coeff( 76)    *x21*x31*x41*x51
     5  +coeff( 77)*x11*x22    *x41    
     6  +coeff( 78)        *x31*x43    
     7  +coeff( 79)    *x24        *x51
     8  +coeff( 80)        *x31*x42*x51
      x_sr6_largex0_q3ex   =x_sr6_largex0_q3ex   
     9  +coeff( 81)    *x23*x31*x41    
     1  +coeff( 82)*x11*x21    *x42    
     2  +coeff( 83)    *x23    *x41*x51
     3  +coeff( 84)    *x23        *x52
     4  +coeff( 85)*x12*x22            
     5  +coeff( 86)*x11        *x43    
     6  +coeff( 87)    *x22        *x53
     7  +coeff( 88)    *x21*x32*x41*x51
     8  +coeff( 89)    *x21*x31*x42*x51
      x_sr6_largex0_q3ex   =x_sr6_largex0_q3ex   
     9  +coeff( 90)*x12        *x42    
c
      return
      end
      function t_sr6_largex0_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2203062E-02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.45980E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49875E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.61585116E-02,-0.23162574E-01, 0.15732811E-05, 0.28237160E-04,
     +  0.10730727E+00,-0.36552469E-02, 0.58145951E-02,-0.97912531E-02,
     +  0.79246244E-03,-0.19698632E-02, 0.44975476E-03,-0.50927163E-03,
     + -0.13155497E-02,-0.51920931E-03,-0.58997929E-03, 0.69754827E-03,
     +  0.35917087E-03, 0.19878328E-02, 0.48885533E-04, 0.65138222E-04,
     + -0.55805122E-03,-0.21212600E-03, 0.20446301E-03, 0.30250934E-04,
     +  0.11258817E-02, 0.60196151E-03,-0.10624320E-03,-0.90086098E-04,
     + -0.80376289E-04, 0.37388419E-03,-0.21546854E-03, 0.13985607E-03,
     +  0.59977308E-03, 0.87864295E-03, 0.32819898E-03, 0.92666472E-04,
     + -0.40988045E-03,-0.27886956E-03,-0.76370960E-03, 0.63249748E-03,
     +  0.16086281E-03, 0.10083892E-03, 0.66280598E-04, 0.57817695E-04,
     + -0.19657504E-03, 0.20416020E-03, 0.19193266E-03, 0.51271496E-03,
     + -0.39786176E-03,-0.36124481E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sr6_largex0_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      t_sr6_largex0_q3ex   =t_sr6_largex0_q3ex   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)                *x53
     8  +coeff( 17)    *x23    *x41    
      t_sr6_largex0_q3ex   =t_sr6_largex0_q3ex   
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x21*x31        
     2  +coeff( 20)        *x32        
     3  +coeff( 21)        *x31*x41    
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)    *x23            
     6  +coeff( 24)            *x42*x51
     7  +coeff( 25)    *x21    *x43    
     8  +coeff( 26)*x11*x23            
      t_sr6_largex0_q3ex   =t_sr6_largex0_q3ex   
     9  +coeff( 27)        *x31    *x51
     1  +coeff( 28)*x12                
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)    *x22    *x41    
     4  +coeff( 31)        *x31*x41*x51
     5  +coeff( 32)*x11            *x52
     6  +coeff( 33)    *x22*x31*x41    
     7  +coeff( 34)    *x21*x31*x42    
     8  +coeff( 35)    *x22    *x41*x51
      t_sr6_largex0_q3ex   =t_sr6_largex0_q3ex   
     9  +coeff( 36)    *x22*x33        
     1  +coeff( 37)    *x22*x31*x41*x51
     2  +coeff( 38)        *x32    *x53
     3  +coeff( 39)    *x23    *x42*x51
     4  +coeff( 40)*x11*x22    *x43    
     5  +coeff( 41)    *x22*x31        
     6  +coeff( 42)*x11        *x42    
     7  +coeff( 43)*x11        *x41*x51
     8  +coeff( 44)*x12            *x51
      t_sr6_largex0_q3ex   =t_sr6_largex0_q3ex   
     9  +coeff( 45)    *x21    *x41*x52
     1  +coeff( 46)        *x31*x41*x52
     2  +coeff( 47)*x11*x22        *x51
     3  +coeff( 48)    *x23    *x42    
     4  +coeff( 49)    *x21*x32*x42    
     5  +coeff( 50)    *x21*x32*x41*x51
c
      return
      end
      function y_sr6_largex0_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3718517E-02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.45980E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49875E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.19449459E-02,-0.42023785E-01,-0.24543574E-01, 0.12667474E-01,
     + -0.12880679E-02, 0.12963012E-03, 0.14798305E-03, 0.10544024E-02,
     +  0.57857637E-02, 0.52077249E-01, 0.24844534E-03, 0.62973825E-02,
     + -0.53759632E-02,-0.50852571E-01,-0.63127768E-02,-0.69049033E-02,
     + -0.10017086E-01,-0.38264741E-03,-0.18006887E-02,-0.39218669E-02,
     +  0.54815141E-02, 0.63995778E-03, 0.84728701E-02, 0.10149334E-01,
     +  0.74504584E-03, 0.31768813E-03,-0.45271711E-02,-0.38689016E-02,
     + -0.16954064E-02,-0.67361241E-03,-0.19868435E-02,-0.18937660E-02,
     + -0.66049141E-02,-0.40778529E-01,-0.71913977E-02, 0.32349417E-02,
     +  0.10054708E-02,-0.79778982E-02,-0.88945543E-03,-0.13055117E-01,
     +  0.10136641E-02,-0.60421909E-03, 0.64130374E-02,-0.96438834E-02,
     + -0.45480891E-02, 0.61060407E-03,-0.32792785E-02, 0.70346817E-02,
     + -0.15867158E-02, 0.71397307E-02, 0.39714076E-01,-0.14802313E+00,
     +  0.45389980E-01, 0.11995018E+00, 0.76404045E-03, 0.17954319E-02,
     + -0.70236400E-02,-0.70423088E-02, 0.13694108E-01, 0.10024659E-01,
     + -0.19204007E-02,-0.15986833E-02,-0.40000957E-02,-0.13747687E-01,
     + -0.42688060E-02, 0.24001990E-02,-0.58346610E-02,-0.18526529E-02,
     + -0.11577667E-01,-0.10345009E-01, 0.82506420E-04,-0.32546073E-01,
     + -0.79878410E-02,-0.85510965E-02, 0.16816428E-01,-0.31425962E-02,
     + -0.33090841E-01,-0.52505080E-01,-0.92255753E-02,-0.31705366E-02,
     + -0.19413701E+00, 0.58775835E-01,-0.97557321E-01,-0.23012672E+00,
     + -0.56158509E-01,-0.32310176E-02, 0.51027453E-02, 0.29040752E-02,
     +  0.29879063E-02,-0.24982600E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sr6_largex0_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sr6_largex0_q3ex   =y_sr6_largex0_q3ex   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)        *x31*x42    
     8  +coeff( 17)            *x43    
      y_sr6_largex0_q3ex   =y_sr6_largex0_q3ex   
     9  +coeff( 18)        *x31*x41*x51
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)    *x22            
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)        *x31    *x52
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x21*x31    *x51
      y_sr6_largex0_q3ex   =y_sr6_largex0_q3ex   
     9  +coeff( 27)        *x32*x42    
     1  +coeff( 28)            *x44    
     2  +coeff( 29)*x11*x21            
     3  +coeff( 30)        *x32*x41*x51
     4  +coeff( 31)        *x31*x42*x51
     5  +coeff( 32)            *x43*x51
     6  +coeff( 33)    *x22*x31        
     7  +coeff( 34)    *x22    *x41    
     8  +coeff( 35)    *x22        *x51
      y_sr6_largex0_q3ex   =y_sr6_largex0_q3ex   
     9  +coeff( 36)    *x23            
     1  +coeff( 37)        *x34*x41    
     2  +coeff( 38)    *x22    *x42    
     3  +coeff( 39)*x11*x21        *x51
     4  +coeff( 40)    *x22    *x41*x51
     5  +coeff( 41)    *x21*x32*x41*x51
     6  +coeff( 42)    *x21    *x43*x51
     7  +coeff( 43)    *x23    *x41    
     8  +coeff( 44)*x11*x21    *x42    
      y_sr6_largex0_q3ex   =y_sr6_largex0_q3ex   
     9  +coeff( 45)*x11*x21    *x41*x51
     1  +coeff( 46)    *x21    *x45    
     2  +coeff( 47)    *x22*x33*x41    
     3  +coeff( 48)    *x23*x31*x42    
     4  +coeff( 49)*x11*x21*x33*x41    
     5  +coeff( 50)    *x22*x35*x41    
     6  +coeff( 51)*x12*x22    *x42*x51
     7  +coeff( 52)*x12*x22    *x43*x52
     8  +coeff( 53)*x12*x22    *x41*x54
      y_sr6_largex0_q3ex   =y_sr6_largex0_q3ex   
     9  +coeff( 54)*x12*x22*x32*x43*x52
     1  +coeff( 55)*x11    *x31        
     2  +coeff( 56)    *x21*x32        
     3  +coeff( 57)    *x21    *x41*x51
     4  +coeff( 58)        *x31*x43    
     5  +coeff( 59)    *x21*x31*x42    
     6  +coeff( 60)    *x21    *x43    
     7  +coeff( 61)*x11*x21*x31        
     8  +coeff( 62)*x11*x21    *x41    
      y_sr6_largex0_q3ex   =y_sr6_largex0_q3ex   
     9  +coeff( 63)    *x22*x32        
     1  +coeff( 64)    *x22*x31*x41    
     2  +coeff( 65)    *x22*x31    *x51
     3  +coeff( 66)    *x23*x31        
     4  +coeff( 67)*x11*x21*x31*x41    
     5  +coeff( 68)*x11*x21*x31    *x51
     6  +coeff( 69)    *x22*x31*x41*x51
     7  +coeff( 70)    *x22    *x42*x51
     8  +coeff( 71)    *x21*x34*x41    
      y_sr6_largex0_q3ex   =y_sr6_largex0_q3ex   
     9  +coeff( 72)    *x22    *x44    
     1  +coeff( 73)*x11*x21*x31*x41*x51
     2  +coeff( 74)*x11*x21    *x42*x51
     3  +coeff( 75)    *x23    *x43    
     4  +coeff( 76)    *x22*x34    *x51
     5  +coeff( 77)    *x22*x33*x43    
     6  +coeff( 78)    *x22*x31*x44*x52
     7  +coeff( 79)*x12*x22*x31*x42    
     8  +coeff( 80)*x12*x21        *x54
      y_sr6_largex0_q3ex   =y_sr6_largex0_q3ex   
     9  +coeff( 81)*x12*x22    *x44*x51
     1  +coeff( 82)    *x22*x34*x44*x53
     2  +coeff( 83)*x12*x22*x32*x44*x51
     3  +coeff( 84)*x12*x22*x31*x45*x51
     4  +coeff( 85)*x12*x22*x32*x41*x54
     5  +coeff( 86)        *x32*x41    
     6  +coeff( 87)    *x21*x32*x41    
     7  +coeff( 88)    *x21*x31*x41*x51
     8  +coeff( 89)    *x21    *x42*x51
      y_sr6_largex0_q3ex   =y_sr6_largex0_q3ex   
     9  +coeff( 90)*x12                
c
      return
      end
      function p_sr6_largex0_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2452171E-03/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.45980E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49875E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.52933948E-03, 0.14454763E-01,-0.53398167E-02,-0.68709934E-02,
     +  0.19491377E-02, 0.81867108E-03, 0.17627799E-03, 0.19601185E-01,
     +  0.57804753E-03, 0.27376644E-02,-0.47058058E-02,-0.21830162E-01,
     + -0.18669275E-02, 0.91616763E-03,-0.19660925E-02,-0.14072662E-02,
     +  0.15927080E-01,-0.24997091E-02,-0.33122420E-02, 0.30645591E-02,
     +  0.23259867E-02, 0.33743740E-02, 0.11638580E-01,-0.43538780E-03,
     +  0.14269141E-01, 0.27460738E-02, 0.57202745E-02, 0.12926925E-03,
     +  0.23600829E-02,-0.10131441E-03, 0.29254644E-03,-0.33447403E-03,
     +  0.40725372E-02, 0.25588898E-02, 0.12582383E-02, 0.82062441E-03,
     +  0.11755774E-02, 0.75358595E-03, 0.20068695E-02,-0.50943764E-02,
     + -0.17080681E-02,-0.52871178E-02, 0.18165874E-02, 0.33979076E-02,
     +  0.22400706E-02, 0.78094252E-02, 0.54342719E-02, 0.73788371E-02,
     + -0.19952825E-02,-0.79587707E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      p_sr6_largex0_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      p_sr6_largex0_q3ex   =p_sr6_largex0_q3ex   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22    *x41    
      p_sr6_largex0_q3ex   =p_sr6_largex0_q3ex   
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)            *x43    
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)            *x41*x52
     5  +coeff( 23)    *x22*x31*x41    
     6  +coeff( 24)        *x33*x41    
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)    *x22    *x41*x51
      p_sr6_largex0_q3ex   =p_sr6_largex0_q3ex   
     9  +coeff( 27)*x11*x21    *x42    
     1  +coeff( 28)    *x22*x33        
     2  +coeff( 29)    *x22*x31    *x52
     3  +coeff( 30)    *x22*x33    *x52
     4  +coeff( 31)    *x21            
     5  +coeff( 32)*x11            *x51
     6  +coeff( 33)    *x22*x31        
     7  +coeff( 34)        *x31*x42    
     8  +coeff( 35)    *x21    *x41*x51
      p_sr6_largex0_q3ex   =p_sr6_largex0_q3ex   
     9  +coeff( 36)            *x42*x51
     1  +coeff( 37)*x11*x21    *x41    
     2  +coeff( 38)*x11*x21        *x51
     3  +coeff( 39)    *x22*x32        
     4  +coeff( 40)    *x23    *x41    
     5  +coeff( 41)        *x32*x42    
     6  +coeff( 42)    *x21    *x43    
     7  +coeff( 43)    *x22*x31    *x51
     8  +coeff( 44)*x11*x21*x31*x41    
      p_sr6_largex0_q3ex   =p_sr6_largex0_q3ex   
     9  +coeff( 45)*x11*x21    *x41*x51
     1  +coeff( 46)    *x22    *x43    
     2  +coeff( 47)    *x22*x31*x41*x51
     3  +coeff( 48)    *x22    *x42*x51
     4  +coeff( 49)*x11*x22    *x42    
     5  +coeff( 50)    *x23*x31*x42    
c
      return
      end
      function l_sr6_largex0_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.2135496E+02/
      data xmin/
     1 -0.15981E-01,-0.52463E-01,-0.19936E-01,-0.29294E-01,-0.45980E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15995E-01, 0.53602E-01, 0.19901E-01, 0.29979E-01, 0.49875E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12805233E-01, 0.31410015E+00, 0.31332985E-01,-0.33618424E-01,
     +  0.33018421E-01,-0.42657219E-03,-0.68102479E-01,-0.22312826E-02,
     +  0.12973417E-01, 0.77458322E-02,-0.59065539E-02,-0.44227626E-01,
     +  0.53923246E-02, 0.71282363E-02, 0.29851040E-02,-0.11419739E-02,
     + -0.90863453E-02, 0.14020738E-02,-0.75824715E-01,-0.62383182E-01,
     + -0.15872981E-02, 0.21624921E-02,-0.20705413E-02,-0.24862096E-02,
     +  0.61237724E-02,-0.32726988E-01,-0.12916302E-01, 0.13709455E-02,
     + -0.45065368E-02,-0.43007624E-02, 0.32647062E-01,-0.66497959E-02,
     +  0.17834137E-02, 0.37761605E-02,-0.17022347E-02, 0.16082476E-02,
     + -0.51335325E-02,-0.72123948E-02,-0.84151700E-02,-0.25581494E-02,
     +  0.15800816E-02,-0.63884594E-02,-0.31406511E-03, 0.29917839E-02,
     + -0.13624256E-01,-0.81417281E-02,-0.95233740E-02, 0.39782785E-03,
     + -0.48787892E-02,-0.33428568E-01,-0.15631475E-01, 0.15846955E-01,
     +  0.96032774E-04,-0.21985765E-02,-0.22302702E-03,-0.43991003E-02,
     + -0.14536667E-02,-0.36134997E-02,-0.15880659E-02,-0.54182932E-02,
     + -0.11622461E-01, 0.45104648E-03,-0.12742706E-02,-0.17974448E-02,
     + -0.68406120E-03, 0.25251543E-01,-0.92848046E-02,-0.11882992E-01,
     + -0.16013773E-01,-0.50664740E-02,-0.89130970E-02, 0.20188313E-01,
     +  0.11205490E-01,-0.29224845E-01,-0.22155991E-01, 0.14971737E-01,
     + -0.76634497E-02,-0.24472198E-01, 0.37796460E-01,-0.15036559E-01,
     + -0.50988916E-01,-0.17406380E-01,-0.35238288E-01, 0.31414232E-03,
     + -0.58254047E-03,-0.16640924E-02,-0.39583340E-03,-0.21036621E-02,
     + -0.21501549E-02,-0.83008222E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sr6_largex0_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x21*x31        
      l_sr6_largex0_q3ex   =l_sr6_largex0_q3ex   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)            *x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)    *x21*x32        
      l_sr6_largex0_q3ex   =l_sr6_largex0_q3ex   
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)    *x21*x31*x42    
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)    *x23*x32*x41    
     4  +coeff( 22)        *x31        
     5  +coeff( 23)*x11        *x41    
     6  +coeff( 24)    *x22    *x41    
     7  +coeff( 25)    *x23    *x41    
     8  +coeff( 26)    *x21*x32*x41    
      l_sr6_largex0_q3ex   =l_sr6_largex0_q3ex   
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)*x11*x22*x31        
     2  +coeff( 29)            *x42*x53
     3  +coeff( 30)*x11*x22    *x42    
     4  +coeff( 31)    *x21*x32*x43    
     5  +coeff( 32)*x11*x22*x33*x41    
     6  +coeff( 33)        *x31*x41    
     7  +coeff( 34)            *x41*x51
     8  +coeff( 35)*x11    *x31        
      l_sr6_largex0_q3ex   =l_sr6_largex0_q3ex   
     9  +coeff( 36)    *x22        *x51
     1  +coeff( 37)*x11    *x31*x41    
     2  +coeff( 38)*x11        *x42    
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)    *x22    *x41*x51
     5  +coeff( 41)    *x21    *x42*x51
     6  +coeff( 42)            *x43*x51
     7  +coeff( 43)    *x23*x31*x41    
     8  +coeff( 44)    *x21*x32*x42    
      l_sr6_largex0_q3ex   =l_sr6_largex0_q3ex   
     9  +coeff( 45)    *x22    *x43    
     1  +coeff( 46)    *x21*x31*x43    
     2  +coeff( 47)    *x23    *x41*x51
     3  +coeff( 48)        *x33*x41*x51
     4  +coeff( 49)*x11*x22*x32        
     5  +coeff( 50)    *x23    *x43    
     6  +coeff( 51)*x11*x22*x31*x42*x52
     7  +coeff( 52)*x11*x23    *x43*x53
     8  +coeff( 53)        *x31    *x51
      l_sr6_largex0_q3ex   =l_sr6_largex0_q3ex   
     9  +coeff( 54)        *x31*x42    
     1  +coeff( 55)*x11        *x41*x51
     2  +coeff( 56)    *x21*x33        
     3  +coeff( 57)    *x21*x31*x41*x51
     4  +coeff( 58)        *x31*x42*x51
     5  +coeff( 59)            *x41*x53
     6  +coeff( 60)*x11        *x43    
     7  +coeff( 61)    *x22*x31*x42    
     8  +coeff( 62)        *x33*x42    
      l_sr6_largex0_q3ex   =l_sr6_largex0_q3ex   
     9  +coeff( 63)        *x32*x43    
     1  +coeff( 64)        *x31*x41*x53
     2  +coeff( 65)*x12*x21*x31*x41    
     3  +coeff( 66)    *x21*x33*x42    
     4  +coeff( 67)*x12*x21    *x41*x51
     5  +coeff( 68)    *x23    *x42*x51
     6  +coeff( 69)*x11*x22*x31*x42    
     7  +coeff( 70)*x11    *x33*x42    
     8  +coeff( 71)*x12*x21    *x42*x51
      l_sr6_largex0_q3ex   =l_sr6_largex0_q3ex   
     9  +coeff( 72)    *x23    *x42*x52
     1  +coeff( 73)*x12*x21*x33*x41    
     2  +coeff( 74)    *x23*x33*x42    
     3  +coeff( 75)*x12*x21*x31*x43    
     4  +coeff( 76)*x12*x23    *x41*x51
     5  +coeff( 77)    *x21*x31*x43*x53
     6  +coeff( 78)*x11*x22*x32*x43    
     7  +coeff( 79)*x11*x23    *x43*x51
     8  +coeff( 80)*x11*x23    *x41*x53
      l_sr6_largex0_q3ex   =l_sr6_largex0_q3ex   
     9  +coeff( 81)    *x23*x32*x42*x52
     1  +coeff( 82)*x11*x22    *x43*x53
     2  +coeff( 83)*x11*x23*x32*x43*x51
     3  +coeff( 84)*x12                
     4  +coeff( 85)    *x23            
     5  +coeff( 86)        *x32*x41    
     6  +coeff( 87)    *x21*x31    *x51
     7  +coeff( 88)        *x31*x41*x51
     8  +coeff( 89)            *x42*x51
      l_sr6_largex0_q3ex   =l_sr6_largex0_q3ex   
     9  +coeff( 90)                *x53
c
      return
      end
      function x_sr6_largex0_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.7177365E-02/
      data xmin/
     1 -0.15997E-01,-0.54731E-01,-0.19992E-01,-0.29286E-01,-0.46196E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15975E-01, 0.53055E-01, 0.19994E-01, 0.29237E-01, 0.49732E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15519223E-01, 0.33404874E-02, 0.33097487E-03, 0.62524176E+00,
     + -0.36002669E-01, 0.52811999E-01,-0.48912626E-01, 0.21338204E-01,
     + -0.10852262E-01,-0.36974780E-01,-0.64234814E-03,-0.72118390E-03,
     + -0.19316271E-01,-0.61379541E-02, 0.49725333E-02,-0.29898926E-02,
     + -0.28850723E-02, 0.44212127E-03,-0.13249891E-02, 0.44148834E-02,
     + -0.46933284E-02,-0.20368849E-02,-0.86380746E-02,-0.53862837E-03,
     + -0.68353402E-03, 0.36315324E-02, 0.38918303E-02,-0.65439183E-03,
     +  0.37310645E-02,-0.25474355E-01,-0.19114476E-01,-0.85550966E-02,
     + -0.14070253E-02,-0.20619179E-02, 0.25068101E-03, 0.23864952E-02,
     + -0.40023707E-03, 0.16395423E-02,-0.28537650E-03,-0.15192491E-02,
     + -0.73851604E-03, 0.28537202E-03,-0.27208799E-02,-0.26948277E-02,
     + -0.99823875E-02, 0.74693601E-03,-0.56587211E-02,-0.72813168E-03,
     +  0.98474463E-03, 0.13561606E-01, 0.82436111E-02,-0.18490726E-01,
     + -0.43845523E-03,-0.42050099E-03,-0.12702325E-03, 0.26002212E-02,
     +  0.43813014E-03,-0.22959381E-02,-0.15436839E-02,-0.15170017E-02,
     +  0.10560072E-02, 0.10926549E-02, 0.12056521E-02, 0.87310653E-03,
     +  0.20057836E-02,-0.30437862E-02,-0.11970084E-02,-0.44663236E-02,
     +  0.17538922E-02,-0.32356102E-02,-0.12418665E-02,-0.96706010E-03,
     + -0.73415474E-02, 0.76910283E-03,-0.22727433E-02,-0.83735166E-02,
     + -0.49485117E-02, 0.10728075E-02,-0.27823942E-02,-0.20134666E-02,
     + -0.80144237E-03,-0.36224369E-02,-0.10565006E-01,-0.27413974E-02,
     + -0.25927825E-02, 0.11654785E-02, 0.24366709E-02, 0.36235051E-02,
     +  0.44583292E-02, 0.11655616E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x_sr6_largex0_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x22            
      x_sr6_largex0_fp     =x_sr6_largex0_fp     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)                *x53
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)            *x41*x51
      x_sr6_largex0_fp     =x_sr6_largex0_fp     
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x21*x32        
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)    *x24            
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)    *x21*x31    *x51
     8  +coeff( 26)    *x23    *x41    
      x_sr6_largex0_fp     =x_sr6_largex0_fp     
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)    *x22*x31    *x51
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)    *x21*x31*x42    
     4  +coeff( 31)    *x21    *x43    
     5  +coeff( 32)    *x21    *x42*x51
     6  +coeff( 33)    *x23*x32*x41    
     7  +coeff( 34)    *x21*x31*x41*x53
     8  +coeff( 35)        *x32        
      x_sr6_largex0_fp     =x_sr6_largex0_fp     
     9  +coeff( 36)    *x23            
     1  +coeff( 37)*x11    *x31        
     2  +coeff( 38)    *x22    *x41    
     3  +coeff( 39)*x12                
     4  +coeff( 40)*x11*x22            
     5  +coeff( 41)        *x31*x41*x51
     6  +coeff( 42)        *x31    *x52
     7  +coeff( 43)*x11    *x31*x41    
     8  +coeff( 44)*x11        *x42    
      x_sr6_largex0_fp     =x_sr6_largex0_fp     
     9  +coeff( 45)    *x21*x32*x41    
     1  +coeff( 46)    *x24*x31        
     2  +coeff( 47)    *x24        *x51
     3  +coeff( 48)*x12*x22    *x41    
     4  +coeff( 49)    *x24        *x53
     5  +coeff( 50)    *x24    *x43*x52
     6  +coeff( 51)    *x24    *x42*x53
     7  +coeff( 52)*x11*x24    *x43*x51
     8  +coeff( 53)        *x31    *x51
      x_sr6_largex0_fp     =x_sr6_largex0_fp     
     9  +coeff( 54)*x11    *x32        
     1  +coeff( 55)*x11    *x31    *x51
     2  +coeff( 56)    *x22    *x41*x51
     3  +coeff( 57)*x11            *x52
     4  +coeff( 58)    *x22        *x52
     5  +coeff( 59)    *x21*x33        
     6  +coeff( 60)    *x21    *x41*x52
     7  +coeff( 61)*x11*x22    *x41    
     8  +coeff( 62)*x11*x22        *x51
      x_sr6_largex0_fp     =x_sr6_largex0_fp     
     9  +coeff( 63)        *x31*x41*x52
     1  +coeff( 64)            *x42*x52
     2  +coeff( 65)    *x23    *x41*x51
     3  +coeff( 66)    *x23        *x52
     4  +coeff( 67)*x12*x22            
     5  +coeff( 68)*x11    *x31*x42    
     6  +coeff( 69)    *x22*x31*x42    
     7  +coeff( 70)*x11        *x43    
     8  +coeff( 71)*x11    *x31*x41*x51
      x_sr6_largex0_fp     =x_sr6_largex0_fp     
     9  +coeff( 72)*x11        *x42*x51
     1  +coeff( 73)    *x21*x31*x43    
     2  +coeff( 74)    *x21*x33    *x51
     3  +coeff( 75)    *x21*x32*x41*x51
     4  +coeff( 76)    *x21*x31*x42*x51
     5  +coeff( 77)    *x21    *x43*x51
     6  +coeff( 78)*x11*x22    *x41*x51
     7  +coeff( 79)    *x24    *x41*x51
     8  +coeff( 80)    *x23*x32    *x51
      x_sr6_largex0_fp     =x_sr6_largex0_fp     
     9  +coeff( 81)    *x23        *x53
     1  +coeff( 82)    *x21*x33*x41*x51
     2  +coeff( 83)    *x21*x31*x43*x51
     3  +coeff( 84)*x11*x22*x32*x41    
     4  +coeff( 85)*x11*x22*x31*x42    
     5  +coeff( 86)*x12*x22*x32        
     6  +coeff( 87)*x11    *x33*x42    
     7  +coeff( 88)    *x21*x33*x43    
     8  +coeff( 89)*x11*x22*x32*x41*x51
      x_sr6_largex0_fp     =x_sr6_largex0_fp     
     9  +coeff( 90)    *x21*x33*x43*x51
c
      return
      end
      function t_sr6_largex0_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.6056797E-04/
      data xmin/
     1 -0.15997E-01,-0.54731E-01,-0.19992E-01,-0.29286E-01,-0.46196E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15975E-01, 0.53055E-01, 0.19994E-01, 0.29237E-01, 0.49732E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.40919916E-02,-0.23652343E-01, 0.93496201E-04, 0.10737766E+00,
     + -0.36734908E-02, 0.59218369E-02,-0.97843157E-02, 0.83670515E-03,
     + -0.18206316E-02, 0.59338205E-03,-0.88851183E-03, 0.58478803E-04,
     + -0.56859036E-03,-0.50747179E-03,-0.11562496E-03,-0.12214467E-02,
     + -0.70003938E-03, 0.10424824E-02, 0.17954644E-02, 0.47177651E-04,
     +  0.32382392E-04,-0.13393170E-03,-0.17527382E-03, 0.25819227E-03,
     + -0.31798679E-03, 0.23301609E-03, 0.93987543E-03,-0.62730373E-03,
     + -0.10367906E-03, 0.21568784E-02, 0.10836763E-03, 0.10304329E-03,
     + -0.73941068E-04,-0.16263132E-03,-0.75310716E-04, 0.67975605E-03,
     +  0.11172013E-02, 0.59233775E-03,-0.33578969E-03,-0.24351705E-03,
     +  0.29368079E-03, 0.35609966E-03, 0.79759462E-04,-0.32778911E-03,
     +  0.87848259E-03,-0.45524101E-03, 0.11953469E-02,-0.25815975E-02,
     +  0.29138285E-02, 0.12009541E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sr6_largex0_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x22            
      t_sr6_largex0_fp     =t_sr6_largex0_fp     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x21        *x52
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x21    *x41*x51
      t_sr6_largex0_fp     =t_sr6_largex0_fp     
     9  +coeff( 18)                *x53
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)        *x31        
     3  +coeff( 21)    *x21*x31        
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)    *x23            
     7  +coeff( 25)        *x31*x41*x51
     8  +coeff( 26)*x11        *x42    
      t_sr6_largex0_fp     =t_sr6_largex0_fp     
     9  +coeff( 27)    *x21    *x43    
     1  +coeff( 28)    *x21*x32*x42    
     2  +coeff( 29)    *x21*x33*x42    
     3  +coeff( 30)    *x22    *x43*x52
     4  +coeff( 31)        *x32        
     5  +coeff( 32)*x11        *x41    
     6  +coeff( 33)*x12                
     7  +coeff( 34)        *x32    *x51
     8  +coeff( 35)    *x21*x33        
      t_sr6_largex0_fp     =t_sr6_largex0_fp     
     9  +coeff( 36)    *x22*x31*x41    
     1  +coeff( 37)    *x21*x31*x42    
     2  +coeff( 38)    *x22    *x41*x51
     3  +coeff( 39)    *x21    *x42*x51
     4  +coeff( 40)    *x22        *x52
     5  +coeff( 41)        *x31*x41*x52
     6  +coeff( 42)*x11*x23            
     7  +coeff( 43)*x11    *x33        
     8  +coeff( 44)    *x22*x31*x41*x51
      t_sr6_largex0_fp     =t_sr6_largex0_fp     
     9  +coeff( 45)    *x23    *x43    
     1  +coeff( 46)    *x23*x31*x42*x51
     2  +coeff( 47)    *x22*x31*x42*x52
     3  +coeff( 48)    *x22    *x43*x53
     4  +coeff( 49)    *x22*x32*x43*x53
     5  +coeff( 50)    *x22    *x41    
c
      return
      end
      function y_sr6_largex0_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.2754166E-02/
      data xmin/
     1 -0.15997E-01,-0.54731E-01,-0.19992E-01,-0.29286E-01,-0.46196E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15975E-01, 0.53055E-01, 0.19994E-01, 0.29237E-01, 0.49732E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24165474E-02,-0.43055127E-03,-0.40170141E-01,-0.69118035E-02,
     + -0.38493101E-02,-0.72026621E-02,-0.92359185E-02, 0.10082520E-02,
     + -0.13606400E-02, 0.39999317E-02, 0.15964782E-02, 0.92148368E-03,
     + -0.21489041E-02,-0.27118607E-02, 0.19793832E-02, 0.79024760E-02,
     + -0.67806475E-04,-0.19496267E-03, 0.78655622E-03,-0.38282403E-02,
     + -0.42102774E-03,-0.29891231E-02, 0.43372260E-02, 0.92108911E-02,
     +  0.70109073E-03,-0.90181560E-03, 0.47154957E-02, 0.91958204E-02,
     + -0.43428806E-02,-0.70557115E-03, 0.80181379E-02, 0.80321766E-02,
     +  0.13510721E-02, 0.14588614E-01, 0.17038593E-01,-0.73211256E-03,
     +  0.35202492E-02,-0.52600065E-02, 0.19869986E-02, 0.44153947E-02,
     +  0.63909516E-02,-0.22242347E-01,-0.37758183E-02, 0.55819657E-01,
     +  0.54583766E-01, 0.66187404E-01, 0.85904228E-03,-0.10087312E-02,
     + -0.42553917E-02,-0.54100171E-01, 0.77664125E-03,-0.32292953E-03,
     + -0.67815511E-03, 0.25993199E-02, 0.13633752E-01, 0.73316967E-03,
     +  0.32171102E-02,-0.16507454E-02, 0.81294058E-02,-0.20673737E-04,
     +  0.10073680E-01, 0.69398931E-02, 0.27181976E-02, 0.66548530E-02,
     +  0.36874851E-02,-0.27593502E-02, 0.25226854E-03,-0.19136235E-02,
     +  0.27945198E-02,-0.56945166E-03,-0.20859500E-02,-0.32903028E-02,
     + -0.38763479E-03,-0.10300075E-02, 0.21849504E-03,-0.32882474E-02,
     + -0.12769568E-02, 0.34442256E-03,-0.20741718E-02, 0.17434987E-02,
     +  0.28642118E-02, 0.80769714E-02, 0.24508804E-02, 0.12247362E-02,
     + -0.18046538E-02,-0.20855076E-02,-0.43660165E-02, 0.21941951E-02,
     +  0.14561652E-02, 0.32227777E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sr6_largex0_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)            *x42    
     6  +coeff(  6)        *x31    *x51
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)                *x52
      y_sr6_largex0_fp     =y_sr6_largex0_fp     
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x22            
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)    *x21*x32        
      y_sr6_largex0_fp     =y_sr6_largex0_fp     
     9  +coeff( 18)*x11            *x51
     1  +coeff( 19)    *x21*x31    *x51
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)        *x34        
     4  +coeff( 22)        *x33*x41    
     5  +coeff( 23)        *x31*x43    
     6  +coeff( 24)            *x44    
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)    *x21        *x52
      y_sr6_largex0_fp     =y_sr6_largex0_fp     
     9  +coeff( 27)    *x22*x31        
     1  +coeff( 28)    *x22    *x41    
     2  +coeff( 29)            *x41*x53
     3  +coeff( 30)    *x23            
     4  +coeff( 31)        *x32*x43    
     5  +coeff( 32)            *x45    
     6  +coeff( 33)*x11*x21    *x41    
     7  +coeff( 34)    *x22*x31*x41    
     8  +coeff( 35)    *x22    *x42    
      y_sr6_largex0_fp     =y_sr6_largex0_fp     
     9  +coeff( 36)*x11    *x32*x41    
     1  +coeff( 37)        *x33    *x52
     2  +coeff( 38)    *x22    *x41*x51
     3  +coeff( 39)    *x21    *x44    
     4  +coeff( 40)        *x31    *x54
     5  +coeff( 41)*x11*x21    *x42    
     6  +coeff( 42)    *x22    *x43    
     7  +coeff( 43)*x11*x21*x33*x41    
     8  +coeff( 44)    *x22*x32*x43    
      y_sr6_largex0_fp     =y_sr6_largex0_fp     
     9  +coeff( 45)    *x22*x31*x44    
     1  +coeff( 46)    *x22    *x45    
     2  +coeff( 47)*x11*x23*x31        
     3  +coeff( 48)*x11*x23*x31*x41    
     4  +coeff( 49)        *x35    *x54
     5  +coeff( 50)    *x22*x32*x45    
     6  +coeff( 51)            *x42*x51
     7  +coeff( 52)    *x21*x31*x41    
     8  +coeff( 53)    *x22        *x51
      y_sr6_largex0_fp     =y_sr6_largex0_fp     
     9  +coeff( 54)        *x33*x42    
     1  +coeff( 55)        *x31*x44    
     2  +coeff( 56)        *x33*x41*x51
     3  +coeff( 57)    *x22*x32        
     4  +coeff( 58)    *x22*x31    *x51
     5  +coeff( 59)*x11*x21*x31*x41    
     6  +coeff( 60)    *x21    *x45    
     7  +coeff( 61)*x11*x21*x31*x42    
     8  +coeff( 62)*x11*x21    *x43    
      y_sr6_largex0_fp     =y_sr6_largex0_fp     
     9  +coeff( 63)*x11*x21    *x43*x51
     1  +coeff( 64)    *x23*x32*x41*x51
     2  +coeff( 65)    *x22    *x41*x54
     3  +coeff( 66)*x11*x21*x35*x42    
     4  +coeff( 67)*x11    *x31        
     5  +coeff( 68)    *x21    *x42    
     6  +coeff( 69)        *x32*x42    
     7  +coeff( 70)            *x42*x52
     8  +coeff( 71)    *x21*x31*x42    
      y_sr6_largex0_fp     =y_sr6_largex0_fp     
     9  +coeff( 72)    *x21    *x43    
     1  +coeff( 73)*x11    *x31    *x51
     2  +coeff( 74)        *x31    *x53
     3  +coeff( 75)*x12                
     4  +coeff( 76)        *x34*x41    
     5  +coeff( 77)    *x23*x31        
     6  +coeff( 78)*x12        *x41    
     7  +coeff( 79)    *x23    *x41    
     8  +coeff( 80)        *x35*x41    
      y_sr6_largex0_fp     =y_sr6_largex0_fp     
     9  +coeff( 81)        *x32*x44    
     1  +coeff( 82)        *x31*x45    
     2  +coeff( 83)            *x41*x54
     3  +coeff( 84)*x11*x21*x32        
     4  +coeff( 85)    *x21    *x42*x52
     5  +coeff( 86)    *x23*x31*x41    
     6  +coeff( 87)        *x35*x42    
     7  +coeff( 88)*x11*x21*x32*x41    
     8  +coeff( 89)    *x23*x31    *x51
      y_sr6_largex0_fp     =y_sr6_largex0_fp     
     9  +coeff( 90)*x11            *x54
c
      return
      end
      function p_sr6_largex0_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.5383719E-03/
      data xmin/
     1 -0.15997E-01,-0.54731E-01,-0.19992E-01,-0.29286E-01,-0.46196E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15975E-01, 0.53055E-01, 0.19994E-01, 0.29237E-01, 0.49732E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12668234E-02,-0.16924087E-03, 0.14190837E-01,-0.61549768E-02,
     + -0.65711075E-02,-0.90102418E-04, 0.17806503E-02, 0.42915944E-03,
     +  0.16603874E-01, 0.44243116E-03, 0.30421363E-02,-0.42192307E-02,
     + -0.21524107E-01,-0.20564022E-02, 0.15288853E-02,-0.31398027E-03,
     + -0.18050198E-02,-0.30112927E-03,-0.13857365E-02, 0.17517492E-01,
     +  0.13334283E-02,-0.31634214E-03, 0.45767762E-02, 0.49517970E-02,
     +  0.15415601E-02, 0.11533208E-02, 0.22465170E-02, 0.19485989E-02,
     +  0.87120623E-03,-0.20504877E-03, 0.11454385E-01,-0.26250808E-03,
     +  0.14271175E-01, 0.23125042E-02, 0.32884060E-03,-0.74634878E-02,
     + -0.94411159E-02, 0.91537256E-02, 0.19332959E-02, 0.14928274E-01,
     + -0.25494689E-01,-0.14698378E-01, 0.16233881E-03, 0.66791184E-03,
     +  0.34385722E-03, 0.40964503E-02, 0.13594663E-02, 0.10623023E-02,
     +  0.11063013E-02, 0.18490975E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sr6_largex0_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      p_sr6_largex0_fp     =p_sr6_largex0_fp     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11    *x31        
     8  +coeff( 17)*x11        *x41    
      p_sr6_largex0_fp     =p_sr6_largex0_fp     
     9  +coeff( 18)*x11            *x51
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)        *x32*x41    
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)        *x31*x42    
     6  +coeff( 24)            *x43    
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)            *x42*x51
      p_sr6_largex0_fp     =p_sr6_largex0_fp     
     9  +coeff( 27)            *x41*x52
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)*x11*x21        *x51
     3  +coeff( 30)    *x23*x31        
     4  +coeff( 31)    *x22*x31*x41    
     5  +coeff( 32)        *x33*x41    
     6  +coeff( 33)    *x22    *x42    
     7  +coeff( 34)    *x22    *x41*x51
     8  +coeff( 35)    *x22*x33        
      p_sr6_largex0_fp     =p_sr6_largex0_fp     
     9  +coeff( 36)    *x23*x31*x41    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)    *x22*x33*x41*x51
     4  +coeff( 40)    *x22    *x43*x52
     5  +coeff( 41)    *x23    *x43*x52
     6  +coeff( 42)    *x23    *x42*x53
     7  +coeff( 43)        *x32        
     8  +coeff( 44)        *x31*x41    
      p_sr6_largex0_fp     =p_sr6_largex0_fp     
     9  +coeff( 45)*x12                
     1  +coeff( 46)    *x22*x31        
     2  +coeff( 47)        *x31*x41*x51
     3  +coeff( 48)        *x31    *x52
     4  +coeff( 49)*x11*x21*x31        
     5  +coeff( 50)    *x22*x32        
c
      return
      end
      function l_sr6_largex0_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2209553E-01/
      data xmin/
     1 -0.15997E-01,-0.54731E-01,-0.19992E-01,-0.29286E-01,-0.46196E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15975E-01, 0.53055E-01, 0.19994E-01, 0.29237E-01, 0.49732E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25713673E-01,-0.31662652E+00,-0.31892922E-01, 0.32449625E-01,
     + -0.34903929E-01,-0.83691982E-03,-0.24533965E-01, 0.66045500E-01,
     + -0.49765245E-02, 0.11895704E-02, 0.57513742E-02, 0.43160066E-01,
     + -0.49247537E-02,-0.27705540E-02,-0.84384112E-02,-0.55399127E-02,
     + -0.26329624E-03, 0.91159930E-02,-0.33959761E-03, 0.50034104E-02,
     +  0.78396857E-01, 0.62125847E-01, 0.10336152E-01,-0.20941044E-02,
     + -0.21093492E-02, 0.80038805E-03,-0.17386887E-02, 0.14328888E-02,
     + -0.42910851E-02, 0.44886721E-02, 0.12735718E-01,-0.71870163E-02,
     +  0.30998213E-01, 0.11916257E-01,-0.39609021E-03,-0.24385394E-02,
     + -0.15164145E-01,-0.64858136E-03, 0.10676009E-02,-0.18800950E-02,
     +  0.30119070E-02, 0.77453507E-02, 0.49635484E-02, 0.73652337E-02,
     +  0.44291103E-02, 0.18140623E-01, 0.14159530E-01, 0.11519653E-01,
     +  0.99576963E-03,-0.65344088E-02,-0.17099483E-02,-0.72853472E-02,
     + -0.90316701E-03, 0.22286228E-02, 0.15383479E-02, 0.16414820E-02,
     +  0.78568729E-02,-0.24614201E-02,-0.23179548E-02, 0.38701352E-02,
     +  0.14049131E-02, 0.58401008E-02,-0.91244588E-02, 0.75094872E-02,
     + -0.48594293E-03, 0.40710294E-02,-0.30055549E-02, 0.78133540E-03,
     + -0.79131899E-02, 0.19009768E-02, 0.68704705E-02,-0.14839772E-01,
     +  0.19337017E-01,-0.25480427E-02, 0.12317301E-01,-0.41127363E-02,
     + -0.36983998E-02, 0.65470524E-02, 0.37732781E-02, 0.57332111E-02,
     + -0.58036246E-02,-0.16786017E-02, 0.86751534E-02, 0.77615899E-03,
     + -0.27284466E-01,-0.51435594E-01,-0.84955823E-02, 0.12774769E-01,
     +  0.80891550E-02, 0.53324204E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (abs(xmin(i)-xmax(i))<1.0E-08) go to 5
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sr6_largex0_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x42    
      l_sr6_largex0_fp     =l_sr6_largex0_fp     
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)            *x41    
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)            *x42    
     7  +coeff( 16)    *x21        *x51
     8  +coeff( 17)*x11        *x41    
      l_sr6_largex0_fp     =l_sr6_largex0_fp     
     9  +coeff( 18)    *x21*x32        
     1  +coeff( 19)        *x33        
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)    *x21*x31*x42    
     4  +coeff( 22)    *x21    *x43    
     5  +coeff( 23)    *x23*x32*x41    
     6  +coeff( 24)        *x31        
     7  +coeff( 25)            *x41*x51
     8  +coeff( 26)*x11    *x31        
      l_sr6_largex0_fp     =l_sr6_largex0_fp     
     9  +coeff( 27)*x11            *x51
     1  +coeff( 28)    *x22    *x41    
     2  +coeff( 29)    *x21        *x52
     3  +coeff( 30)                *x53
     4  +coeff( 31)*x11        *x42    
     5  +coeff( 32)    *x23    *x41    
     6  +coeff( 33)    *x21*x32*x41    
     7  +coeff( 34)    *x22    *x42    
     8  +coeff( 35)*x11    *x33*x41    
      l_sr6_largex0_fp     =l_sr6_largex0_fp     
     9  +coeff( 36)    *x22*x33*x41    
     1  +coeff( 37)    *x21*x32*x43    
     2  +coeff( 38)        *x32        
     3  +coeff( 39)    *x22*x31        
     4  +coeff( 40)    *x22        *x51
     5  +coeff( 41)*x11*x22            
     6  +coeff( 42)*x11    *x31*x41    
     7  +coeff( 43)    *x21*x33        
     8  +coeff( 44)    *x21    *x42*x51
      l_sr6_largex0_fp     =l_sr6_largex0_fp     
     9  +coeff( 45)*x11*x21    *x42    
     1  +coeff( 46)*x11    *x31*x42    
     2  +coeff( 47)*x11        *x43    
     3  +coeff( 48)    *x21*x31*x43    
     4  +coeff( 49)    *x23    *x41*x52
     5  +coeff( 50)*x11*x21*x32*x42    
     6  +coeff( 51)*x11    *x33*x42    
     7  +coeff( 52)    *x21    *x43*x53
     8  +coeff( 53)        *x32    *x51
      l_sr6_largex0_fp     =l_sr6_largex0_fp     
     9  +coeff( 54)            *x42*x51
     1  +coeff( 55)*x11    *x32        
     2  +coeff( 56)*x11*x21    *x41    
     3  +coeff( 57)    *x22*x31*x41    
     4  +coeff( 58)        *x32*x42    
     5  +coeff( 59)        *x31*x43    
     6  +coeff( 60)    *x21*x31*x41*x51
     7  +coeff( 61)            *x41*x53
     8  +coeff( 62)*x11    *x32*x41    
      l_sr6_largex0_fp     =l_sr6_largex0_fp     
     9  +coeff( 63)    *x23    *x42    
     1  +coeff( 64)    *x22    *x43    
     2  +coeff( 65)*x12            *x52
     3  +coeff( 66)    *x23        *x52
     4  +coeff( 67)            *x42*x53
     5  +coeff( 68)*x11*x21*x33        
     6  +coeff( 69)*x11*x22    *x42    
     7  +coeff( 70)*x11    *x32*x42    
     8  +coeff( 71)*x11    *x31*x43    
      l_sr6_largex0_fp     =l_sr6_largex0_fp     
     9  +coeff( 72)    *x21*x33*x42    
     1  +coeff( 73)    *x23    *x43    
     2  +coeff( 74)    *x22*x32*x41*x51
     3  +coeff( 75)    *x22    *x43*x51
     4  +coeff( 76)        *x31*x42*x53
     5  +coeff( 77)            *x43*x53
     6  +coeff( 78)    *x22*x33*x42    
     7  +coeff( 79)*x12*x22*x31    *x51
     8  +coeff( 80)*x12*x22    *x41*x51
      l_sr6_largex0_fp     =l_sr6_largex0_fp     
     9  +coeff( 81)    *x22*x32*x42*x51
     1  +coeff( 82)        *x31*x43*x53
     2  +coeff( 83)*x11*x22*x32*x42    
     3  +coeff( 84)    *x23*x33*x42    
     4  +coeff( 85)    *x23*x32*x43    
     5  +coeff( 86)    *x22    *x43*x53
     6  +coeff( 87)*x12*x22    *x42*x52
     7  +coeff( 88)    *x22*x32*x42*x53
     8  +coeff( 89)    *x22*x33*x42*x53
      l_sr6_largex0_fp     =l_sr6_largex0_fp     
     9  +coeff( 90)    *x22*x32*x43*x53
c
      return
      end
