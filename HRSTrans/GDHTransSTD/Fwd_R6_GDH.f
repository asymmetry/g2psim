C Forward transfer functions for right hrs with septum based on Rtest_dir.dat
c                     -JJL 8/06/03


c typical call: answer = function(x,5)
c INPUTS: x = 5 or more element array 
c              x(1)=x0  (meters)
c              x(2)=theta0 (really tan(theta0))
c              x(3)=y0   (meters)
c              x(4)=phi0 (really tan(phi0))
c              x(5)=delta (fractional value NOT percent)
c         M=5
c
c OUTPUT: units are the same as inputs
c 
c NOMENCLATURE: function name = prefix + _sr_ +suffix
c           prefixes:     x means xfinal
c                         t means thetafinal
c                         y means yfinal
c                         p means phifinal
c                         l means path length
c     
c           suffixes:     fp means target to focus
c                         q1ex means target to Q1 exit
c                         dent means target to Dipole entrance
c                         dext means target to dipole exit
c                         q3en means target to Q3 entrance
c                         q3ex means target to Q3 exit
c                         ep3  means septum entrance
c                         ep4  means one quarter the way into the septum
c                         ep5  means halfway through the septum
c                         ep6  means three quarters the way through the septum
c                         ep7  means septum exit
c
c          _sr_ is for right hrs with septum
c
c APERTURES:
c    Coordinate systems are different in the spectrometer with septum model compared to the spectrometer
c     without septum. So some apertures even in the body of the spectrometer appear to be different. Numbers
c    given here supercede any numbers given in regard to transfer functions for the spectrometers alone.
c
c     ep3: -0.1486233 < x < -0.08869672
c          -0.110 < y < 0.110
c     ep4: -0.1792231 < x < -0.1089169
c          -0.110 < y < 0.110
c     ep5: -0.2209211 < x < -0.1353789
c          -0.110 < y < 0.110
c     ep6: -0.2763536 < x < -0.1697464
c          -0.110 < y < 0.110
c     ep7: -0.3485396 < x < -0.2156404
c          -0.110 < y < 0.110
c     q1ex is a circle of radius 0.1492 m
c     dent is a trapazoid:
c                                   -5.22008 < x < -4.98099
c             -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784  
c     dext is also a trapazoid: 
c                                   -0.46188 < x < 0.46188
c                   -(-0.01610808*x + 0.125) < y < -0.01610808*x + 0.125
c     q3en is a circle of radius 0.3 m
c     q3ex is a circle of radius 0.3 m
c


      function x_sr_ep3     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.1164706E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29922E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25726720E-02, 0.32485950E-01, 0.75389570E-06, 0.19876383E-01,
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
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      x_sr_ep3     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x33        
     4  +coeff(  4)        *x31        
c
      return
      end
      function t_sr_ep3     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.1069279E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29922E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11965680E-02,-0.43500770E-06, 0.13630892E-04, 0.29982750E-01,
     +  0.44200543E-04, 0.38676461E-03,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      t_sr_ep3     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
c
      return
      end
      function y_sr_ep3     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/  0.2605673E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29922E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.27264560E-03,-0.59839100E-01,-0.35577011E-06, 0.18868391E-06,
     + -0.19964510E-02,-0.60292681E-04,
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
      x41 = x4
c
c                  function
c
      y_sr_ep3     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x23            
c
      return
      end
      function p_sr_ep3     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/  0.2267454E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29922E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23819940E-03,-0.54912430E-01, 0.19500340E-03, 0.15295280E-03,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
c
c                  function
c
      p_sr_ep3     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23*x31        
c
      return
      end
      function l_sr_ep3     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/  0.1089599E+01/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29922E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.70162540E-03, 0.16447261E-02, 0.52682054E-03,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_sr_ep3     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x22            
     3  +coeff(  3)            *x42    
c
      return
      end
      function x_sr_ep4     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.1415698E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29580480E-02, 0.38780430E-01, 0.19813382E-01,
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
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x31 = x3
      x41 = x4
c
c                  function
c
      x_sr_ep4     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31        
c
      return
      end
      function t_sr_ep4     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.1268557E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15345194E-02,-0.13169590E-04, 0.19664601E-03, 0.30169720E-01,
     +  0.10404641E-02,-0.29002650E-03, 0.11662253E-04,-0.10712320E-03,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x51 = x5
c
c                  function
c
      t_sr_ep4     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
c
      return
      end
      function y_sr_ep4     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/  0.2217465E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23665800E-03,-0.71826520E-01, 0.56318140E-07, 0.26447780E-06,
     + -0.19940764E-02, 0.20651500E-03, 0.95339082E-05, 0.98902620E-04,
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
      x41 = x4
c
c                  function
c
      y_sr_ep4     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23*x31        
     8  +coeff(  8)    *x21*x31        
c
      return
      end
      function p_sr_ep4     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/  0.1496922E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16031574E-03,-0.54008760E-01, 0.62500410E-03, 0.13130090E-02,
     +  0.28057102E-03, 0.36220284E-03,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      p_sr_ep4     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21*x31        
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)    *x23*x31*x41    
c
      return
      end
      function l_sr_ep4     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/  0.1308974E+01/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.64883050E-03, 0.10586522E-02, 0.19794914E-02,-0.39006231E-03,
     +  0.29099650E-03, 0.68285514E-03,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      l_sr_ep4     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22            
     4  +coeff(  4)        *x31*x41    
     5  +coeff(  5)    *x22*x31        
     6  +coeff(  6)        *x33        
c
      return
      end
      function x_sr_ep5     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.1722790E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.34938312E-02, 0.44696263E-01, 0.19683990E-01,
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
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x31 = x3
      x41 = x4
c
c                  function
c
      x_sr_ep5     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31        
c
      return
      end
      function t_sr_ep5     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.1636614E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21210840E-02,-0.23799480E-04,-0.46609140E-04, 0.29467193E-01,
     +  0.28855800E-02,-0.87882754E-05,-0.15406074E-02, 0.43505911E-04,
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
      x41 = x4
      x51 = x5
c
c                  function
c
      t_sr_ep5     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
c
      return
      end
      function y_sr_ep5     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/  0.2530935E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.26953584E-03,-0.83175830E-01,-0.45140510E-06, 0.15610883E-05,
     + -0.19908980E-02, 0.34872000E-03, 0.10074813E-04, 0.15939974E-03,
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
      x41 = x4
c
c                  function
c
      y_sr_ep5     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23*x31        
     8  +coeff(  8)    *x21*x31        
c
      return
      end
      function p_sr_ep5     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/  0.1507681E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16207400E-03,-0.53917821E-01,-0.20272022E-05, 0.49577991E-03,
     +  0.49906870E-03,-0.47538971E-04, 0.41928260E-03, 0.15201260E-02,
     +  0.12505110E-02, 0.22944050E-02, 0.46725091E-02, 0.43516140E-02,
     +  0.21824940E-02, 0.47294432E-02, 0.29093330E-02, 0.13238752E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      p_sr_ep5     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21*x32        
     8  +coeff(  8)    *x21*x31*x41    
      p_sr_ep5     =p_sr_ep5     
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x21*x32*x41    
     2  +coeff( 11)    *x21*x31*x42    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)    *x21*x32*x42    
     5  +coeff( 14)    *x21*x31*x43    
     6  +coeff( 15)    *x21    *x44    
     7  +coeff( 16)    *x21*x31*x44    
c
      return
      end
      function l_sr_ep5     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/  0.1519782E+01/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29732E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.54557324E-03, 0.33695092E-02, 0.22884630E-02, 0.88453580E-03,
     +  0.19322330E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      l_sr_ep5     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x22*x31        
     5  +coeff(  5)        *x33        
c
      return
      end
      function x_sr_ep6     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.2110171E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35949250E-02, 0.49929410E-01, 0.19489790E-01,
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
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x31 = x3
      x41 = x4
c
c                  function
c
      x_sr_ep6     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31        
c
      return
      end
      function t_sr_ep6     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.2007603E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25733170E-02,-0.81933630E-04,-0.45774950E-03, 0.29046530E-01,
     +  0.46237830E-02,-0.12769383E-02,-0.47392424E-03, 0.66262400E-04,
     + -0.13738660E-02,-0.46794560E-02,-0.10728694E-01,-0.87227150E-02,
     + -0.29809151E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      t_sr_ep6     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x33*x41    
      t_sr_ep6     =t_sr_ep6     
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x24*x31*x41    
     4  +coeff( 13)    *x24*x33        
c
      return
      end
      function y_sr_ep6     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 20)
      data ncoeff/ 19/
      data avdat/  0.2892061E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.30865200E-03,-0.94699442E-01, 0.46971020E-06, 0.48645442E-05,
     + -0.19834480E-02, 0.48355000E-03, 0.10524840E-02, 0.44425243E-04,
     +  0.17614143E-03,-0.55013670E-03, 0.12433550E-03, 0.15002170E-03,
     +  0.67469052E-03, 0.86995562E-04,-0.10277700E-04, 0.13386090E-02,
     +  0.85970730E-03, 0.21762274E-03,-0.39471942E-03,
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
      y_sr_ep6     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31        
      y_sr_ep6     =y_sr_ep6     
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)    *x21*x33*x42    
     2  +coeff( 11)    *x23    *x43    
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x21*x33        
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)    *x21    *x43    
      y_sr_ep6     =y_sr_ep6     
     9  +coeff( 18)    *x21*x34        
     1  +coeff( 19)    *x21*x32*x42    
c
      return
      end
      function p_sr_ep6     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/  0.1641316E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.17517313E-03,-0.54103010E-01,-0.34445790E-05, 0.13715540E-04,
     +  0.12104370E-03,-0.75596430E-03, 0.25957204E-04,-0.42141910E-04,
     + -0.71577363E-04, 0.56062132E-03, 0.22496990E-02,-0.41677160E-04,
     +  0.15101630E-02,-0.11270600E-04, 0.39476351E-02, 0.83422030E-02,
     + -0.60501050E-04, 0.84579354E-02, 0.50489143E-02, 0.93695390E-02,
     +  0.65457530E-02, 0.30863800E-03, 0.44516450E-04, 0.26603310E-02,
     +  0.12761140E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sr_ep6     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sr_ep6     =p_sr_ep6     
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)        *x32*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)    *x21*x32*x41    
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)        *x32*x42    
      p_sr_ep6     =p_sr_ep6     
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x21*x32*x42    
     2  +coeff( 20)    *x21*x31*x43    
     3  +coeff( 21)    *x21    *x44    
     4  +coeff( 22)    *x23*x33        
     5  +coeff( 23)    *x21*x32*x43    
     6  +coeff( 24)    *x21*x31*x44    
     7  +coeff( 25)    *x21*x33*x43    
c
      return
      end
      function l_sr_ep6     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/  0.1732274E+01/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.52108831E-03, 0.54939980E-02, 0.25783370E-02, 0.13462450E-02,
     +  0.28910690E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      l_sr_ep6     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22            
     4  +coeff(  4)    *x22*x31        
     5  +coeff(  5)        *x33        
c
      return
      end
      function x_sr_ep7     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.2604133E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.45087570E-02, 0.54825901E-01, 0.78956620E-03, 0.24934150E-03,
     + -0.90968930E-03, 0.18677140E-01, 0.23844621E-02,-0.38922110E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
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
      x_sr_ep7     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22*x31        
     4  +coeff(  4)        *x33        
     5  +coeff(  5)    *x22*x33        
     6  +coeff(  6)        *x31        
     7  +coeff(  7)                *x51
     8  +coeff(  8)    *x22    *x42    
c
      return
      end
      function t_sr_ep7     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.2220073E+00/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.31429620E-02,-0.42625630E-04, 0.31703170E-05, 0.29546981E-01,
     +  0.56824300E-02,-0.20011013E-02, 0.41428440E-04,-0.65158060E-02,
     + -0.10814744E-01,-0.14070883E-01,-0.25723173E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      t_sr_ep7     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x22    *x41    
      t_sr_ep7     =t_sr_ep7     
     9  +coeff(  9)    *x22*x31*x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x22*x31        
c
      return
      end
      function y_sr_ep7     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/  0.4561184E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.47642710E-03,-0.10728790E+00, 0.30013690E-05,-0.19756610E-02,
     + -0.16721343E-03, 0.53434702E-03, 0.21754880E-02, 0.39046080E-04,
     + -0.35328532E-05, 0.89380494E-03, 0.27535611E-02,-0.22835130E-04,
     + -0.10193834E-02,-0.19453621E-03,-0.14468390E-02, 0.15013440E-03,
     +  0.19373423E-03, 0.30433392E-03,-0.28372242E-04, 0.75604191E-04,
     + -0.54864543E-04, 0.23967630E-02, 0.30857560E-02, 0.82873490E-03,
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
      y_sr_ep7     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21*x32        
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31        
      y_sr_ep7     =y_sr_ep7     
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x21*x32*x41    
     2  +coeff( 11)    *x21    *x43    
     3  +coeff( 12)    *x23*x31*x41    
     4  +coeff( 13)    *x21*x33*x41    
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)    *x21*x32*x42    
     7  +coeff( 16)    *x23*x31*x42    
     8  +coeff( 17)    *x21*x33*x42    
      y_sr_ep7     =y_sr_ep7     
     9  +coeff( 18)    *x23*x34*x41    
     1  +coeff( 19)    *x23*x34*x42    
     2  +coeff( 20)    *x21*x31        
     3  +coeff( 21)    *x21        *x51
     4  +coeff( 22)    *x21*x31*x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x44    
c
      return
      end
      function p_sr_ep7     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.2245670E-03/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23666140E-03,-0.53710170E-01, 0.21232370E-04, 0.37418340E-05,
     +  0.23260770E-03,-0.46780082E-03, 0.41252322E-04,-0.10180024E-03,
     +  0.16137411E-04,-0.50375281E-04, 0.81289804E-03, 0.32430970E-02,
     + -0.67656663E-04, 0.31188153E-02,-0.41030380E-04,-0.26458974E-04,
     + -0.13430021E-05,-0.31441260E-03, 0.80059384E-04, 0.39861993E-02,
     +  0.10717120E-01,-0.77098630E-04, 0.86225224E-02,-0.53018800E-04,
     +  0.37536950E-02, 0.87895030E-02, 0.58236964E-02, 0.19951321E-03,
     +  0.58506301E-03,
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
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sr_ep7     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sr_ep7     =p_sr_ep7     
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)        *x32*x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)        *x31*x41*x51
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)                *x53
      p_sr_ep7     =p_sr_ep7     
     9  +coeff( 18)    *x23*x31        
     1  +coeff( 19)*x11    *x32        
     2  +coeff( 20)    *x21*x32*x41    
     3  +coeff( 21)    *x21*x31*x42    
     4  +coeff( 22)        *x32*x42    
     5  +coeff( 23)    *x21    *x43    
     6  +coeff( 24)    *x21    *x41*x52
     7  +coeff( 25)    *x21*x32*x42    
     8  +coeff( 26)    *x21*x31*x43    
      p_sr_ep7     =p_sr_ep7     
     9  +coeff( 27)    *x21    *x44    
     1  +coeff( 28)    *x23    *x41*x51
     2  +coeff( 29)    *x23*x33        
c
      return
      end
      function l_sr_ep7     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/  0.1962772E+01/
      data xmin/
     1 -0.19975E-02,-0.54936E-01,-0.19979E-01,-0.29994E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54959E-01, 0.19967E-01, 0.29220E-01, 0.49960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.78466650E-03, 0.10955550E-01, 0.20659870E-02, 0.52067660E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
c
c                  function
c
      l_sr_ep7     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x22*x31        
     4  +coeff(  4)        *x33        
c
      return
      end
      function x_sr_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.2574444E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30085991E-02, 0.14510660E+00,-0.35522150E-02, 0.32514860E-02,
     + -0.60296720E-02,-0.19442860E-02,-0.73312770E-02,-0.77480940E-02,
     +  0.10292650E-02,
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
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      x_sr_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x23    *x42    
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)    *x21    *x42    
      x_sr_q1ex    =x_sr_q1ex    
     9  +coeff(  9)*x11                
c
      return
      end
      function t_sr_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.1792072E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21492030E-03,-0.91769583E-02,-0.82579744E-03,-0.32899670E-03,
     +  0.29217420E-02,-0.31644100E-02,-0.27577872E-03,-0.22071290E-02,
     + -0.34522500E-03,-0.37294061E-03,-0.20241080E-02,-0.15643793E-02,
     + -0.16467650E-02,
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
      t_sr_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21*x31*x41    
      t_sr_q1ex    =t_sr_q1ex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21    *x43    
     3  +coeff( 12)    *x23*x31*x42    
     4  +coeff( 13)    *x21*x33*x42    
c
      return
      end
      function y_sr_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 43)
      data ncoeff/ 42/
      data avdat/ -0.8513198E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.44633880E-02,-0.36919920E-04, 0.30132240E-01, 0.15707650E+00,
     +  0.17367031E-01,-0.68896161E-02, 0.41438190E-03,-0.19374421E-02,
     + -0.64761172E-02, 0.11590830E-03,-0.88661123E-03,-0.96695290E-03,
     + -0.90681420E-02, 0.11994220E-02, 0.11783424E-01, 0.66695640E-04,
     + -0.38983260E-01, 0.13444874E-02,-0.12198751E+00,-0.13163261E+00,
     + -0.27367090E-01,-0.10684760E-02,-0.40516260E-02,-0.47314334E-01,
     + -0.95891490E-01, 0.11629020E-01,-0.62598690E-03,-0.51901913E-02,
     + -0.22618500E-01, 0.70611713E-03,-0.25704770E-01,-0.93377614E-02,
     +  0.20769860E-01,-0.51547610E-02, 0.12818680E+00,-0.58523672E-02,
     + -0.17220070E-01,-0.56133340E-02, 0.10203640E+00, 0.55143260E-01,
     + -0.60996674E-01, 0.17670660E-01,
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
      y_sr_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)            *x41*x51
      y_sr_q1ex    =y_sr_q1ex    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)        *x32*x41    
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x22*x31*x41    
     5  +coeff( 14)        *x33*x41    
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)    *x24*x31        
     8  +coeff( 17)    *x22    *x43    
      y_sr_q1ex    =y_sr_q1ex    
     9  +coeff( 18)    *x24*x32        
     1  +coeff( 19)    *x22*x31*x43    
     2  +coeff( 20)    *x22    *x44    
     3  +coeff( 21)    *x22*x31*x44    
     4  +coeff( 22)                *x52
     5  +coeff( 23)    *x22*x32        
     6  +coeff( 24)    *x22*x31*x42    
     7  +coeff( 25)    *x22*x32*x42    
     8  +coeff( 26)    *x24*x32*x41    
      y_sr_q1ex    =y_sr_q1ex    
     9  +coeff( 27)        *x31    *x51
     1  +coeff( 28)    *x22*x31        
     2  +coeff( 29)    *x22*x32*x41    
     3  +coeff( 30)    *x22    *x42*x51
     4  +coeff( 31)    *x22*x33*x41    
     5  +coeff( 32)    *x22    *x43*x51
     6  +coeff( 33)    *x24*x31*x42    
     7  +coeff( 34)*x11*x23    *x41*x51
     8  +coeff( 35)    *x22*x32*x44    
      y_sr_q1ex    =y_sr_q1ex    
     9  +coeff( 36)*x11*x24    *x41*x51
     1  +coeff( 37)    *x22*x34*x43    
     2  +coeff( 38)*x12*x23    *x41*x51
     3  +coeff( 39)    *x24*x33*x43    
     4  +coeff( 40)    *x24*x32*x44    
     5  +coeff( 41)    *x22*x34*x44    
     6  +coeff( 42)*x11*x23*x33*x43    
c
      return
      end
      function p_sr_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 18)
      data ncoeff/ 17/
      data avdat/ -0.3824610E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.18801010E-02, 0.92157530E-02, 0.68227440E-01, 0.92090950E-02,
     + -0.27860074E-02, 0.14326070E-03,-0.82084452E-02,-0.84861350E-02,
     +  0.87969694E-02,-0.22103330E-02,-0.37900132E-02,-0.14807210E-01,
     + -0.17869620E-01,-0.74947632E-01,-0.46234980E-01, 0.17540090E-03,
     + -0.44643480E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      p_sr_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x22    *x42    
      p_sr_q1ex    =p_sr_q1ex    
     9  +coeff(  9)    *x24*x31*x41    
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x22*x31*x41    
     4  +coeff( 13)    *x22    *x44    
     5  +coeff( 14)    *x24*x31*x43    
     6  +coeff( 15)    *x24    *x44    
     7  +coeff( 16)        *x32        
     8  +coeff( 17)    *x22*x32        
c
      return
      end
      function l_sr_q1ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/  0.3943433E+01/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19958E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16952070E-02, 0.21620572E-02, 0.48303930E-02, 0.38695900E-02,
     +  0.33554554E-03, 0.33365320E-02, 0.74823520E-03,-0.54391071E-03,
     + -0.26938480E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      l_sr_q1ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22    *x41    
      l_sr_q1ex    =l_sr_q1ex    
     9  +coeff(  9)    *x24    *x44    
c
      return
      end
      function x_sr_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.5098287E+01/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21937360E-02,-0.12602543E+00, 0.34738490E-03,-0.69291090E-02,
     +  0.21979310E-01, 0.36465651E-02, 0.11515430E-02, 0.16899714E-01,
     +  0.14750410E-02, 0.32325820E-02, 0.28014670E-01, 0.25127730E-01,
     +  0.95651233E-02,
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
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      x_sr_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x21    *x42    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21*x31*x41    
      x_sr_dent    =x_sr_dent    
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31*x42    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)    *x21*x34*x41    
c
      return
      end
      function t_sr_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/  0.1298947E+01/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.53297064E-03, 0.55919900E-01,-0.55810520E-02, 0.35877470E-02,
     +  0.30242344E-02,-0.73934951E-02, 0.73339883E-03,-0.89417460E-03,
     + -0.16815784E-02,-0.28431040E-02,-0.16110253E-03,-0.93725612E-02,
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
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      t_sr_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x23    *x42    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21*x31*x41    
     7  +coeff(  7)                *x51
     8  +coeff(  8)*x11                
      t_sr_dent    =t_sr_dent    
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x42    
c
      return
      end
      function y_sr_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 81)
      data ncoeff/ 80/
      data avdat/ -0.5781254E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.33770443E-02,-0.53993810E-03, 0.66151130E-02, 0.11161850E+00,
     +  0.17902344E-01,-0.70137400E-02, 0.28497680E-03, 0.56922701E-02,
     +  0.29516974E-02, 0.25242790E-02, 0.16551030E-01, 0.12690121E-02,
     + -0.13261680E-01,-0.22553620E-02,-0.33018660E-02,-0.36436773E-02,
     + -0.14762711E-02,-0.22325640E-01,-0.73844194E-02,-0.74120662E-02,
     + -0.21701222E-01,-0.49779582E-02, 0.14263550E-02, 0.22941510E-03,
     + -0.82602490E-01, 0.26137900E-01,-0.95567882E-01, 0.19965292E-01,
     + -0.91685093E-02, 0.11017541E-01,-0.15123290E+00,-0.26528320E-02,
     +  0.10629723E-02, 0.13696890E-02, 0.38959140E-02,-0.74797272E-02,
     + -0.14541460E-02,-0.54806380E-02,-0.16219930E-01,-0.92123580E-01,
     + -0.17839441E-01,-0.89974940E-01,-0.77013190E-01,-0.13390061E-01,
     +  0.59387722E-03,-0.31293820E-03,-0.30689220E-03,-0.12907792E-02,
     + -0.36974630E-02, 0.24035380E-03,-0.18573520E-02,-0.60852400E-03,
     +  0.30289310E-02, 0.47103760E-03, 0.24000970E-01, 0.39845321E-01,
     + -0.20750374E-02, 0.46576750E-02,-0.40387650E-02,-0.19914310E-01,
     + -0.16899921E-01, 0.66932400E-01, 0.25745402E-03,-0.12524850E-01,
     +  0.57548941E-02, 0.47138541E-01, 0.22995094E-02,-0.64773983E-02,
     + -0.94764650E-01, 0.94450123E-01,-0.11453174E-01, 0.10085950E-01,
     +  0.74121281E-01,-0.28021184E-02, 0.11840360E+00,-0.10894041E-01,
     + -0.18348911E-01, 0.11603950E-01,-0.36307130E-01, 0.11875580E+00,
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
      y_sr_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      y_sr_dent    =y_sr_dent    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)        *x32*x41    
     6  +coeff( 15)        *x31*x42    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)            *x42*x51
      y_sr_dent    =y_sr_dent    
     9  +coeff( 18)    *x22*x31*x41    
     1  +coeff( 19)        *x33*x41    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)        *x32*x42    
     4  +coeff( 22)            *x44    
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)    *x24*x31        
     7  +coeff( 25)    *x22    *x43    
     8  +coeff( 26)    *x24    *x42    
      y_sr_dent    =y_sr_dent    
     9  +coeff( 27)    *x22    *x44    
     1  +coeff( 28)        *x32*x44    
     2  +coeff( 29)    *x24*x31*x42    
     3  +coeff( 30)    *x22*x32*x43    
     4  +coeff( 31)    *x22*x31*x44    
     5  +coeff( 32)    *x22    *x44*x51
     6  +coeff( 33)    *x24*x34        
     7  +coeff( 34)    *x21*x31        
     8  +coeff( 35)        *x31*x41    
      y_sr_dent    =y_sr_dent    
     9  +coeff( 36)    *x22*x31        
     1  +coeff( 37)        *x31*x41*x51
     2  +coeff( 38)    *x22*x32        
     3  +coeff( 39)        *x31*x43    
     4  +coeff( 40)    *x22*x31*x43    
     5  +coeff( 41)    *x22    *x43*x51
     6  +coeff( 42)    *x24*x32*x42    
     7  +coeff( 43)    *x24    *x44    
     8  +coeff( 44)*x12*x21    *x43*x51
      y_sr_dent    =y_sr_dent    
     9  +coeff( 45)    *x21        *x51
     1  +coeff( 46)*x11*x21            
     2  +coeff( 47)        *x33        
     3  +coeff( 48)    *x22        *x51
     4  +coeff( 49)            *x41*x52
     5  +coeff( 50)*x11*x21        *x51
     6  +coeff( 51)        *x31*x42*x51
     7  +coeff( 52)*x11*x22*x31        
     8  +coeff( 53)            *x41*x54
      y_sr_dent    =y_sr_dent    
     9  +coeff( 54)*x11*x22*x31*x41    
     1  +coeff( 55)        *x34*x42    
     2  +coeff( 56)        *x33*x43    
     3  +coeff( 57)    *x24*x31    *x51
     4  +coeff( 58)*x11*x22    *x41*x51
     5  +coeff( 59)*x11*x21    *x42*x51
     6  +coeff( 60)    *x24*x32*x41    
     7  +coeff( 61)    *x22*x33*x42    
     8  +coeff( 62)    *x24    *x43    
      y_sr_dent    =y_sr_dent    
     9  +coeff( 63)        *x34*x43    
     1  +coeff( 64)*x11*x23    *x41*x51
     2  +coeff( 65)*x11*x21*x32*x41*x51
     3  +coeff( 66)    *x22    *x43*x52
     4  +coeff( 67)    *x23*x31    *x53
     5  +coeff( 68)    *x22    *x41*x54
     6  +coeff( 69)    *x22*x34*x42    
     7  +coeff( 70)    *x22*x32*x44    
     8  +coeff( 71)*x11*x24    *x41*x51
      y_sr_dent    =y_sr_dent    
     9  +coeff( 72)    *x22*x34*x43    
     1  +coeff( 73)    *x24*x31*x44    
     2  +coeff( 74)*x11*x21*x32*x44    
     3  +coeff( 75)    *x22*x33*x44    
     4  +coeff( 76)*x12*x21*x31*x42*x51
     5  +coeff( 77)    *x22*x32*x43*x52
     6  +coeff( 78)*x11*x23    *x41*x53
     7  +coeff( 79)    *x22    *x43*x54
     8  +coeff( 80)    *x24*x34*x42    
c
      return
      end
      function p_sr_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.1392198E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.71128120E-03,-0.67318151E-02,-0.23744640E-01,-0.18127660E-02,
     +  0.85982480E-03,-0.35031030E-02, 0.44305240E-02, 0.73141430E-04,
     +  0.20338030E-03,-0.11086163E-04,-0.46981751E-03, 0.70281820E-03,
     +  0.62512060E-03, 0.32635760E-02, 0.99758581E-04,-0.16703880E-02,
     +  0.33586603E-03,-0.69298333E-03, 0.18663860E-02, 0.30733920E-02,
     +  0.18149350E-02, 0.85733542E-02, 0.84199090E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sr_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x24    *x42    
      p_sr_dent    =p_sr_dent    
     9  +coeff(  9)    *x21            
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x22*x31*x41    
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21    *x42    
      p_sr_dent    =p_sr_dent    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)    *x22    *x44    
     4  +coeff( 22)    *x24*x31*x42    
     5  +coeff( 23)    *x24    *x43    
c
      return
      end
      function l_sr_dent    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/  0.1075976E+02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.71754311E-02,-0.23651260E+00, 0.16966061E-01,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
c
c                  function
c
      l_sr_dent    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x22            
c
      return
      end
      function x_sr_dext     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 35)
      data ncoeff/ 34/
      data avdat/ -0.2336794E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.40092320E-02, 0.34526693E+00, 0.13303080E+00, 0.66025620E-02,
     +  0.45855360E-01,-0.36773380E-01, 0.75839660E-02, 0.10891420E-01,
     + -0.10282120E-02,-0.29183880E-01,-0.54791810E-02,-0.33827380E-02,
     + -0.76843863E-02,-0.82020033E-02,-0.81015750E-01,-0.93387514E-01,
     + -0.84693040E-03,-0.39410950E-02,-0.48579390E-02,-0.11765431E-02,
     + -0.10713454E-02,-0.19210143E-02,-0.36073650E-01,-0.51151350E-02,
     + -0.11977761E+00,-0.95100753E-01,-0.19154380E-02, 0.20599694E-01,
     + -0.22464390E-01,-0.54260270E-01,-0.51626130E-02,-0.34355000E-01,
     + -0.62016160E-01,-0.36275330E-01,
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
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sr_dext     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x23*x31*x41    
     8  +coeff(  8)    *x22            
      x_sr_dext     =x_sr_dext     
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21*x31*x42    
     7  +coeff( 16)    *x21    *x43    
     8  +coeff( 17)    *x21*x34*x41    
      x_sr_dext     =x_sr_dext     
     9  +coeff( 18)            *x42    
     1  +coeff( 19)    *x21*x33        
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)            *x41*x51
     4  +coeff( 22)    *x21*x31    *x51
     5  +coeff( 23)    *x21*x32*x41    
     6  +coeff( 24)    *x21    *x42*x51
     7  +coeff( 25)    *x21*x31*x43    
     8  +coeff( 26)    *x21    *x44    
      x_sr_dext     =x_sr_dext     
     9  +coeff( 27)    *x23        *x52
     1  +coeff( 28)    *x23    *x43    
     2  +coeff( 29)    *x21*x32*x43    
     3  +coeff( 30)    *x21*x31*x44    
     4  +coeff( 31)    *x21*x31*x41*x53
     5  +coeff( 32)    *x23*x33*x41    
     6  +coeff( 33)    *x23*x32*x42    
     7  +coeff( 34)    *x21*x34*x42    
c
      return
      end
      function t_sr_dext     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/  0.5336839E+00/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.53904360E-03,-0.78800424E-01, 0.22719563E-01,-0.33301864E-02,
     +  0.40684253E-02,-0.38496410E-02, 0.14043630E-01, 0.88294290E-03,
     +  0.22366500E-02,-0.12483622E-02, 0.87312562E-02,
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
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sr_dext     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)*x11                
      t_sr_dext     =t_sr_dext     
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21*x31*x41    
c
      return
      end
      function y_sr_dext     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.4796520E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.55530230E-03,-0.36494333E-01, 0.13115750E-01, 0.17326800E-01,
     + -0.70642520E-02,-0.52925031E-02, 0.77735880E-03,-0.51475290E-01,
     +  0.20992590E-02, 0.16386650E-03,-0.68741482E-02, 0.88135283E-02,
     +  0.54991020E-01, 0.59559860E-02, 0.37386552E-02,-0.94941280E-02,
     + -0.84876391E-03,-0.41887111E-01,-0.30594242E-02, 0.59717064E-02,
     + -0.63011883E-02,-0.81857964E-02,-0.16247543E-02,-0.40236022E-02,
     + -0.53426390E-02,-0.47146340E-02, 0.24162363E-02,-0.14003442E-01,
     +  0.62645480E-02, 0.13198540E-01,-0.79383890E-02, 0.10903120E-01,
     +  0.15805703E-02,-0.45070740E-02,-0.17060970E-02, 0.71103712E-02,
     +  0.69973524E-02,-0.14088303E-01, 0.74969190E-03,-0.65403774E-01,
     + -0.15634470E-02,-0.87189972E-02, 0.25956850E-02, 0.11222950E-01,
     +  0.37422490E-02, 0.24871290E-01, 0.32116010E-01,-0.12568610E+00,
     + -0.39474400E-01, 0.46223890E-01, 0.17731593E-02, 0.76051882E-03,
     +  0.50001140E-02,-0.73558674E-02,-0.10861502E-02,-0.10591790E-02,
     + -0.10359772E-02, 0.18380692E-02,-0.36019750E-02, 0.31306530E-03,
     + -0.41466820E-02,-0.32537630E-02, 0.23359990E-02,-0.78247264E-02,
     +  0.25247340E-02,-0.33821470E-02, 0.42612560E-02,-0.10517240E+00,
     + -0.36064420E-01, 0.62812264E-02,-0.16217300E-01,-0.22357760E-02,
     +  0.65760140E-01,-0.80917462E-01, 0.18220524E-02,-0.57132500E-02,
     + -0.10121642E+00,-0.80271914E-01, 0.20550280E-01, 0.11642792E+00,
     +  0.28414253E-02,-0.58733550E-02, 0.52988223E-01, 0.15439870E+00,
     + -0.50604690E-01,-0.69869750E-01, 0.18300880E-01,-0.10741911E+00,
     + -0.12114670E+00,-0.28346033E-03, 0.20806463E-02,-0.75945150E-02,
     +  0.16358693E-02,-0.32721790E-02, 0.82275290E-03,
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
      x34 = x33*x3
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
      y_sr_dext     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      y_sr_dext     =y_sr_dext     
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)        *x33        
      y_sr_dext     =y_sr_dext     
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)        *x32*x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)            *x43    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x21*x31    *x51
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)            *x41*x52
      y_sr_dext     =y_sr_dext     
     9  +coeff( 27)    *x23    *x41    
     1  +coeff( 28)    *x22*x31*x41    
     2  +coeff( 29)    *x21*x32*x41    
     3  +coeff( 30)    *x22    *x42    
     4  +coeff( 31)        *x32*x42    
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)            *x44    
     7  +coeff( 34)    *x22    *x41*x51
     8  +coeff( 35)            *x43*x51
      y_sr_dext     =y_sr_dext     
     9  +coeff( 36)    *x23*x31*x41    
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x22*x31*x42    
     3  +coeff( 39)        *x33*x42    
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)        *x32*x43    
     6  +coeff( 42)    *x22*x31*x41*x51
     7  +coeff( 43)        *x33*x41*x51
     8  +coeff( 44)    *x22    *x42*x51
      y_sr_dext     =y_sr_dext     
     9  +coeff( 45)    *x24*x32        
     1  +coeff( 46)    *x23*x31*x42    
     2  +coeff( 47)    *x23    *x43    
     3  +coeff( 48)    *x22    *x44    
     4  +coeff( 49)    *x22    *x44*x51
     5  +coeff( 50)    *x24*x31*x43    
     6  +coeff( 51)    *x21*x32        
     7  +coeff( 52)*x11        *x41    
     8  +coeff( 53)    *x21*x31*x41    
      y_sr_dext     =y_sr_dext     
     9  +coeff( 54)        *x31*x42    
     1  +coeff( 55)        *x31*x41*x51
     2  +coeff( 56)        *x31    *x52
     3  +coeff( 57)                *x53
     4  +coeff( 58)    *x24            
     5  +coeff( 59)        *x33*x41    
     6  +coeff( 60)*x11        *x42    
     7  +coeff( 61)        *x31*x43    
     8  +coeff( 62)    *x22*x31    *x51
      y_sr_dext     =y_sr_dext     
     9  +coeff( 63)        *x31*x44    
     1  +coeff( 64)        *x31*x43*x51
     2  +coeff( 65)    *x21*x31    *x53
     3  +coeff( 66)    *x24*x31*x41    
     4  +coeff( 67)    *x23*x32*x41    
     5  +coeff( 68)    *x22*x31*x43    
     6  +coeff( 69)    *x21*x32*x43    
     7  +coeff( 70)    *x23    *x42*x51
     8  +coeff( 71)    *x22    *x43*x51
      y_sr_dext     =y_sr_dext     
     9  +coeff( 72)        *x31*x44*x51
     1  +coeff( 73)    *x24    *x43    
     2  +coeff( 74)    *x22*x31*x44    
     3  +coeff( 75)*x12*x22        *x51
     4  +coeff( 76)*x12*x21    *x41*x51
     5  +coeff( 77)    *x24*x32*x42    
     6  +coeff( 78)    *x22*x34*x42    
     7  +coeff( 79)    *x21*x34*x43    
     8  +coeff( 80)    *x22*x32*x44    
      y_sr_dext     =y_sr_dext     
     9  +coeff( 81)*x12*x21*x31*x41*x51
     1  +coeff( 82)    *x24*x34*x41    
     2  +coeff( 83)    *x24*x31*x44    
     3  +coeff( 84)    *x24*x34*x42    
     4  +coeff( 85)    *x22*x34*x44    
     5  +coeff( 86)    *x24*x34*x43    
     6  +coeff( 87)*x12*x23*x32*x41*x51
     7  +coeff( 88)    *x24*x34*x42*x51
     8  +coeff( 89)    *x24*x33*x43*x51
      y_sr_dext     =y_sr_dext     
     9  +coeff( 90)*x11*x21            
     1  +coeff( 91)    *x23*x31        
     2  +coeff( 92)    *x22*x32        
     3  +coeff( 93)    *x21*x31*x41*x51
     4  +coeff( 94)        *x31*x42*x51
     5  +coeff( 95)    *x21*x31    *x52
c
      return
      end
      function p_sr_dext     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.9551189E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.37555190E-03,-0.92669430E-02,-0.13063563E-01, 0.12240763E-02,
     + -0.37229070E-03,-0.73650590E-03,-0.56217580E-04,-0.97131934E-02,
     + -0.28184714E-03,-0.12421134E-02, 0.11366330E-02, 0.99695594E-02,
     +  0.12707200E-02, 0.91840770E-03,-0.68805470E-02, 0.18410441E-02,
     + -0.17953274E-02,-0.60153024E-03, 0.15872140E-02,-0.87826041E-03,
     + -0.11972230E-02, 0.24790011E-03,-0.22137491E-02, 0.72861790E-02,
     +  0.69327550E-02,-0.29661363E-03,-0.86757872E-03,-0.77018390E-03,
     + -0.41799770E-03,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
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
      p_sr_dext     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      p_sr_dext     =p_sr_dext     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x22        *x51
      p_sr_dext     =p_sr_dext     
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)            *x42*x51
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)    *x23    *x43    
     7  +coeff( 25)    *x23*x31*x44    
     8  +coeff( 26)        *x31*x41    
      p_sr_dext     =p_sr_dext     
     9  +coeff( 27)            *x43    
     1  +coeff( 28)        *x31*x41*x51
     2  +coeff( 29)            *x41*x52
c
      return
      end
      function l_sr_dext     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/  0.1734885E+02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21031250E-02, 0.49520020E+00, 0.98460420E-01, 0.31669530E-01,
     + -0.16785961E-02, 0.44648122E-01,-0.10148034E+00, 0.14825254E-01,
     + -0.74184370E-02,-0.47976700E-02,-0.80369800E-01,-0.13958903E-01,
     + -0.88777740E-02,-0.11788070E+00,-0.10050160E+00,-0.44408261E-01,
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
      x51 = x5
c
c                  function
c
      l_sr_dext     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31*x41    
      l_sr_dext     =l_sr_dext     
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x21*x31*x42    
     6  +coeff( 15)    *x21    *x43    
     7  +coeff( 16)    *x21*x34*x41    
c
      return
      end
      function x_sr_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 39)
      data ncoeff/ 38/
      data avdat/ -0.5595049E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.80596350E-02, 0.21889270E+00,-0.37562483E-03, 0.13866850E+00,
     +  0.18050160E-01, 0.22287773E-02, 0.35813421E-01,-0.23674300E-01,
     +  0.24772910E-02,-0.11296550E-03,-0.58004893E-02,-0.18413914E-01,
     + -0.38600354E-02,-0.42140200E-02,-0.44866753E-02,-0.53588070E-02,
     + -0.67028291E-01,-0.58318010E-01, 0.10215210E-01,-0.16491230E-02,
     + -0.92757411E-01,-0.62358530E-01,-0.26809473E-03,-0.10997240E-02,
     +  0.22667541E-02, 0.22727470E-02,-0.18355700E-02,-0.39635080E-02,
     +  0.63137500E-02,-0.34202670E-01,-0.17277203E-02,-0.45775990E-02,
     + -0.10158820E-01,-0.50502900E-01, 0.10604600E-01,-0.26793092E-01,
     + -0.22188110E-01, 0.15742640E-01,
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
      x24 = x23*x2
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
      x_sr_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x21    *x42    
      x_sr_q3en   =x_sr_q3en   
     9  +coeff(  9)    *x23*x31*x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)*x11                
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x21*x31*x42    
      x_sr_q3en   =x_sr_q3en   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x21*x34*x41    
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)    *x21*x31*x43    
     4  +coeff( 22)    *x21    *x44    
     5  +coeff( 23)        *x31        
     6  +coeff( 24)            *x41*x51
     7  +coeff( 25)    *x22    *x41    
     8  +coeff( 26)    *x22        *x51
      x_sr_q3en   =x_sr_q3en   
     9  +coeff( 27)    *x21*x31    *x51
     1  +coeff( 28)    *x21*x33        
     2  +coeff( 29)    *x23    *x41    
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)    *x21    *x42*x51
     5  +coeff( 32)    *x24    *x41    
     6  +coeff( 33)    *x21*x33*x41    
     7  +coeff( 34)    *x21*x32*x42    
     8  +coeff( 35)    *x23*x31*x42    
      x_sr_q3en   =x_sr_q3en   
     9  +coeff( 36)    *x21*x31*x44    
     1  +coeff( 37)    *x21*x34*x43    
     2  +coeff( 38)    *x23*x31*x44*x51
c
      return
      end
      function t_sr_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.6737196E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.20929780E-03,-0.79241600E-01, 0.22705630E-01,-0.45621283E-02,
     +  0.41371220E-02,-0.51315520E-02, 0.14413760E-01, 0.88699190E-03,
     +  0.21576392E-02,-0.14410692E-02, 0.94350170E-02,
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
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sr_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)*x11                
      t_sr_q3en   =t_sr_q3en   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21*x31*x41    
c
      return
      end
      function y_sr_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.1449694E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.93785460E-04, 0.89909640E-03,-0.46059620E-01, 0.21850820E-02,
     +  0.19415610E-01,-0.80955080E-02,-0.33588272E-02, 0.98945530E-03,
     + -0.61476431E-01, 0.52051700E-02, 0.50074132E-02,-0.89644873E-02,
     +  0.10138043E-01, 0.70655220E-01, 0.75132040E-02, 0.28026732E-02,
     + -0.12654421E-01,-0.23342594E-02,-0.52386503E-01,-0.96558680E-02,
     +  0.94336430E-04,-0.16370630E-01,-0.98218042E-02,-0.13281410E-02,
     + -0.55601401E-02,-0.53595510E-02,-0.62436354E-02, 0.61651510E-02,
     + -0.30573340E-01, 0.78057893E-02,-0.88579320E-02,-0.19765650E-01,
     +  0.12258974E-01,-0.13500580E-01,-0.11289110E-01,-0.81402110E-02,
     +  0.16547180E-01, 0.14844173E-01,-0.73073930E-02, 0.14039811E-01,
     + -0.25105122E-01,-0.16491271E-02, 0.16199050E-01,-0.15934740E-02,
     +  0.28335010E-01, 0.30499570E-01,-0.61588510E-02,-0.83888390E-01,
     +  0.91638310E-03, 0.96083470E-03,-0.22308000E-03,-0.19313580E-01,
     + -0.88515100E-03,-0.14937960E-02, 0.33657120E-02,-0.44936700E-02,
     +  0.79337780E-03,-0.27286590E-01, 0.68183990E-02, 0.83414213E-02,
     + -0.26132920E-02,-0.10129192E-01, 0.97907613E-03,-0.18343850E-01,
     + -0.60211290E-01, 0.14366170E-02, 0.41569970E-02,-0.12893943E-01,
     +  0.14385112E-04,-0.20055022E-01, 0.25140710E-02,-0.34826763E-02,
     + -0.64094790E-01,-0.54766014E-02,-0.15511980E-01, 0.86635250E-01,
     +  0.22868970E-02,-0.12977532E-01, 0.12654810E+00,-0.82652300E-02,
     + -0.10482062E-01,-0.47948720E-03,-0.22959813E-02,-0.64027630E-02,
     +  0.65778603E-03,-0.30432364E-02, 0.32783010E-02, 0.31584701E-02,
     +  0.58339770E-02, 0.41856760E-02, 0.12777301E-01, 0.12465610E-01,
     + -0.11827900E-02,-0.69136070E-02, 0.20334450E-02,
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
      y_sr_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      y_sr_q3en   =y_sr_q3en   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22*x31        
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 18)        *x33        
     1  +coeff( 19)    *x22    *x41    
     2  +coeff( 20)        *x32*x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)            *x43    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x21*x31    *x51
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)            *x42*x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 27)            *x41*x52
     1  +coeff( 28)    *x23    *x41    
     2  +coeff( 29)    *x22*x31*x41    
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)    *x22    *x42    
     5  +coeff( 32)        *x32*x42    
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)            *x44    
     8  +coeff( 35)    *x22    *x41*x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 36)            *x43*x51
     1  +coeff( 37)    *x23*x31*x41    
     2  +coeff( 38)    *x23    *x42    
     3  +coeff( 39)    *x22*x31*x42    
     4  +coeff( 40)        *x33*x42    
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)    *x22*x31*x41*x51
     7  +coeff( 43)    *x22    *x42*x51
     8  +coeff( 44)    *x24*x32        
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 45)    *x23*x31*x42    
     1  +coeff( 46)    *x23    *x43    
     2  +coeff( 47)    *x24*x31    *x51
     3  +coeff( 48)    *x22    *x44*x51
     4  +coeff( 49)    *x21*x32        
     5  +coeff( 50)*x11        *x41    
     6  +coeff( 51)    *x21*x31*x41    
     7  +coeff( 52)        *x31*x42    
     8  +coeff( 53)        *x31*x41*x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 54)                *x53
     1  +coeff( 55)    *x24            
     2  +coeff( 56)        *x33*x41    
     3  +coeff( 57)*x11        *x42    
     4  +coeff( 58)        *x31*x43    
     5  +coeff( 59)    *x21    *x42*x51
     6  +coeff( 60)        *x32*x43    
     7  +coeff( 61)    *x22*x32    *x51
     8  +coeff( 62)        *x31*x43*x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 63)        *x31    *x54
     1  +coeff( 64)    *x21*x32*x43    
     2  +coeff( 65)    *x22    *x44    
     3  +coeff( 66)        *x32*x44    
     4  +coeff( 67)    *x21*x33*x41*x51
     5  +coeff( 68)        *x31*x44*x51
     6  +coeff( 69)    *x22*x32*x43    
     7  +coeff( 70)    *x22*x31*x44    
     8  +coeff( 71)*x12*x22        *x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 72)*x12*x21    *x41*x51
     1  +coeff( 73)    *x22*x31*x43*x51
     2  +coeff( 74)*x12*x21    *x43    
     3  +coeff( 75)    *x24*x31*x43    
     4  +coeff( 76)    *x22*x33*x43    
     5  +coeff( 77)*x12        *x44    
     6  +coeff( 78)    *x24    *x44    
     7  +coeff( 79)    *x22*x32*x44    
     8  +coeff( 80)*x12*x21    *x42*x51
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 81)    *x22*x34*x42*x51
     1  +coeff( 82)*x11*x21            
     2  +coeff( 83)        *x31    *x52
     3  +coeff( 84)    *x22*x32        
     4  +coeff( 85)*x11    *x31*x41    
     5  +coeff( 86)        *x32*x41*x51
     6  +coeff( 87)    *x23*x32        
     7  +coeff( 88)    *x22*x33        
     8  +coeff( 89)    *x24    *x41    
      y_sr_q3en   =y_sr_q3en   
     9  +coeff( 90)        *x34*x41    
     1  +coeff( 91)    *x21*x31*x43    
     2  +coeff( 92)    *x21    *x44    
     3  +coeff( 93)*x11*x21    *x41*x51
     4  +coeff( 94)        *x32*x42*x51
     5  +coeff( 95)    *x21*x31    *x53
c
      return
      end
      function p_sr_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.1029573E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.30900631E-03,-0.11264130E-01,-0.12536560E-01, 0.18889630E-02,
     + -0.85437064E-03,-0.11294924E-02,-0.12183130E-03,-0.11806994E-01,
     + -0.67604380E-03,-0.15118750E-02, 0.12700321E-02, 0.12946540E-01,
     +  0.16520984E-02, 0.90901240E-03,-0.10064742E-01, 0.55086380E-03,
     + -0.17493340E-02,-0.32788760E-02, 0.41467220E-03,-0.62088161E-03,
     + -0.11612680E-02,-0.13988020E-02,-0.11917350E-02, 0.69990080E-03,
     + -0.35822254E-02, 0.98596200E-03, 0.16280020E-01, 0.13609220E-01,
     +  0.59933972E-03,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sr_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      p_sr_q3en   =p_sr_q3en   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x22        *x51
      p_sr_q3en   =p_sr_q3en   
     9  +coeff( 18)    *x22    *x41*x51
     1  +coeff( 19)    *x21*x33*x41    
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)            *x43    
     5  +coeff( 23)    *x21    *x41*x51
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x22    *x42*x51
     8  +coeff( 26)    *x23    *x43    
      p_sr_q3en   =p_sr_q3en   
     9  +coeff( 27)    *x23*x31*x43    
     1  +coeff( 28)    *x23    *x44    
     2  +coeff( 29)    *x21*x32        
c
      return
      end
      function l_sr_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 20)
      data ncoeff/ 19/
      data avdat/  0.1836018E+02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.49888E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54848E-01, 0.19967E-01, 0.29220E-01, 0.49885E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21838440E-02, 0.32147800E+00, 0.32073870E-01, 0.28518160E-01,
     + -0.32873900E-03, 0.20327840E-01,-0.63310110E-01, 0.83724950E-02,
     + -0.44358642E-02,-0.29988550E-02,-0.50610363E-01, 0.22466250E-02,
     +  0.34943650E-02,-0.85392724E-02,-0.52119302E-02,-0.77154280E-01,
     + -0.65913330E-01,-0.85361780E-02,-0.25100184E-01,
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
      x51 = x5
c
c                  function
c
      l_sr_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x23*x31*x41    
      l_sr_q3en   =l_sr_q3en   
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)            *x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)    *x21    *x43    
      l_sr_q3en   =l_sr_q3en   
     9  +coeff( 18)    *x23*x32*x41    
     1  +coeff( 19)    *x21*x34*x41    
c
      return
      end
      function x_sr_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 29)
      data ncoeff/ 28/
      data avdat/ -0.1112230E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.46635E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54444E-01, 0.19967E-01, 0.29220E-01, 0.49077E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.95429340E-02, 0.72408670E-01, 0.38153870E-03, 0.31645820E+00,
     +  0.35488083E-01,-0.20998971E-01,-0.36397820E-02, 0.14655410E-02,
     + -0.24533020E-01,-0.32348360E-02, 0.15286350E-01,-0.11285392E-02,
     + -0.57281200E-02,-0.14591872E-01,-0.38719470E-02,-0.45049450E-02,
     + -0.17596630E-02,-0.12771530E-02, 0.15431200E-02,-0.10544710E-02,
     +  0.24967773E-02,-0.12099382E-01,-0.32203060E-01,-0.25489762E-01,
     + -0.31488882E-02,-0.29090870E-01,-0.20850520E-01,-0.11660343E-01,
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
      x24 = x23*x2
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
c
c                  function
c
      x_sr_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)    *x24            
     8  +coeff(  8)    *x21    *x41    
      x_sr_q3ex    =x_sr_q3ex    
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x22            
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)        *x31*x41    
      x_sr_q3ex    =x_sr_q3ex    
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)                *x53
     4  +coeff( 22)    *x21*x32*x41    
     5  +coeff( 23)    *x21*x31*x42    
     6  +coeff( 24)    *x21    *x43    
     7  +coeff( 25)    *x21    *x42*x51
     8  +coeff( 26)    *x21*x31*x43    
      x_sr_q3ex    =x_sr_q3ex    
     9  +coeff( 27)    *x21    *x44    
     1  +coeff( 28)    *x21*x34*x42    
c
      return
      end
      function t_sr_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.4451506E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.46635E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54444E-01, 0.19967E-01, 0.29220E-01, 0.49077E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30094722E-02,-0.24234930E-01, 0.39379490E-03, 0.10775150E+00,
     +  0.56202071E-02,-0.97521350E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sr_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
c
      return
      end
      function y_sr_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.3068465E-02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.46635E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54444E-01, 0.19967E-01, 0.29220E-01, 0.49077E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.69063420E-03, 0.10333930E-02,-0.41686080E-01,-0.25507312E-01,
     +  0.12113420E-01,-0.56943000E-02,-0.43119490E-02,-0.18809823E-04,
     + -0.53683303E-01, 0.64017571E-03, 0.99747490E-04,-0.64158720E-02,
     +  0.54643102E-02, 0.52159383E-01, 0.63373770E-02, 0.18448153E-02,
     + -0.63189550E-02,-0.65734220E-03,-0.36228653E-01, 0.16297480E-02,
     + -0.63473230E-02,-0.68743624E-02,-0.37967740E-03,-0.78249853E-02,
     + -0.34793743E-02, 0.52977560E-02,-0.91808550E-02, 0.18108040E-01,
     + -0.15401490E-02, 0.15110281E-01,-0.18306214E-01,-0.65617280E-03,
     +  0.73894471E-02, 0.85552763E-02,-0.12795930E-01,-0.30718792E-01,
     + -0.10095360E-01,-0.55171260E-02, 0.30812384E-01, 0.31207420E-01,
     +  0.40406980E-02,-0.49808580E-01, 0.42360020E-02, 0.26277634E-02,
     + -0.30022021E-02,-0.16778721E-02,-0.61100730E-03,-0.90541522E-02,
     + -0.38728753E-02,-0.52623100E-01,-0.63004190E-01,-0.60393850E-01,
     +  0.45022390E-02, 0.65243830E-02, 0.55177211E-02, 0.28840710E-02,
     + -0.30182240E-01,-0.36746822E-02,-0.13925153E-01,-0.34913800E-02,
     + -0.33789260E-01, 0.43893270E-01,-0.35301250E-01, 0.11684960E-01,
     + -0.50335124E-01,-0.24416801E-03, 0.85943300E-03,-0.46679770E-03,
     + -0.53734844E-02,-0.42245432E-03,-0.10290161E-02, 0.25331211E-03,
     + -0.19901320E-02, 0.46845851E-03,-0.30972480E-02,-0.11492600E-02,
     +  0.82644550E-03, 0.24916520E-02, 0.37417600E-02, 0.51533612E-02,
     + -0.22607810E-02, 0.31132580E-02,-0.11754130E-01, 0.18186680E-02,
     +  0.12123640E-01,-0.19645111E-02,-0.66623762E-02,-0.16556943E-02,
     + -0.88286632E-03, 0.11754783E-01,-0.92584690E-02,-0.83127210E-02,
     + -0.19074860E-02,-0.47894220E-02, 0.70462810E-02,
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
      y_sr_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22*x31        
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 18)        *x33        
     1  +coeff( 19)    *x22    *x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)            *x43    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x21*x31    *x51
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)    *x23    *x41    
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 27)    *x22*x31*x41    
     1  +coeff( 28)    *x21*x32*x41    
     2  +coeff( 29)    *x22    *x42    
     3  +coeff( 30)    *x21    *x43    
     4  +coeff( 31)    *x22    *x41*x51
     5  +coeff( 32)*x11*x22    *x41    
     6  +coeff( 33)    *x23*x31*x41    
     7  +coeff( 34)    *x23    *x42    
     8  +coeff( 35)    *x22*x31*x42    
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 36)    *x22    *x43    
     1  +coeff( 37)    *x22*x31*x41*x51
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)    *x23*x31*x42    
     4  +coeff( 40)    *x23    *x43    
     5  +coeff( 41)    *x23*x34        
     6  +coeff( 42)    *x24    *x44*x51
     7  +coeff( 43)    *x21*x31*x41    
     8  +coeff( 44)    *x24            
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 45)    *x22*x32        
     1  +coeff( 46)    *x22*x31    *x51
     2  +coeff( 47)    *x21    *x42*x51
     3  +coeff( 48)    *x22*x32*x41    
     4  +coeff( 49)    *x22*x32    *x51
     5  +coeff( 50)    *x22*x31*x43    
     6  +coeff( 51)    *x21*x32*x43    
     7  +coeff( 52)    *x22    *x44    
     8  +coeff( 53)    *x21*x33*x41*x51
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 54)    *x23    *x42*x51
     1  +coeff( 55)    *x21    *x42*x53
     2  +coeff( 56)    *x22*x32*x43    
     3  +coeff( 57)    *x22*x31*x44    
     4  +coeff( 58)*x12*x21    *x41*x51
     5  +coeff( 59)    *x22*x31*x43*x51
     6  +coeff( 60)*x12*x21        *x52
     7  +coeff( 61)    *x24*x32*x42    
     8  +coeff( 62)    *x21*x34*x43    
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 63)    *x22*x31*x43*x52
     1  +coeff( 64)*x12*x24    *x41*x51
     2  +coeff( 65)    *x22*x32*x42*x54
     3  +coeff( 66)*x11*x21            
     4  +coeff( 67)*x11        *x41    
     5  +coeff( 68)        *x32*x41    
     6  +coeff( 69)        *x31*x42    
     7  +coeff( 70)            *x41*x52
     8  +coeff( 71)                *x53
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 72)*x11        *x42    
     1  +coeff( 73)        *x31*x43    
     2  +coeff( 74)    *x21*x31    *x52
     3  +coeff( 75)    *x21    *x41*x52
     4  +coeff( 76)            *x41*x53
     5  +coeff( 77)*x12*x21            
     6  +coeff( 78)    *x24    *x41    
     7  +coeff( 79)        *x33*x42    
     8  +coeff( 80)    *x21    *x44    
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 81)        *x31*x43*x51
     1  +coeff( 82)    *x21    *x41*x53
     2  +coeff( 83)    *x21*x34*x41    
     3  +coeff( 84)*x11*x22    *x42    
     4  +coeff( 85)    *x24    *x42    
     5  +coeff( 86)*x11*x21*x31*x42    
     6  +coeff( 87)    *x21*x33*x42    
     7  +coeff( 88)        *x34*x42    
     8  +coeff( 89)*x12*x21        *x51
      y_sr_q3ex    =y_sr_q3ex    
     9  +coeff( 90)    *x24    *x41*x51
     1  +coeff( 91)    *x22*x32*x41*x51
     2  +coeff( 92)    *x22*x31*x42*x51
     3  +coeff( 93)*x11*x21    *x41*x52
     4  +coeff( 94)    *x22*x31    *x53
     5  +coeff( 95)    *x21    *x41*x54
c
      return
      end
      function p_sr_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.2974516E-03/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.46635E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54444E-01, 0.19967E-01, 0.29220E-01, 0.49077E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.35036190E-03,-0.33103200E-03, 0.14955300E-01,-0.46445052E-02,
     + -0.60064680E-02, 0.24376020E-02, 0.41709500E-03, 0.51488482E-03,
     +  0.17519583E-01,-0.75080040E-04, 0.27539850E-02,-0.42016454E-02,
     + -0.21942640E-01,-0.19007534E-02,-0.83464011E-03, 0.17053920E-01,
     +  0.38959292E-03, 0.17704420E-02, 0.22879564E-02, 0.23134970E-02,
     + -0.26260891E-02, 0.12194603E-01, 0.10900420E-01,-0.41460423E-03,
     + -0.80567104E-02,-0.12913210E-01, 0.11728240E-01, 0.15601334E-01,
     +  0.43717500E-02,
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
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
      p_sr_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sr_q3ex    =p_sr_q3ex    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21    *x42    
      p_sr_q3ex    =p_sr_q3ex    
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)    *x24*x31        
     7  +coeff( 25)    *x23*x31*x41    
     8  +coeff( 26)    *x23    *x42    
      p_sr_q3ex    =p_sr_q3ex    
     9  +coeff( 27)    *x22    *x43*x51
     1  +coeff( 28)    *x24    *x44    
     2  +coeff( 29)    *x22*x31        
c
      return
      end
      function l_sr_q3ex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10)!,xmean(x)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.2133930E+02/
      data xmin/
     1 -0.19972E-02,-0.54501E-01,-0.19979E-01,-0.29899E-01,-0.46635E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19884E-02, 0.54444E-01, 0.19967E-01, 0.29220E-01, 0.49077E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.63321064E-03, 0.32105153E+00, 0.30394660E-01, 0.32575763E-01,
     + -0.65469330E-03,-0.62573430E-01, 0.76433440E-02,-0.44968780E-02,
     + -0.30751670E-02, 0.12338270E-01, 0.78210070E-02,-0.49449924E-01,
     +  0.22396490E-02, 0.35996413E-02,-0.84933360E-02,-0.52896010E-02,
     + -0.75872080E-01,-0.65512600E-01,-0.65620220E-02,-0.26024540E-01,
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
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sr_q3ex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21    *x42    
     7  +coeff(  7)    *x23*x31*x41    
     8  +coeff(  8)*x11                
      l_sr_q3ex    =l_sr_q3ex    
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)            *x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x21*x31*x42    
      l_sr_q3ex    =l_sr_q3ex    
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x23*x32*x41    
     2  +coeff( 20)    *x21*x34*x41    
c
      return
      end
      function x_sr_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.9414225E-02/
      data xmin/
     1 -0.19972E-02,-0.54829E-01,-0.99828E-02,-0.29963E-01,-0.45933E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19993E-02, 0.54184E-01, 0.99907E-02, 0.27775E-01, 0.49990E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17082154E-01, 0.31706710E-02, 0.68658840E-04, 0.10155600E-02,
     +  0.62519884E+00, 0.54543150E-01,-0.48334920E-01, 0.22326970E-01,
     + -0.16351813E-01,-0.10864420E-01,-0.45405332E-02,-0.38314824E-02,
     + -0.23479554E-02,-0.26231980E-02, 0.39855030E-02,-0.29205900E-02,
     +  0.52753314E-02,-0.17898970E-02,-0.14897742E-01,-0.26909130E-01,
     + -0.15818921E-02, 0.45236372E-02,-0.98071983E-02,-0.84522934E-02,
     + -0.31411331E-02,-0.30925171E-01,-0.34642620E-03, 0.23394774E-02,
     + -0.54579292E-03,-0.58766650E-03,-0.66701460E-03,-0.93222810E-03,
     + -0.46441360E-02,-0.14731180E-02, 0.33061301E-02,-0.24017423E-01,
     + -0.56848980E-02,-0.11956120E-01,-0.22474360E-05, 0.22076530E-01,
     +  0.52531860E-03, 0.77435730E-02,-0.23693722E-03,-0.35973920E-02,
     +  0.88762930E-03,-0.34134523E-02,-0.62853324E-03, 0.26523340E-02,
     + -0.80291640E-02,-0.19219690E-03,-0.92119520E-02,-0.42994124E-02,
     + -0.23060990E-03, 0.13056790E-01, 0.59540863E-02,-0.48790540E-02,
     +  0.18265140E-02, 0.10298480E-01,-0.74360840E-02, 0.18213400E-02,
     +  0.50777030E-03, 0.54054064E-02,-0.17876170E-02,-0.20231870E-03,
     +  0.71311442E-04, 0.65126013E-03, 0.80317780E-03, 0.11720004E-02,
     + -0.73045010E-03,-0.92823060E-04,-0.13380510E-02, 0.90851910E-03,
     + -0.55457110E-02, 0.28780780E-02,-0.57601680E-03,-0.15075670E-03,
     +  0.15163222E-02, 0.19927604E-02, 0.59074494E-02,-0.98433904E-03,
     + -0.18264450E-02, 0.87062020E-03,-0.10765620E-01,-0.31017400E-03,
     +  0.38332590E-02, 0.30737340E-02, 0.13160043E-02, 0.50572930E-02,
     + -0.93627130E-02, 0.65799980E-03,
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
      x24 = x23*x2
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
      x_sr_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x22            
      x_sr_fp     =x_sr_fp     
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)                *x53
      x_sr_fp     =x_sr_fp     
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x21*x31*x42    
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)        *x31*x41    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x24            
     6  +coeff( 24)    *x21    *x42*x51
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)    *x21    *x44    
      x_sr_fp     =x_sr_fp     
     9  +coeff( 27)        *x31    *x51
     1  +coeff( 28)    *x23            
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)        *x31*x41*x51
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)    *x22*x31*x41    
     6  +coeff( 33)    *x21*x32*x41    
     7  +coeff( 34)    *x21    *x41*x52
     8  +coeff( 35)    *x22    *x43    
      x_sr_fp     =x_sr_fp     
     9  +coeff( 36)    *x21*x31*x43    
     1  +coeff( 37)    *x24        *x51
     2  +coeff( 38)    *x21    *x43*x51
     3  +coeff( 39)    *x21*x31    *x53
     4  +coeff( 40)    *x23    *x43    
     5  +coeff( 41)*x11*x23            
     6  +coeff( 42)    *x23*x31*x42*x51
     7  +coeff( 43)            *x43    
     8  +coeff( 44)    *x23        *x51
      x_sr_fp     =x_sr_fp     
     9  +coeff( 45)    *x22    *x41*x51
     1  +coeff( 46)    *x22        *x52
     2  +coeff( 47)    *x21*x33*x41    
     3  +coeff( 48)    *x22*x31*x42    
     4  +coeff( 49)    *x21*x32*x42    
     5  +coeff( 50)*x11            *x51
     6  +coeff( 51)    *x21*x31*x42*x51
     7  +coeff( 52)    *x23        *x52
     8  +coeff( 53)        *x32*x41*x52
      x_sr_fp     =x_sr_fp     
     9  +coeff( 54)    *x23*x31*x42    
     1  +coeff( 55)    *x22    *x44    
     2  +coeff( 56)    *x21*x31*x44    
     3  +coeff( 57)    *x23*x31*x41*x51
     4  +coeff( 58)    *x23    *x42*x51
     5  +coeff( 59)    *x21*x31*x43*x51
     6  +coeff( 60)    *x24    *x43    
     7  +coeff( 61)    *x23*x32*x41*x51
     8  +coeff( 62)    *x21    *x43*x54
      x_sr_fp     =x_sr_fp     
     9  +coeff( 63)    *x22*x33*x43*x51
     1  +coeff( 64)    *x21*x31        
     2  +coeff( 65)        *x32        
     3  +coeff( 66)            *x41*x52
     4  +coeff( 67)        *x31*x41*x52
     5  +coeff( 68)            *x42*x52
     6  +coeff( 69)                *x54
     7  +coeff( 70)*x11        *x41    
     8  +coeff( 71)    *x23*x31*x41    
      x_sr_fp     =x_sr_fp     
     9  +coeff( 72)    *x22*x32*x41    
     1  +coeff( 73)    *x21*x32*x41*x51
     2  +coeff( 74)    *x22    *x42*x51
     3  +coeff( 75)            *x44*x51
     4  +coeff( 76)*x11*x22            
     5  +coeff( 77)    *x24*x31*x41    
     6  +coeff( 78)    *x23*x32*x41    
     7  +coeff( 79)    *x24    *x42    
     8  +coeff( 80)    *x22*x31*x43    
      x_sr_fp     =x_sr_fp     
     9  +coeff( 81)    *x21*x32*x42*x51
     1  +coeff( 82)    *x22        *x54
     2  +coeff( 83)    *x23    *x44    
     3  +coeff( 84)*x11    *x31*x41*x51
     4  +coeff( 85)    *x21*x34*x41*x51
     5  +coeff( 86)    *x21    *x44*x52
     6  +coeff( 87)    *x23        *x54
     7  +coeff( 88)    *x24*x31*x43    
     8  +coeff( 89)    *x23*x31*x44    
      x_sr_fp     =x_sr_fp     
     9  +coeff( 90)*x11*x21*x31*x41*x51
c
      return
      end
      function t_sr_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 18)
      data ncoeff/ 17/
      data avdat/  0.1201644E-02/
      data xmin/
     1 -0.19972E-02,-0.54829E-01,-0.99828E-02,-0.29963E-01,-0.45933E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19993E-02, 0.54184E-01, 0.99907E-02, 0.27775E-01, 0.49990E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35918470E-02,-0.23732670E-01, 0.45985943E-04, 0.25194630E-03,
     +  0.10735630E+00, 0.59513971E-02,-0.98453790E-02,-0.45717820E-03,
     +  0.76593214E-03,-0.16449372E-02,-0.91194454E-03, 0.56492350E-03,
     + -0.40243414E-03,-0.11544590E-02,-0.71767780E-03, 0.10187060E-02,
     +  0.16581024E-02,
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
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sr_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      t_sr_fp     =t_sr_fp     
     9  +coeff(  9)    *x22            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x52
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)                *x53
     8  +coeff( 17)    *x22    *x42    
c
      return
      end
      function y_sr_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.2276455E-02/
      data xmin/
     1 -0.19972E-02,-0.54829E-01,-0.99828E-02,-0.29963E-01,-0.45933E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19993E-02, 0.54184E-01, 0.99907E-02, 0.27775E-01, 0.49990E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.95784681E-03,-0.33802680E-03,-0.40157270E-01,-0.64089990E-02,
     +  0.97572483E-04,-0.11466610E-02, 0.31979414E-02,-0.35303083E-02,
     + -0.89742941E-02, 0.22768660E-05, 0.18785974E-02, 0.17067680E-02,
     +  0.10742174E-02, 0.26171433E-03, 0.15914350E-02, 0.13919570E-05,
     + -0.15519261E-03, 0.29308714E-05,-0.90009394E-04, 0.10099550E-01,
     + -0.29143970E-02, 0.98826773E-02,-0.54039760E-03,-0.12078762E-03,
     +  0.69012492E-03, 0.68239540E-04, 0.88765122E-03, 0.40843570E-02,
     +  0.86962210E-02,-0.25360864E-05,-0.13152883E-02,-0.59783593E-02,
     + -0.38525070E-02, 0.79826562E-03,-0.13632790E-02, 0.58312254E-03,
     + -0.84705870E-02,-0.31221930E-01,-0.37916231E-05, 0.44491280E-03,
     + -0.14268523E-02, 0.24850584E-01, 0.72541020E-01, 0.30094090E-01,
     +  0.28210854E-01,-0.88411560E-03, 0.73762121E-03, 0.20029700E-02,
     +  0.54223270E-03, 0.16750190E-02, 0.25745140E-03, 0.32980470E-03,
     + -0.68720530E-03, 0.71179580E-04,-0.44603740E-03, 0.20383340E-02,
     +  0.50597481E-01,-0.18246073E-01, 0.10454940E-01,-0.20445514E-01,
     + -0.86400471E-01,-0.27747321E-03, 0.18702753E-03,-0.72165880E-03,
     +  0.11036790E-03,-0.46096794E-03,-0.40666990E-03,-0.42840890E-03,
     + -0.76572000E-03,-0.11639300E-02,-0.18298634E-02, 0.46527520E-02,
     +  0.29046750E-02, 0.92531510E-02,-0.30948991E-04,-0.84364524E-03,
     +  0.11562842E-01,-0.13706500E-01,-0.36157890E-02,-0.18688340E-02,
     + -0.21617960E-01, 0.60362811E-02,-0.54608490E-02,-0.45875290E-02,
     +  0.12458612E-01, 0.52659963E-02,-0.88400654E-02, 0.20802980E-01,
     + -0.31342470E-02, 0.10645600E-01,-0.13245920E-01,-0.11726690E-01,
     +  0.64544421E-02,-0.16593330E-02, 0.67299650E-02,
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
      x34 = x33*x3
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
      y_sr_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)        *x32        
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31    *x51
      y_sr_fp     =y_sr_fp     
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x22            
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x32*x41    
     5  +coeff( 14)        *x31*x42    
     6  +coeff( 15)            *x43    
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)    *x21    *x42    
      y_sr_fp     =y_sr_fp     
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)            *x41*x52
     5  +coeff( 23)    *x23            
     6  +coeff( 24)        *x33*x41    
     7  +coeff( 25)            *x44    
     8  +coeff( 26)    *x21*x33        
      y_sr_fp     =y_sr_fp     
     9  +coeff( 27)    *x21    *x43    
     1  +coeff( 28)    *x22*x31*x41    
     2  +coeff( 29)    *x22    *x42    
     3  +coeff( 30)        *x32    *x52
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)    *x22    *x41*x51
     6  +coeff( 33)            *x41*x53
     7  +coeff( 34)    *x24            
     8  +coeff( 35)            *x45    
      y_sr_fp     =y_sr_fp     
     9  +coeff( 36)    *x22*x33        
     1  +coeff( 37)    *x22*x31*x42    
     2  +coeff( 38)    *x22    *x43    
     3  +coeff( 39)        *x33    *x52
     4  +coeff( 40)        *x31    *x54
     5  +coeff( 41)    *x21        *x54
     6  +coeff( 42)    *x22    *x44    
     7  +coeff( 43)    *x22    *x45    
     8  +coeff( 44)    *x22*x32*x44    
      y_sr_fp     =y_sr_fp     
     9  +coeff( 45)    *x22*x31*x45    
     1  +coeff( 46)    *x21*x31        
     2  +coeff( 47)                *x52
     3  +coeff( 48)    *x22*x31        
     4  +coeff( 49)    *x21*x31    *x51
     5  +coeff( 50)        *x31    *x52
     6  +coeff( 51)    *x21*x31*x42    
     7  +coeff( 52)*x11*x21    *x41    
     8  +coeff( 53)        *x31*x44    
      y_sr_fp     =y_sr_fp     
     9  +coeff( 54)    *x21*x32*x42    
     1  +coeff( 55)            *x43*x52
     2  +coeff( 56)    *x22    *x42*x52
     3  +coeff( 57)    *x22*x31*x44    
     4  +coeff( 58)    *x23    *x45    
     5  +coeff( 59)    *x22*x34*x43    
     6  +coeff( 60)    *x22    *x44*x54
     7  +coeff( 61)*x12*x24*x33*x44*x52
     8  +coeff( 62)    *x21            
      y_sr_fp     =y_sr_fp     
     9  +coeff( 63)*x11*x21            
     1  +coeff( 64)    *x22        *x51
     2  +coeff( 65)                *x53
     3  +coeff( 66)    *x23*x31        
     4  +coeff( 67)    *x22*x31    *x51
     5  +coeff( 68)    *x21    *x41*x52
     6  +coeff( 69)        *x31    *x53
     7  +coeff( 70)        *x34*x41    
     8  +coeff( 71)    *x21    *x41*x53
      y_sr_fp     =y_sr_fp     
     9  +coeff( 72)    *x22*x31*x43    
     1  +coeff( 73)    *x23*x31*x42    
     2  +coeff( 74)    *x22    *x44*x51
     3  +coeff( 75)*x11*x21*x31*x42*x51
     4  +coeff( 76)*x12*x22    *x41    
     5  +coeff( 77)    *x23*x32*x43    
     6  +coeff( 78)    *x23*x31*x44    
     7  +coeff( 79)    *x21    *x45*x52
     8  +coeff( 80)            *x45*x53
      y_sr_fp     =y_sr_fp     
     9  +coeff( 81)    *x24    *x44    
     1  +coeff( 82)    *x23    *x44*x51
     2  +coeff( 83)    *x21    *x44*x53
     3  +coeff( 84)*x11*x22*x31*x42*x51
     4  +coeff( 85)    *x23    *x43*x52
     5  +coeff( 86)*x11*x23*x31*x41*x51
     6  +coeff( 87)    *x22*x34*x42*x51
     7  +coeff( 88)    *x24    *x45    
     8  +coeff( 89)*x11*x22    *x44*x51
      y_sr_fp     =y_sr_fp     
     9  +coeff( 90)    *x22*x32*x42*x53
     1  +coeff( 91)    *x22*x34*x44    
     2  +coeff( 92)    *x23*x34*x43    
     3  +coeff( 93)*x11*x24*x31*x41*x52
     4  +coeff( 94)*x12*x24        *x52
     5  +coeff( 95)*x11*x24*x32*x44    
c
      return
      end
      function p_sr_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1205423E-02/
      data xmin/
     1 -0.19972E-02,-0.54829E-01,-0.99828E-02,-0.29963E-01,-0.45933E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19993E-02, 0.54184E-01, 0.99907E-02, 0.27775E-01, 0.49990E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12374633E-02, 0.71601140E-02,-0.52271530E-02,-0.58831010E-02,
     + -0.10946390E-04, 0.19251704E-02, 0.51379553E-03,-0.81407630E-04,
     +  0.19652660E-01, 0.81200472E-04, 0.19306174E-03, 0.31013851E-02,
     + -0.23663500E-02,-0.21109400E-01,-0.19770840E-02, 0.19247272E-03,
     + -0.11227640E-02, 0.17909130E-02,-0.31851960E-03, 0.10988142E-01,
     +  0.77473954E-03, 0.10193000E-02, 0.21010304E-02, 0.70989710E-03,
     +  0.33705640E-02, 0.10172662E-03,-0.38080913E-02,-0.32292670E-02,
     + -0.52930023E-02,-0.16251640E-02,-0.18599034E-02,-0.80861700E-03,
     + -0.48463040E-04,-0.21423120E-02, 0.24644981E-02, 0.16324280E-01,
     +  0.23135700E-01, 0.52649430E-02,-0.31246590E-02,-0.45259950E-02,
     +  0.33335670E-01,-0.94169620E-03, 0.16139302E-02, 0.76720920E-04,
     +  0.70846593E-03, 0.39422000E-02,-0.19208950E-02,-0.14255141E-02,
     + -0.52252173E-03,-0.42531434E-02, 0.39810100E-02,-0.28948270E-02,
     +  0.10876842E-01, 0.42360710E-02, 0.18071403E-01,-0.26216961E-01,
     +  0.22318041E-03, 0.18900963E-02,-0.31122500E-03,-0.62389241E-03,
     + -0.10317372E-03, 0.53762984E-02, 0.22289822E-02,-0.39512193E-02,
     + -0.17491530E-02, 0.32354280E-03,-0.16680350E-02,-0.13115560E-02,
     +  0.20009460E-02, 0.71540991E-02,-0.91117020E-03,-0.22412512E-02,
     +  0.13932010E-02,-0.24071434E-01,-0.14067370E-01,-0.10235590E-01,
     +  0.13131853E-01, 0.19364714E-01, 0.32903372E-02, 0.73544770E-02,
     + -0.43624532E-02,-0.23195890E-02, 0.31371161E-01,-0.15571710E-02,
     +  0.18590260E-02, 0.85943294E-02, 0.19535442E-01,-0.17395073E-01,
     + -0.95191262E-02,-0.23855680E-01,
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
      p_sr_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sr_fp     =p_sr_fp     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)                *x52
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x23            
      p_sr_fp     =p_sr_fp     
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)            *x43    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)            *x42*x51
     7  +coeff( 25)            *x41*x52
     8  +coeff( 26)*x11*x21    *x41    
      p_sr_fp     =p_sr_fp     
     9  +coeff( 27)    *x23    *x41    
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)            *x44    
     4  +coeff( 31)    *x22    *x41*x51
     5  +coeff( 32)    *x21    *x42*x51
     6  +coeff( 33)            *x43*x51
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x23    *x42    
      p_sr_fp     =p_sr_fp     
     9  +coeff( 36)    *x22*x31*x42    
     1  +coeff( 37)    *x22    *x43    
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)    *x24*x31*x41    
     4  +coeff( 40)    *x23    *x43    
     5  +coeff( 41)    *x22    *x44    
     6  +coeff( 42)    *x21            
     7  +coeff( 43)    *x21    *x41*x51
     8  +coeff( 44)        *x31*x41*x51
      p_sr_fp     =p_sr_fp     
     9  +coeff( 45)        *x31    *x52
     1  +coeff( 46)    *x22*x31*x41    
     2  +coeff( 47)    *x21*x32*x41    
     3  +coeff( 48)    *x21*x31*x42    
     4  +coeff( 49)        *x31*x43    
     5  +coeff( 50)    *x21    *x44    
     6  +coeff( 51)    *x22*x31*x41*x51
     7  +coeff( 52)    *x23*x31*x42    
     8  +coeff( 53)    *x22*x31*x43    
      p_sr_fp     =p_sr_fp     
     9  +coeff( 54)    *x21*x32*x43    
     1  +coeff( 55)    *x22    *x43*x51
     2  +coeff( 56)    *x24    *x43    
     3  +coeff( 57)    *x21*x31*x41    
     4  +coeff( 58)    *x22*x32        
     5  +coeff( 59)        *x33*x41    
     6  +coeff( 60)    *x23        *x51
     7  +coeff( 61)    *x21    *x41*x52
     8  +coeff( 62)    *x24    *x41    
      p_sr_fp     =p_sr_fp     
     9  +coeff( 63)    *x22*x32*x41    
     1  +coeff( 64)    *x21*x31*x43    
     2  +coeff( 65)    *x21    *x43*x51
     3  +coeff( 66)        *x32    *x53
     4  +coeff( 67)    *x24*x32        
     5  +coeff( 68)    *x24    *x42    
     6  +coeff( 69)*x11*x21*x31*x41*x51
     7  +coeff( 70)    *x22*x31*x42*x51
     8  +coeff( 71)    *x22*x31*x41*x52
      p_sr_fp     =p_sr_fp     
     9  +coeff( 72)            *x43*x53
     1  +coeff( 73)    *x23*x33*x41    
     2  +coeff( 74)    *x24*x31*x42    
     3  +coeff( 75)    *x23    *x44    
     4  +coeff( 76)    *x22*x31*x44    
     5  +coeff( 77)    *x24*x31*x43    
     6  +coeff( 78)    *x24    *x44    
     7  +coeff( 79)*x11*x22*x31*x41*x52
     8  +coeff( 80)    *x23    *x43*x52
      p_sr_fp     =p_sr_fp     
     9  +coeff( 81)    *x21*x32*x42*x53
     1  +coeff( 82)*x12*x24    *x41    
     2  +coeff( 83)    *x24*x31*x44    
     3  +coeff( 84)*x12*x22*x31    *x52
     4  +coeff( 85)*x11*x23*x31*x41*x52
     5  +coeff( 86)    *x21*x34*x44*x51
     6  +coeff( 87)    *x22*x31*x43*x54
     7  +coeff( 88)*x11*x23*x31*x43*x52
     8  +coeff( 89)    *x24*x33*x41*x53
      p_sr_fp     =p_sr_fp     
     9  +coeff( 90)    *x23*x32*x44*x53
c
      return
      end
      function l_sr_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 19)
      data ncoeff/ 18/
      data avdat/ -0.1519356E-01/
      data xmin/
     1 -0.19972E-02,-0.54829E-01,-0.99828E-02,-0.29963E-01,-0.45933E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.19993E-02, 0.54184E-01, 0.99907E-02, 0.27775E-01, 0.49990E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15668330E-01,-0.31905782E+00,-0.33070933E-01, 0.44130400E-02,
     + -0.35053590E-01,-0.33357230E-02,-0.24291200E-01, 0.31800480E-01,
     + -0.33933222E-02,-0.90786663E-04, 0.21768214E-01, 0.60267850E-01,
     + -0.65980880E-02,-0.43015424E-02, 0.54756603E-02, 0.38182243E-01,
     +  0.36959763E-01, 0.17033780E-01,
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
      l_sr_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x42    
      l_sr_fp     =l_sr_fp     
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21    *x43    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21        *x51
     6  +coeff( 15)                *x53
     7  +coeff( 16)    *x21*x31*x42    
     8  +coeff( 17)    *x21    *x44    
      l_sr_fp     =l_sr_fp     
     9  +coeff( 18)    *x22    *x44    
c
      return
      end

