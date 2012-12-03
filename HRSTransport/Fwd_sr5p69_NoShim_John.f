c
c Spectrometer transfer functions based on hifield_g2p_dir.dat 
c with target field off.
c                                 7/6/2011 -jjl
c
c
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
c NOMENCLATURE: function name = prefix + _g2_ +suffix
c           prefixes:     x means xfinal
c                         t means thetafinal
c                         y means yfinal
c                         p means phifinal
c                         l means pathlength difference from central trajectory
c     
c           suffixes:     fp means target to focus
c                         sen  means target to septum entrance
c                         sex  means target to septum exit
c                         q1en means target to Q1 entrance
c                         q1ex means target to Q1 exit
c                         dent means target to Dipole entrance (see note below)
c                         dext means target to dipole exit
c                         q3en means target to Q3 entrance
c                         q3ex means target to Q3 exit
c
c
C re the dipole entrance x and theta are in the coordinate system whose origin is at the dipole
c     exit with x pointing outward (+delta gives +x) 30 degrees from perpendicular to the optic 
c     axis (like the exit face of the dipole itself). This is typical of dent functions.
c    see http://hallaweb.jlab.org/news/minutes/g2p/dipole_coords.JPG
      function x_g2_fp(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.8781402E-02/
      data xmin/
     1 -0.49994E-02,-0.47786E-01,-0.29982E-02,-0.26350E-01,-0.42989E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49996E-02, 0.47984E-01, 0.29994E-02, 0.15376E-01, 0.44930E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17705151E-02, 0.96471356E-02, 0.83660870E-03, 0.62333137E+00,
     + -0.26974026E-01, 0.43578647E-01,-0.45374513E-01,-0.15657034E-01,
     +  0.16087025E-01,-0.69484059E-02, 0.54258271E-02,-0.10415395E-01,
     + -0.78732437E-02,-0.26973265E-02,-0.21786441E-02, 0.37640468E-02,
     + -0.15160234E-02, 0.47811139E-02, 0.14042771E-01, 0.14425187E-02,
     + -0.28840003E-02,-0.20466580E-02,-0.79003494E-03, 0.22711938E-02,
     + -0.46452615E-03,-0.15880682E-02,-0.95344346E-03, 0.13008371E-02,
     +  0.79090748E-03,-0.17467159E-03,-0.69730979E-03,-0.45565161E-03,
     +  0.94825192E-03, 0.22708932E-02,-0.16100667E-02, 0.10454981E-02,
     + -0.20447386E-03, 0.13320737E-02,-0.54224571E-02,-0.19723205E-02,
     +  0.84520160E-03,-0.15578985E-02, 0.74052747E-03,-0.80300734E-03,
     + -0.37819514E-03, 0.41773714E-03, 0.13063604E-02, 0.10101206E-02,
     + -0.30601205E-03,-0.28182354E-02, 0.80961245E-03, 0.10008396E-02,
     + -0.13529840E-02, 0.38131056E-03, 0.10531243E-02, 0.75737044E-03,
     + -0.23584313E-04,-0.26763181E-03,-0.55463493E-05,-0.34547604E-04,
     + -0.93555034E-04,-0.15694227E-03, 0.13312448E-02, 0.25254008E-03,
     + -0.56575116E-03, 0.46965102E-03,-0.19057323E-02,-0.10470173E-02,
     + -0.61888096E-03, 0.22799848E-02,-0.31611318E-04,-0.22323368E-03,
     +  0.17962199E-03,-0.26380719E-03, 0.90034649E-03, 0.15102951E-03,
     +  0.16154234E-03,-0.43063847E-03,-0.33089254E-03, 0.41995387E-03,
     + -0.35600053E-03, 0.18424408E-03, 0.23346760E-03, 0.17846485E-03,
     +  0.10037120E-03,-0.18137853E-03, 0.20901265E-03,-0.38471047E-03,
     +  0.16270929E-03, 0.18537132E-03,
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
      x31 = x3
      x32 = x31*x3
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
      x_g2_fp=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      x_g2_fp=x_g2_fp    
     9  +coeff(  9)    *x22            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)    *x21        *x52
      x_g2_fp = x_g2_fp    
     9  +coeff( 18)                *x53
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x24            
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)*x11*x22    *x41    
     7  +coeff( 25)        *x31*x41    
     8  +coeff( 26)    *x21*x31*x41    
      x_g2_fp =x_g2_fp    
     9  +coeff( 27)            *x42*x51
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)        *x31        
     4  +coeff( 31)    *x21*x31    *x51
     5  +coeff( 32)            *x41*x52
     6  +coeff( 33)    *x23*x31        
     7  +coeff( 34)    *x21    *x43    
     8  +coeff( 35)    *x22        *x52
      x_g2_fp =x_g2_fp    
     9  +coeff( 36)            *x42*x52
     1  +coeff( 37)*x11    *x31        
     2  +coeff( 38)    *x24    *x41    
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)    *x24        *x51
     5  +coeff( 41)    *x23    *x41*x51
     6  +coeff( 42)    *x23        *x52
     7  +coeff( 43)*x11*x23            
     8  +coeff( 44)            *x43    
      x_g2_fp =x_g2_fp    
     9  +coeff( 45)        *x31*x41*x51
     1  +coeff( 46)    *x22    *x42    
     2  +coeff( 47)    *x21    *x41*x52
     3  +coeff( 48)            *x41*x53
     4  +coeff( 49)*x11*x21            
     5  +coeff( 50)    *x21    *x44    
     6  +coeff( 51)    *x22    *x42*x51
     7  +coeff( 52)    *x22    *x41*x52
     8  +coeff( 53)            *x42*x53
      x_g2_fp =x_g2_fp    
     9  +coeff( 54)*x11*x21    *x41    
     1  +coeff( 55)    *x23    *x42*x51
     2  +coeff( 56)*x11*x22        *x51
     3  +coeff( 57)    *x23*x31*x42*x51
     4  +coeff( 58)    *x23    *x43*x51
     5  +coeff( 59)*x11        *x43*x51
     6  +coeff( 60)    *x22*x31        
     7  +coeff( 61)    *x21*x32        
     8  +coeff( 62)        *x31*x42    
      x_g2_fp =x_g2_fp    
     9  +coeff( 63)    *x21    *x42*x51
     1  +coeff( 64)        *x31*x41*x52
     2  +coeff( 65)                *x54
     3  +coeff( 66)    *x22    *x43    
     4  +coeff( 67)    *x21    *x43*x51
     5  +coeff( 68)    *x21    *x42*x52
     6  +coeff( 69)            *x41*x54
     7  +coeff( 70)    *x23    *x43    
     8  +coeff( 71)*x11*x21        *x51
      x_g2_fp =x_g2_fp    
     9  +coeff( 72)*x11        *x41*x51
     1  +coeff( 73)*x11            *x52
     2  +coeff( 74)    *x21    *x42*x53
     3  +coeff( 75)            *x42*x54
     4  +coeff( 76)*x11*x22*x31        
     5  +coeff( 77)*x11*x21        *x52
     6  +coeff( 78)*x11*x23    *x41    
     7  +coeff( 79)*x11*x22    *x42    
     8  +coeff( 80)*x11*x23        *x51
      x_g2_fp =x_g2_fp    
     9  +coeff( 81)*x11*x23    *x41*x51
     1  +coeff( 82)    *x21*x31*x42    
     2  +coeff( 83)    *x21*x31*x41*x51
     3  +coeff( 84)            *x43*x51
     4  +coeff( 85)    *x21*x31    *x52
     5  +coeff( 86)    *x21        *x53
     6  +coeff( 87)    *x24*x31        
     7  +coeff( 88)    *x21*x31*x43    
     8  +coeff( 89)    *x23*x31    *x51
      x_g2_fp =x_g2_fp    
     9  +coeff( 90)    *x22*x31*x41*x51
c
      return
      end
      function t_g2_fp(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/  0.1619678E-02/
      data xmin/
     1 -0.49994E-02,-0.47786E-01,-0.29982E-02,-0.26350E-01,-0.42989E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49996E-02, 0.47984E-01, 0.29994E-02, 0.15376E-01, 0.44930E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.43251200E-03,-0.14244510E-01, 0.33313234E-03, 0.11021952E+00,
     +  0.55802925E-02,-0.93246186E-02,-0.17155901E-02, 0.17718530E-02,
     + -0.84743422E-03,-0.11558997E-02, 0.87412406E-03, 0.26888252E-03,
     + -0.14073157E-02,-0.12521374E-02,-0.49343082E-03,-0.28452199E-03,
     +  0.32850876E-04, 0.56113006E-03, 0.49630634E-03, 0.95794688E-03,
     + -0.65872577E-04,-0.88918016E-04,-0.30298973E-03,-0.17993357E-03,
     + -0.47525082E-03, 0.29372890E-03,-0.29359877E-03, 0.36898223E-03,
     + -0.10790992E-03,-0.49504692E-04,-0.14382736E-03,-0.84505373E-04,
     + -0.70782451E-04, 0.11544314E-03, 0.41522009E-04, 0.16539359E-03,
     +  0.11388959E-03,-0.34790320E-03,-0.36470484E-03, 0.32530568E-03,
     + -0.23835782E-03,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_g2_fp =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      t_g2_fp =t_g2_fp    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)                *x53
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)*x11*x21            
      t_g2_fp =t_g2_fp    
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)    *x21*x31        
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)            *x42*x51
     6  +coeff( 24)            *x41*x52
     7  +coeff( 25)    *x24            
     8  +coeff( 26)    *x21    *x43    
      t_g2_fp =t_g2_fp    
     9  +coeff( 27)    *x22        *x52
     1  +coeff( 28)            *x42*x52
     2  +coeff( 29)*x11        *x43    
     3  +coeff( 30)        *x31*x41    
     4  +coeff( 31)    *x21*x31*x41    
     5  +coeff( 32)    *x21*x31    *x51
     6  +coeff( 33)        *x31*x41*x51
     7  +coeff( 34)    *x22    *x42    
     8  +coeff( 35)    *x23        *x51
      t_g2_fp =t_g2_fp    
     9  +coeff( 36)            *x41*x53
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)    *x23    *x42    
     3  +coeff( 39)    *x24        *x51
     4  +coeff( 40)    *x23    *x41*x51
     5  +coeff( 41)    *x23        *x52
c
      return
      end
      function y_g2_fp(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 84)
      data ncoeff/ 83/
      data avdat/ -0.1974921E-02/
      data xmin/
     1 -0.49994E-02,-0.47786E-01,-0.29982E-02,-0.26350E-01,-0.42989E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49996E-02, 0.47984E-01, 0.29994E-02, 0.15376E-01, 0.44930E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.72358805E-03, 0.18402042E-02,-0.13331115E-01, 0.34372861E-03,
     + -0.12571955E-02,-0.22155760E-04,-0.58359141E-03,-0.32994151E-02,
     +  0.91356374E-02,-0.17805202E-02,-0.24258647E-01, 0.10250880E-01,
     +  0.23105214E-02,-0.31235686E-02, 0.97663049E-03,-0.39823100E-03,
     +  0.86291367E-02,-0.10826715E-02, 0.87918928E-02, 0.16762364E-02,
     + -0.29595713E-02, 0.30260612E-02,-0.29136110E-03, 0.17426143E-02,
     + -0.41721523E-03,-0.44723758E-02,-0.26889204E-02,-0.17572701E-03,
     + -0.17331428E-02,-0.17355410E-02,-0.40887095E-04, 0.17295287E-03,
     +  0.58114143E-04,-0.32378040E-04, 0.27454658E-02, 0.13998270E-03,
     +  0.16307367E-02,-0.13009052E-02, 0.79472555E-03,-0.25786604E-02,
     +  0.25359081E-03, 0.29839814E-03, 0.78656980E-04, 0.30215624E-05,
     +  0.12828087E-03, 0.64641068E-03,-0.72535797E-04,-0.11657350E-02,
     + -0.14764497E-02,-0.92695537E-03,-0.21212560E-03, 0.15988910E-02,
     +  0.11690884E-02,-0.55204111E-03,-0.97695843E-03, 0.48052039E-04,
     +  0.56509394E-03, 0.14974162E-02, 0.18350717E-03,-0.64379499E-04,
     +  0.12093570E-03,-0.49681211E-03,-0.92789654E-04, 0.18949340E-03,
     + -0.43094062E-03,-0.13546426E-03,-0.22036691E-03, 0.46475245E-04,
     +  0.65840075E-04,-0.47465716E-03, 0.63800981E-03, 0.48258499E-03,
     +  0.10906568E-03, 0.83203946E-03, 0.12253340E-02,-0.56104193E-03,
     + -0.97243319E-03,-0.29870591E-03, 0.16934698E-02, 0.17420939E-02,
     + -0.42866077E-03, 0.11003189E-02, 0.58665022E-03,
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
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_g2_fp=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_g2_fp=y_g2_fp    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x22            
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)            *x42*x51
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x22    *x41    
      y_g2_fp=y_g2_fp    
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)                *x53
     7  +coeff( 25)    *x21    *x41*x52
     8  +coeff( 26)            *x41*x53
      y_g2_fp=y_g2_fp    
     9  +coeff( 27)    *x24            
     1  +coeff( 28)            *x45    
     2  +coeff( 29)    *x22        *x52
     3  +coeff( 30)                *x54
     4  +coeff( 31)    *x21*x32*x42    
     5  +coeff( 32)    *x21    *x44    
     6  +coeff( 33)    *x22*x33        
     7  +coeff( 34)    *x22*x31    *x52
     8  +coeff( 35)            *x41*x54
      y_g2_fp=y_g2_fp    
     9  +coeff( 36)*x11                
     1  +coeff( 37)            *x43    
     2  +coeff( 38)    *x21    *x42    
     3  +coeff( 39)        *x31    *x52
     4  +coeff( 40)    *x22    *x41*x51
     5  +coeff( 41)*x11*x21        *x51
     6  +coeff( 42)    *x23        *x51
     7  +coeff( 43)    *x21*x31        
     8  +coeff( 44)        *x33        
      y_g2_fp=y_g2_fp    
     9  +coeff( 45)        *x31*x41*x51
     1  +coeff( 46)    *x22*x31        
     2  +coeff( 47)            *x44    
     3  +coeff( 48)            *x43*x51
     4  +coeff( 49)    *x22    *x42    
     5  +coeff( 50)            *x42*x52
     6  +coeff( 51)*x11*x22            
     7  +coeff( 52)    *x22    *x42*x51
     8  +coeff( 53)    *x24    *x41    
      y_g2_fp=y_g2_fp    
     9  +coeff( 54)*x11*x23            
     1  +coeff( 55)    *x21        *x54
     2  +coeff( 56)*x11*x21    *x43    
     3  +coeff( 57)    *x22    *x45    
     4  +coeff( 58)    *x22    *x43*x52
     5  +coeff( 59)        *x31*x42    
     6  +coeff( 60)    *x21*x31*x41    
     7  +coeff( 61)    *x21*x31    *x51
     8  +coeff( 62)    *x21    *x43    
      y_g2_fp=y_g2_fp    
     9  +coeff( 63)*x11        *x42    
     1  +coeff( 64)*x11*x21    *x41    
     2  +coeff( 65)    *x23    *x41    
     3  +coeff( 66)*x11        *x41*x51
     4  +coeff( 67)        *x31    *x53
     5  +coeff( 68)*x12                
     6  +coeff( 69)*x11            *x52
     7  +coeff( 70)    *x21        *x53
     8  +coeff( 71)    *x23    *x42    
      y_g2_fp=y_g2_fp    
     9  +coeff( 72)            *x42*x53
     1  +coeff( 73)*x11*x21    *x41*x51
     2  +coeff( 74)    *x23    *x41*x51
     3  +coeff( 75)    *x21    *x41*x53
     4  +coeff( 76)    *x21    *x44*x51
     5  +coeff( 77)    *x22    *x42*x52
     6  +coeff( 78)*x11*x21    *x41*x52
     7  +coeff( 79)    *x22    *x41*x53
     8  +coeff( 80)    *x21    *x41*x54
      y_g2_fp=y_g2_fp    
     9  +coeff( 81)    *x21*x33*x43*x51
     1  +coeff( 82)    *x22    *x41*x54
     2  +coeff( 83)            *x45*x54
c
      return
      end
      function p_g2_fp(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 58)
      data ncoeff/ 57/
      data avdat/ -0.7563119E-02/
      data xmin/
     1 -0.49994E-02,-0.47786E-01,-0.29982E-02,-0.26350E-01,-0.42989E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49996E-02, 0.47984E-01, 0.29994E-02, 0.15376E-01, 0.44930E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.11147902E-02,-0.13490084E-02, 0.29398846E-02, 0.12559831E-01,
     + -0.38873419E-04, 0.52086799E-02, 0.28451139E-03, 0.11521494E-01,
     + -0.24048614E-02, 0.21224823E-02,-0.20245055E-01,-0.29955595E-02,
     +  0.90853637E-03,-0.27969622E-02,-0.64409617E-03, 0.72941175E-02,
     + -0.16067876E-02, 0.41761398E-02, 0.30253665E-02, 0.99031866E-04,
     +  0.15495290E-04,-0.36206807E-03,-0.10271083E-02, 0.13830076E-02,
     +  0.11358854E-02,-0.23417596E-02,-0.22222379E-02, 0.17930477E-03,
     +  0.38680952E-03, 0.32617079E-03, 0.67821506E-03, 0.36609959E-03,
     +  0.49245061E-03,-0.11006669E-02,-0.80719168E-04, 0.21597515E-03,
     + -0.22284978E-03,-0.50646783E-03,-0.35607102E-03, 0.12714438E-02,
     + -0.16719126E-03, 0.15622507E-03,-0.15303397E-03,-0.33362900E-03,
     +  0.64643020E-04,-0.43258612E-03, 0.28003490E-03,-0.41529522E-03,
     +  0.19350433E-03, 0.79939445E-03,-0.30016556E-03, 0.12838266E-03,
     +  0.79983886E-03,-0.17951477E-03, 0.61944791E-03,-0.27101883E-03,
     + -0.39955892E-03,
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
      x54 = x53*x5
c
c                  function
c
      p_g2_fp=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x41    
      p_g2_fp=p_g2_fp    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x23            
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21    *x42    
      p_g2_fp=p_g2_fp    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)    *x22*x31    *x51
     3  +coeff( 21)        *x33    *x51
     4  +coeff( 22)        *x31*x41    
     5  +coeff( 23)        *x31    *x51
     6  +coeff( 24)            *x43    
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)    *x24            
      p_g2_fp=p_g2_fp    
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)*x11                
     2  +coeff( 29)    *x22*x31        
     3  +coeff( 30)            *x42*x51
     4  +coeff( 31)                *x53
     5  +coeff( 32)*x11*x21        *x51
     6  +coeff( 33)    *x22    *x41*x51
     7  +coeff( 34)            *x41*x53
     8  +coeff( 35)*x11            *x51
      p_g2_fp=p_g2_fp    
     9  +coeff( 36)        *x31    *x52
     1  +coeff( 37)*x11*x22            
     2  +coeff( 38)    *x22        *x52
     3  +coeff( 39)*x11*x23            
     4  +coeff( 40)    *x24    *x41    
     5  +coeff( 41)    *x21*x31*x41    
     6  +coeff( 42)        *x31*x42    
     7  +coeff( 43)*x11        *x42    
     8  +coeff( 44)    *x21    *x43    
      p_g2_fp=p_g2_fp    
     9  +coeff( 45)            *x44    
     1  +coeff( 46)            *x43*x51
     2  +coeff( 47)    *x21    *x41*x52
     3  +coeff( 48)                *x54
     4  +coeff( 49)*x11*x22    *x41    
     5  +coeff( 50)    *x22    *x43    
     6  +coeff( 51)    *x24        *x51
     7  +coeff( 52)*x11*x21    *x41*x51
     8  +coeff( 53)    *x22    *x42*x51
      p_g2_fp=p_g2_fp    
     9  +coeff( 54)    *x23        *x52
     1  +coeff( 55)            *x41*x54
     2  +coeff( 56)            *x44*x52
     3  +coeff( 57)    *x23    *x43*x51
c
      return
      end
      function l_g2_fp(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 19)
      data ncoeff/ 18/
      data avdat/ -0.1028843E-01/
      data xmin/
     1 -0.49994E-02,-0.47786E-01,-0.29982E-02,-0.26350E-01,-0.42989E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49996E-02, 0.47984E-01, 0.29994E-02, 0.15376E-01, 0.44930E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.87742964E-02,-0.22733247E+00,-0.28663604E-01, 0.14186012E-01,
     + -0.22951491E-01, 0.54638159E-01,-0.25658684E-01, 0.33944764E-02,
     +  0.48816809E-02,-0.56713559E-02,-0.78382920E-02,-0.92449579E-02,
     + -0.29603059E-02, 0.25327238E-02, 0.42402595E-02, 0.10877755E-01,
     +  0.44719083E-02,-0.19380245E-01,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_g2_fp=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)                *x52
     8  +coeff(  8)            *x41    
      l_g2_fp=l_g2_fp    
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)    *x23            
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)                *x53
      l_g2_fp=l_g2_fp    
     9  +coeff( 18)    *x23    *x41    
c
      return
      end
      function x_g2_sen(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/  0.1162559E+00/
      data xmin/
     1 -0.49974E-02,-0.48003E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.13544101E-02,-0.30128402E-02,-0.22640238E-01,-0.20694439E-04,
     + -0.82778297E-05, 0.22178841E-04,-0.19722074E-05, 0.29583327E-05,
     +  0.10182307E-04, 0.19846143E-05, 0.20186039E-05,-0.20926261E-05,
     +  0.93245899E-06, 0.16061558E-05,-0.33541057E-05, 0.18572381E-05,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_g2_sen=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21            
     8  +coeff(  8)            *x41*x51
      x_g2_sen=x_g2_sen
       
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x24            
     7  +coeff( 16)*x11*x21    *x41    
c
      return
      end
      function t_g2_sen(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/  0.1118258E+00/
      data xmin/
     1 -0.49974E-02,-0.48003E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12438061E-02,-0.95010277E-04,-0.21648774E-01,-0.20725423E-03,
     +  0.17191263E-03, 0.12398297E-04, 0.32329870E-04, 0.26479796E-04,
     + -0.18789333E-04,-0.31135380E-04, 0.24833646E-04,-0.53361793E-04,
     +  0.30669216E-04, 0.94892648E-05, 0.34126053E-04,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_g2_sen=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22    *x41    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)            *x41*x51
      t_g2_sen=t_g2_sen
       
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x24            
     4  +coeff( 13)*x11*x21    *x41    
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x22    *x42    
c
      return
      end
      function y_g2_sen(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.8106587E-03/
      data xmin/
     1 -0.49974E-02,-0.48003E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.73274411E-03,-0.15657131E-05, 0.51938623E-01, 0.49805613E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x41 = x4
c
c                  function
c
      y_g2_sen=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x11                
c
      return
      end
      function p_g2_sen(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.1024431E-02/
      data xmin/
     1 -0.49974E-02,-0.48003E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.83870214E-03, 0.46970002E-01,-0.53320831E-03,-0.13088199E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x41 = x4
c
c                  function
c
      p_g2_sen=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)*x11                
c
      return
      end
      function l_g2_sen(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.9967133E-03/
      data xmin/
     1 -0.49974E-02,-0.48003E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.58108964E-03, 0.39134625E-05, 0.17321903E-03, 0.14516882E-02,
     + -0.12403060E-02,-0.23651888E-03,
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
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      l_g2_sen=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
c
      return
      end
      function x_g2_sex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 89)
      data ncoeff/ 88/
      data avdat/  0.2755160E+00/
      data xmin/
     1 -0.49974E-02,-0.48003E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.36580428E-02,-0.45147222E-01,-0.25843868E-02, 0.34727538E-02,
     +  0.16340041E-03,-0.19707766E-05, 0.34838824E-05, 0.52716519E-03,
     + -0.32843149E-02, 0.14147736E-02, 0.52933989E-03, 0.23609080E-03,
     + -0.55857608E-03,-0.14325891E-03, 0.11703598E-03, 0.12073132E-03,
     + -0.37466941E-03, 0.56183897E-03, 0.23223880E-03,-0.11851197E-03,
     + -0.14552867E-03, 0.45330431E-04, 0.26554682E-04,-0.21933154E-04,
     +  0.87559878E-04,-0.78536112E-04, 0.12386676E-03,-0.23715077E-03,
     + -0.11966752E-03,-0.77746583E-04,-0.41537882E-04, 0.21034813E-04,
     + -0.46458397E-04, 0.22533301E-04,-0.25775709E-04, 0.81882070E-04,
     + -0.68941146E-04, 0.11610660E-04, 0.11653516E-04,-0.26387886E-04,
     +  0.23564256E-04,-0.47771398E-04, 0.74617565E-04,-0.16596330E-03,
     + -0.48077177E-05,-0.59790250E-05, 0.11924026E-04,-0.77968543E-05,
     + -0.30319035E-04, 0.31004522E-04, 0.67035842E-04,-0.41574189E-04,
     +  0.19235340E-04,-0.85660085E-05, 0.84570647E-05,-0.35443918E-05,
     + -0.76942752E-05,-0.58129581E-06,-0.29504417E-05,-0.15430725E-04,
     +  0.40612267E-05, 0.28951729E-05, 0.26813802E-05, 0.12703897E-04,
     + -0.15227136E-04, 0.13949886E-04,-0.89061004E-05,-0.12566124E-04,
     + -0.34096745E-04, 0.37217740E-05, 0.22474762E-04,-0.40573918E-05,
     + -0.12203740E-04,-0.17257713E-04, 0.18876335E-04, 0.88584566E-05,
     +  0.70558381E-05, 0.93929229E-05, 0.22069551E-04, 0.33935165E-04,
     + -0.30376174E-04,-0.13571133E-04,-0.92226410E-05, 0.11050524E-04,
     +  0.10148733E-04,-0.20679505E-04, 0.20923115E-04,-0.16251950E-04,
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
      x_g2_sex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x22*x31        
     6  +coeff(  6)        *x33        
     7  +coeff(  7)    *x22*x33        
     8  +coeff(  8)    *x21            
      x_g2_sex=x_g2_sex
       
     9  +coeff(  9)        *x31        
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)                *x52
     8  +coeff( 17)    *x24            
      x_g2_sex=x_g2_sex
       
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)    *x23            
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)*x11                
     5  +coeff( 23)    *x21*x31        
     6  +coeff( 24)    *x21        *x51
     7  +coeff( 25)    *x21    *x42    
     8  +coeff( 26)    *x23    *x41    
      x_g2_sex=x_g2_sex
       
     9  +coeff( 27)    *x22*x31*x41    
     1  +coeff( 28)    *x24    *x41    
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)            *x43    
     4  +coeff( 31)    *x22    *x41*x51
     5  +coeff( 32)*x11        *x41    
     6  +coeff( 33)*x11*x22            
     7  +coeff( 34)*x11*x21*x31        
     8  +coeff( 35)*x11*x21        *x51
      x_g2_sex=x_g2_sex
       
     9  +coeff( 36)*x11*x21    *x42    
     1  +coeff( 37)*x11*x23    *x41    
     2  +coeff( 38)        *x31    *x51
     3  +coeff( 39)    *x21*x31*x41    
     4  +coeff( 40)        *x31*x42    
     5  +coeff( 41)            *x42*x51
     6  +coeff( 42)    *x23    *x42    
     7  +coeff( 43)    *x22    *x43    
     8  +coeff( 44)    *x24    *x42    
      x_g2_sex=x_g2_sex
       
     9  +coeff( 45)        *x32        
     1  +coeff( 46)    *x21    *x41*x51
     2  +coeff( 47)        *x31*x41*x51
     3  +coeff( 48)                *x53
     4  +coeff( 49)    *x24*x31        
     5  +coeff( 50)    *x22*x31*x42    
     6  +coeff( 51)    *x22    *x44    
     7  +coeff( 52)*x11*x22    *x41    
     8  +coeff( 53)*x11*x21*x31*x41    
      x_g2_sex=x_g2_sex
       
     9  +coeff( 54)*x11*x21    *x41*x51
     1  +coeff( 55)        *x31    *x52
     2  +coeff( 56)            *x41*x52
     3  +coeff( 57)    *x23*x31        
     4  +coeff( 58)    *x22*x32        
     5  +coeff( 59)        *x31*x43    
     6  +coeff( 60)            *x44    
     7  +coeff( 61)    *x23        *x51
     8  +coeff( 62)    *x22        *x52
      x_g2_sex=x_g2_sex
       
     9  +coeff( 63)*x11    *x31        
     1  +coeff( 64)    *x24        *x51
     2  +coeff( 65)    *x22    *x42*x51
     3  +coeff( 66)    *x21*x32    *x52
     4  +coeff( 67)        *x31*x41*x53
     5  +coeff( 68)        *x31    *x54
     6  +coeff( 69)    *x24*x31*x41    
     7  +coeff( 70)*x11        *x42    
     8  +coeff( 71)    *x22*x31*x43    
      x_g2_sex=x_g2_sex
       
     9  +coeff( 72)    *x24*x31    *x51
     1  +coeff( 73)        *x34    *x52
     2  +coeff( 74)    *x23*x32    *x52
     3  +coeff( 75)*x11*x24            
     4  +coeff( 76)*x11    *x33*x41    
     5  +coeff( 77)    *x21*x34*x43    
     6  +coeff( 78)*x11*x23        *x51
     7  +coeff( 79)    *x22*x34    *x52
     8  +coeff( 80)*x11*x24    *x41    
      x_g2_sex=x_g2_sex
       
     9  +coeff( 81)*x11*x23    *x42    
     1  +coeff( 82)*x11    *x33*x41*x51
     2  +coeff( 83)*x11*x24*x31*x41    
     3  +coeff( 84)*x11*x21*x32*x43    
     4  +coeff( 85)*x11    *x32*x41*x53
     5  +coeff( 86)*x11*x23*x33*x41    
     6  +coeff( 87)*x11*x22*x31*x43*x51
     7  +coeff( 88)*x11    *x31*x43*x54
c
      return
      end
      function t_g2_sex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 46)
      data ncoeff/ 45/
      data avdat/  0.2258550E+00/
      data xmin/
     1 -0.49974E-02,-0.48003E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.38417920E-02, 0.10530937E-02,-0.36362905E-03,-0.25177527E-01,
     + -0.53919186E-02, 0.80045722E-02, 0.10659258E-02,-0.11664447E-02,
     +  0.17798939E-02, 0.29831726E-03, 0.25663761E-03,-0.78491349E-03,
     +  0.12328482E-02, 0.77617682E-04,-0.25522962E-03, 0.16442355E-03,
     + -0.21586464E-03,-0.34008085E-03, 0.29307310E-03,-0.17135852E-04,
     + -0.41821510E-04, 0.16310650E-03, 0.19424431E-03,-0.22839155E-03,
     +  0.29030884E-04, 0.20470744E-03,-0.93915376E-04, 0.39296472E-04,
     + -0.45549154E-04, 0.13194578E-03, 0.18485405E-04, 0.24348115E-04,
     + -0.28705274E-04, 0.36918307E-04, 0.29978404E-04,-0.74377487E-04,
     + -0.42830405E-04,-0.18543139E-03,-0.24132618E-03,-0.10342539E-04,
     +  0.11753244E-04, 0.24986619E-04,-0.70737464E-04,-0.24678693E-04,
     + -0.56692999E-04,
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
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_g2_sex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)            *x42    
      t_g2_sex=t_g2_sex
       
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x24            
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)*x11                
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)    *x23            
      t_g2_sex=t_g2_sex
       
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)    *x24*x31        
     3  +coeff( 21)    *x21        *x51
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)    *x22*x31*x41    
     6  +coeff( 24)*x11*x23            
     7  +coeff( 25)    *x21*x31        
     8  +coeff( 26)    *x22*x31        
      t_g2_sex=t_g2_sex
       
     9  +coeff( 27)            *x43    
     1  +coeff( 28)*x12                
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)*x11*x21    *x42    
     4  +coeff( 31)*x11        *x41    
     5  +coeff( 32)    *x21*x31*x41    
     6  +coeff( 33)        *x31*x42    
     7  +coeff( 34)            *x42*x51
     8  +coeff( 35)*x11*x21*x31        
      t_g2_sex=t_g2_sex
       
     9  +coeff( 36)    *x23    *x41    
     1  +coeff( 37)*x11*x21        *x51
     2  +coeff( 38)    *x24    *x41    
     3  +coeff( 39)    *x24    *x42    
     4  +coeff( 40)        *x32        
     5  +coeff( 41)        *x31    *x51
     6  +coeff( 42)*x11*x21*x31*x41    
     7  +coeff( 43)    *x23    *x42    
     8  +coeff( 44)*x12*x22            
      t_g2_sex=t_g2_sex
       
     9  +coeff( 45)*x11*x23    *x41    
c
      return
      end
      function y_g2_sex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/ -0.1959047E-02/
      data xmin/
     1 -0.49974E-02,-0.48003E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14073198E-02,-0.18174591E-03, 0.93360923E-01,-0.32692046E-02,
     +  0.46423152E-02, 0.10650293E-02,-0.28981327E-03, 0.25829003E-03,
     + -0.42894177E-03,-0.17919914E-03, 0.72151900E-03,
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
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
c
c                  function
c
      y_g2_sex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x22            
      y_g2_sex=y_g2_sex
       
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x23    *x41    
c
      return
      end
      function p_g2_sex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/ -0.1328049E-02/
      data xmin/
     1 -0.49974E-02,-0.48003E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.66317111E-03, 0.41487791E-01,-0.58794320E-02,-0.61094074E-03,
     + -0.38631895E-03,-0.46972060E-03, 0.13995623E-02, 0.15935197E-02,
     +  0.33713697E-03, 0.19830985E-03,-0.39617694E-03, 0.29621512E-03,
     + -0.59586792E-03, 0.21552861E-03, 0.32333893E-03, 0.30883341E-03,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      p_g2_sex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)            *x41    
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x23    *x41    
      p_g2_sex=p_g2_sex
       
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)*x11        *x41    
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)*x11*x22    *x41    
c
      return
      end
      function l_g2_sex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.2513543E-02/
      data xmin/
     1 -0.49974E-02,-0.48003E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12743601E-02, 0.20611724E-03, 0.42190482E-02, 0.37320625E-03,
     + -0.26369756E-02,-0.42449869E-03,-0.49283106E-04,-0.74925818E-04,
     + -0.53906813E-04,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_g2_sex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21            
     8  +coeff(  8)            *x41*x51
      l_g2_sex=l_g2_sex
       
     9  +coeff(  9)*x11*x21            
c
      return
      end
      function x_g2_q1en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 55)
      data ncoeff/ 54/
      data avdat/ -0.1838974E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11621967E-02, 0.10204443E+00,-0.41760341E-02, 0.44378461E-02,
     + -0.34394386E-03, 0.12723680E-02, 0.35951138E-03,-0.43616106E-03,
     +  0.24507765E-03,-0.71970874E-03, 0.11367291E-02,-0.34924186E-03,
     +  0.33166338E-03, 0.37703052E-03,-0.37394777E-04,-0.13852598E-03,
     +  0.18838860E-03, 0.33123529E-03, 0.21203625E-04,-0.61037936E-04,
     +  0.11122810E-03,-0.14571920E-03, 0.69854133E-04,-0.36857131E-04,
     +  0.28722961E-04,-0.61000428E-04,-0.10289094E-03, 0.15593336E-03,
     +  0.18145196E-04, 0.56392186E-04,-0.51794941E-04,-0.11808797E-04,
     +  0.10644409E-04, 0.17386667E-04, 0.96318465E-04,-0.41100957E-04,
     + -0.37543559E-04,-0.33912147E-04, 0.28050834E-04,-0.83286715E-04,
     + -0.42787848E-04,-0.11242628E-03,-0.87511453E-05,-0.92852952E-05,
     + -0.99024737E-05,-0.76119977E-05, 0.15684791E-04,-0.46027373E-04,
     + -0.77508139E-05,-0.66024368E-04, 0.16605321E-04, 0.10773731E-04,
     +  0.28423759E-04, 0.22636670E-04,
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
      x_g2_q1en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)            *x41    
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      x_g2_q1en=x_g2_q1en
       
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)*x11        *x41    
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)        *x31        
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x21    *x41*x51
      x_g2_q1en=x_g2_q1en
       
     9  +coeff( 18)*x11*x22    *x41    
     1  +coeff( 19)                *x51
     2  +coeff( 20)            *x42    
     3  +coeff( 21)    *x23*x31        
     4  +coeff( 22)    *x21    *x43    
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)*x11    *x31        
     7  +coeff( 25)    *x22*x31        
     8  +coeff( 26)    *x24            
      x_g2_q1en=x_g2_q1en
       
     9  +coeff( 27)    *x24    *x41    
     1  +coeff( 28)    *x23    *x42    
     2  +coeff( 29)*x11            *x51
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)*x11        *x42    
     5  +coeff( 32)        *x31*x41    
     6  +coeff( 33)            *x41*x51
     7  +coeff( 34)    *x21*x31    *x51
     8  +coeff( 35)    *x22    *x42    
      x_g2_q1en=x_g2_q1en
       
     9  +coeff( 36)    *x21*x31*x42    
     1  +coeff( 37)    *x23        *x51
     2  +coeff( 38)*x11*x23            
     3  +coeff( 39)*x11*x22*x31        
     4  +coeff( 40)*x11*x24            
     5  +coeff( 41)*x11*x23    *x41    
     6  +coeff( 42)*x11*x24    *x41    
     7  +coeff( 43)    *x21*x32        
     8  +coeff( 44)            *x43    
      x_g2_q1en=x_g2_q1en
       
     9  +coeff( 45)    *x22        *x51
     1  +coeff( 46)    *x21        *x52
     2  +coeff( 47)    *x21    *x42*x51
     3  +coeff( 48)    *x23    *x41*x51
     4  +coeff( 49)*x11    *x31*x41    
     5  +coeff( 50)    *x24    *x42    
     6  +coeff( 51)    *x22*x31*x43    
     7  +coeff( 52)*x11        *x41*x51
     8  +coeff( 53)*x11*x22    *x42    
      x_g2_q1en=x_g2_q1en
       
     9  +coeff( 54)*x11*x23*x31*x42    
c
      return
      end
      function t_g2_q1en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 36)
      data ncoeff/ 35/
      data avdat/ -0.1152366E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.49797347E-03, 0.40216405E-01,-0.59127626E-02,-0.59622794E-03,
     + -0.42652580E-03,-0.58223028E-03, 0.12482365E-02,-0.41802315E-03,
     +  0.16722473E-02, 0.38178690E-03, 0.23669546E-03, 0.42501179E-03,
     + -0.58224227E-03, 0.27740389E-03, 0.30842228E-03, 0.31773924E-03,
     + -0.43272070E-04, 0.47098249E-04,-0.12959592E-03, 0.15101016E-03,
     + -0.25753857E-03, 0.18593588E-04,-0.37821432E-04,-0.31893116E-04,
     +  0.36948812E-04, 0.14047937E-04, 0.37723519E-04, 0.18084873E-04,
     +  0.24794303E-04,-0.84361534E-04,-0.39797105E-04,-0.61178202E-04,
     +  0.23860033E-04,-0.12307892E-03,-0.77300923E-04,
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      t_g2_q1en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)            *x41    
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x23            
     8  +coeff(  8)*x11        *x41    
      t_g2_q1en=t_g2_q1en
       
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)*x11*x22    *x41    
     8  +coeff( 17)        *x31        
      t_g2_q1en=t_g2_q1en
       
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)    *x21    *x43    
     4  +coeff( 22)                *x51
     5  +coeff( 23)            *x42    
     6  +coeff( 24)*x11    *x31        
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)            *x41*x51
      t_g2_q1en=t_g2_q1en
       
     9  +coeff( 27)    *x22*x31        
     1  +coeff( 28)*x11            *x51
     2  +coeff( 29)    *x21*x31    *x51
     3  +coeff( 30)    *x24            
     4  +coeff( 31)*x11        *x42    
     5  +coeff( 32)    *x21*x31*x42    
     6  +coeff( 33)*x12*x21            
     7  +coeff( 34)    *x24    *x41    
     8  +coeff( 35)    *x23    *x41*x51
c
      return
      end
      function y_g2_q1en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/ -0.2055924E-01/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.40562004E-02, 0.32501528E-02, 0.49601253E-01,-0.65864390E-03,
     +  0.36834758E-02,-0.50372188E-02,-0.15979479E-02, 0.97040931E-03,
     + -0.65220176E-03,-0.23685298E-03, 0.14435267E-03,-0.16917138E-03,
     +  0.17394101E-03,-0.61342446E-03,-0.23133087E-03, 0.57818473E-03,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_g2_q1en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)            *x42    
      y_g2_q1en=y_g2_q1en
       
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)*x11*x21    *x41    
     7  +coeff( 16)    *x24            
c
      return
      end
      function p_g2_q1en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/ -0.8276522E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.32518557E-02,-0.95952157E-03, 0.33831652E-03, 0.23876345E-01,
     +  0.53349445E-02,-0.72197574E-02, 0.13311724E-02,-0.10170753E-02,
     + -0.15414566E-02,-0.24394463E-03,-0.24381671E-03, 0.92733291E-03,
     + -0.10482009E-02,-0.54844480E-04, 0.25411620E-03, 0.23709841E-03,
     +  0.28107458E-03,-0.23766897E-03,-0.76329598E-04,-0.17574590E-03,
     + -0.12969444E-03,-0.20077068E-03, 0.25824201E-03,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_g2_q1en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)*x11*x21            
      p_g2_q1en=p_g2_q1en
       
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x24            
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)*x11                
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22        *x51
      p_g2_q1en=p_g2_q1en
       
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)*x11*x23            
c
      return
      end
      function l_g2_q1en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 12)
      data ncoeff/ 11/
      data avdat/  0.4485654E-03/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.27010E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49991E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44960E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.71436161E-03,-0.36033423E-03,-0.36920626E-02,-0.82758976E-04,
     + -0.22462031E-02,-0.54461235E-03, 0.29359353E-03, 0.37590064E-04,
     + -0.85280015E-04, 0.29930354E-04, 0.32389373E-04,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_g2_q1en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x21            
      l_g2_q1en=l_g2_q1en
       
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)*x11*x21            
c
      return
      end
      function x_g2_q1ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2029334E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.68230735E-03, 0.11677490E+00,-0.10997940E-01, 0.28860855E-02,
     +  0.15895150E-02,-0.81127329E-03,-0.11159887E-02, 0.23060956E-02,
     +  0.73414005E-03, 0.34040872E-02,-0.82422048E-03,-0.16715551E-03,
     +  0.82036655E-03,-0.14151938E-02, 0.73736260E-03,-0.83352490E-04,
     + -0.28520034E-03, 0.12770707E-03, 0.90469036E-03,-0.11299409E-03,
     +  0.27819851E-03,-0.12103504E-03, 0.32066228E-03, 0.13759795E-03,
     + -0.84468018E-04, 0.21477592E-04, 0.67354813E-04,-0.12957356E-03,
     + -0.36560706E-03,-0.23097872E-03, 0.14213266E-03,-0.11226054E-03,
     + -0.29576246E-04, 0.21620855E-04,-0.84918720E-04,-0.62659448E-04,
     +  0.61954204E-04,-0.16084427E-03,-0.11008210E-03,-0.83830230E-06,
     +  0.28012162E-04, 0.10419914E-03,-0.22020057E-03,-0.15916831E-03,
     + -0.66678206E-06,-0.19969204E-04, 0.11752072E-04, 0.12174759E-03,
     + -0.37107100E-04,-0.30472971E-03, 0.25334662E-05, 0.23275621E-04,
     + -0.33158005E-04,-0.58671858E-05, 0.40102899E-04, 0.32703716E-04,
     + -0.18929401E-05,-0.21347307E-04,-0.26698186E-04,-0.76493983E-04,
     + -0.16812330E-04,-0.13849896E-04, 0.21319721E-03,-0.19133176E-04,
     +  0.23519421E-03, 0.28715353E-03,-0.27404901E-05,-0.63663880E-04,
     + -0.13504236E-03,-0.29890917E-03, 0.49999075E-04,-0.10924682E-03,
     +  0.18627272E-03,-0.20102831E-03, 0.74163079E-04,-0.31301845E-03,
     + -0.18527353E-04, 0.14312319E-03, 0.29287400E-03, 0.80437647E-04,
     +  0.27790013E-04, 0.64523476E-04, 0.56827143E-04,-0.36553072E-04,
     + -0.48443413E-04, 0.31566557E-04, 0.62458865E-04,-0.30838535E-03,
     +  0.71389397E-04,-0.13708763E-03,
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
      x_g2_q1ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x23            
      x_g2_q1ex=x_g2_q1ex
       
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)*x11        *x41    
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)        *x31        
     8  +coeff( 17)    *x21*x31*x41    
      x_g2_q1ex=x_g2_q1ex
       
     9  +coeff( 18)*x11            *x51
     1  +coeff( 19)*x11*x22    *x41    
     2  +coeff( 20)            *x42    
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)    *x23*x31        
     6  +coeff( 24)*x11*x21            
     7  +coeff( 25)*x11    *x31        
     8  +coeff( 26)                *x51
      x_g2_q1ex=x_g2_q1ex
       
     9  +coeff( 27)    *x22*x31        
     1  +coeff( 28)    *x24            
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)    *x24    *x41    
     4  +coeff( 31)*x11*x21    *x41    
     5  +coeff( 32)*x11        *x42    
     6  +coeff( 33)        *x31*x41    
     7  +coeff( 34)            *x41*x51
     8  +coeff( 35)    *x21*x31*x42    
      x_g2_q1ex=x_g2_q1ex
       
     9  +coeff( 36)*x11*x23            
     1  +coeff( 37)*x11*x22*x31        
     2  +coeff( 38)*x11*x24            
     3  +coeff( 39)*x11*x23    *x41    
     4  +coeff( 40)    *x21*x32        
     5  +coeff( 41)    *x21*x31    *x51
     6  +coeff( 42)    *x22    *x42    
     7  +coeff( 43)    *x21    *x44    
     8  +coeff( 44)    *x23    *x41*x51
      x_g2_q1ex=x_g2_q1ex
       
     9  +coeff( 45)*x11*x21*x31        
     1  +coeff( 46)*x11    *x31*x41    
     2  +coeff( 47)*x11        *x41*x51
     3  +coeff( 48)    *x23    *x42*x51
     4  +coeff( 49)    *x22*x32*x43    
     5  +coeff( 50)*x11*x24    *x41    
     6  +coeff( 51)        *x31    *x51
     7  +coeff( 52)        *x32*x41    
     8  +coeff( 53)            *x43    
      x_g2_q1ex=x_g2_q1ex
       
     9  +coeff( 54)    *x22*x31*x41    
     1  +coeff( 55)    *x23        *x51
     2  +coeff( 56)    *x21    *x42*x51
     3  +coeff( 57)    *x22*x33        
     4  +coeff( 58)    *x21*x34        
     5  +coeff( 59)        *x34*x41    
     6  +coeff( 60)    *x21*x31*x43    
     7  +coeff( 61)        *x31*x44    
     8  +coeff( 62)    *x22    *x42*x51
      x_g2_q1ex=x_g2_q1ex
       
     9  +coeff( 63)    *x22*x33*x41    
     1  +coeff( 64)    *x21*x34*x41    
     2  +coeff( 65)    *x23    *x43    
     3  +coeff( 66)    *x21*x32*x42*x51
     4  +coeff( 67)*x11            *x52
     5  +coeff( 68)*x11    *x32    *x51
     6  +coeff( 69)    *x23    *x43*x51
     7  +coeff( 70)    *x24*x33*x41    
     8  +coeff( 71)*x11*x21*x31*x42    
      x_g2_q1ex=x_g2_q1ex
       
     9  +coeff( 72)    *x24*x32*x42    
     1  +coeff( 73)    *x24*x31*x43    
     2  +coeff( 74)    *x22*x33*x43    
     3  +coeff( 75)    *x22*x33*x42*x51
     4  +coeff( 76)    *x21*x34*x42*x51
     5  +coeff( 77)*x11*x22        *x52
     6  +coeff( 78)    *x24*x31*x41*x52
     7  +coeff( 79)    *x23*x32*x41*x52
     8  +coeff( 80)    *x22*x32*x42*x52
      x_g2_q1ex=x_g2_q1ex
       
     9  +coeff( 81)    *x21*x33*x41*x53
     1  +coeff( 82)*x11*x22*x31*x42    
     2  +coeff( 83)*x11    *x34    *x51
     3  +coeff( 84)*x11*x21*x32*x41*x51
     4  +coeff( 85)*x11*x22    *x41*x52
     5  +coeff( 86)*x11*x24*x32        
     6  +coeff( 87)*x11*x22    *x44    
     7  +coeff( 88)    *x23*x34*x41*x52
     8  +coeff( 89)    *x24    *x44*x52
      x_g2_q1ex=x_g2_q1ex
       
     9  +coeff( 90)    *x21*x32*x44*x53
c
      return
      end
      function t_g2_q1ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 53)
      data ncoeff/ 52/
      data avdat/ -0.1874139E-03/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12247649E-04,-0.11758026E-01,-0.22146234E-02,-0.25954398E-02,
     +  0.22621911E-02,-0.24441685E-03,-0.16789460E-03, 0.12020129E-03,
     +  0.15618817E-03,-0.16079204E-03,-0.31538017E-03, 0.12897578E-02,
     +  0.18247291E-03, 0.73645489E-04,-0.69223315E-04,-0.10524805E-03,
     +  0.99991485E-04,-0.14600964E-04,-0.70397204E-04, 0.19069585E-03,
     + -0.63454121E-03,-0.13379597E-04, 0.20858732E-04,-0.14889048E-04,
     +  0.86748492E-04, 0.12880584E-03,-0.11542302E-04, 0.12688374E-04,
     + -0.10089369E-04,-0.34310571E-04, 0.18683162E-04,-0.22966340E-04,
     +  0.13965825E-03, 0.11459381E-03, 0.14435655E-04,-0.21920261E-03,
     + -0.18924696E-03,-0.14164168E-03, 0.19801639E-03,-0.31427442E-05,
     + -0.50824065E-05,-0.91971415E-05, 0.24798530E-04, 0.76995275E-05,
     + -0.12876155E-04, 0.12119883E-04,-0.26359628E-04,-0.32425807E-04,
     +  0.16932938E-04,-0.34022327E-04,-0.33882592E-04, 0.11207156E-03,
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_g2_q1ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)            *x41    
     8  +coeff(  8)    *x22            
      t_g2_q1ex=t_g2_q1ex
       
     9  +coeff(  9)    *x23            
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)*x11*x22            
      t_g2_q1ex=t_g2_q1ex
       
     9  +coeff( 18)        *x31        
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)*x11*x22    *x41    
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)                *x51
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)*x11    *x31        
     7  +coeff( 25)    *x23*x31        
     8  +coeff( 26)    *x23        *x51
      t_g2_q1ex=t_g2_q1ex
       
     9  +coeff( 27)            *x42    
     1  +coeff( 28)    *x22        *x51
     2  +coeff( 29)    *x21*x31    *x51
     3  +coeff( 30)    *x24            
     4  +coeff( 31)*x11*x21    *x41    
     5  +coeff( 32)*x11        *x42    
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)    *x21    *x42*x51
     8  +coeff( 35)    *x24*x31        
      t_g2_q1ex=t_g2_q1ex
       
     9  +coeff( 36)    *x21    *x44    
     1  +coeff( 37)    *x23    *x41*x51
     2  +coeff( 38)    *x21    *x43*x51
     3  +coeff( 39)    *x23    *x44    
     4  +coeff( 40)            *x41*x51
     5  +coeff( 41)    *x22    *x42    
     6  +coeff( 42)*x11        *x41*x51
     7  +coeff( 43)    *x21    *x41*x52
     8  +coeff( 44)*x12*x21            
      t_g2_q1ex=t_g2_q1ex
       
     9  +coeff( 45)*x11*x23            
     1  +coeff( 46)*x11*x22*x31        
     2  +coeff( 47)    *x24    *x41    
     3  +coeff( 48)    *x23*x31*x41    
     4  +coeff( 49)*x11*x22        *x51
     5  +coeff( 50)*x11*x24            
     6  +coeff( 51)*x11*x22    *x42    
     7  +coeff( 52)    *x23    *x43*x51
c
      return
      end
      function y_g2_q1ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 22)
      data ncoeff/ 21/
      data avdat/ -0.3170791E-01/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24706847E-02, 0.56421342E-02, 0.12213390E+00,-0.31966465E-02,
     +  0.16814681E-01, 0.40886072E-02,-0.22401031E-01,-0.60014538E-02,
     + -0.32343687E-02,-0.17143163E-02,-0.92659384E-03, 0.68256742E-03,
     + -0.99125225E-03, 0.10529424E-02,-0.26389454E-02, 0.30104807E-02,
     + -0.18684797E-03,-0.77192055E-03, 0.88117883E-03,-0.10847086E-02,
     +  0.97559585E-03,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_g2_q1ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x22    *x41    
      y_g2_q1ex=y_g2_q1ex
       
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)    *x24            
     8  +coeff( 17)*x11                
      y_g2_q1ex=y_g2_q1ex
       
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x23            
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)*x11*x23            
c
      return
      end
      function p_g2_q1ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 23)
      data ncoeff/ 22/
      data avdat/ -0.1323062E-01/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.62697753E-03,-0.16606236E-02, 0.19330748E-02, 0.53740501E-01,
     +  0.90716034E-02,-0.11752377E-01,-0.16682391E-02, 0.21508993E-02,
     + -0.15370152E-02,-0.31657328E-02,-0.48288939E-03,-0.57108066E-03,
     +  0.16505966E-02, 0.35745601E-03, 0.45276870E-03, 0.62310026E-03,
     + -0.56265196E-03,-0.14373109E-02,-0.44401825E-03,-0.94908341E-04,
     +  0.92683484E-04, 0.49444364E-03,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_g2_q1ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)            *x42    
      p_g2_q1ex=p_g2_q1ex
       
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x24            
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)*x11*x21    *x41    
      p_g2_q1ex=p_g2_q1ex
       
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x24*x31        
     2  +coeff( 20)*x11                
     3  +coeff( 21)    *x21        *x51
     4  +coeff( 22)*x11*x23            
c
      return
      end
      function l_g2_q1ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 23)
      data ncoeff/ 22/
      data avdat/ -0.4975448E-03/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15423987E-02, 0.18341512E-04,-0.28580700E-02,-0.32206581E-05,
     + -0.28177735E-02, 0.11426452E-03,-0.10764664E-03,-0.21531091E-02,
     + -0.66618022E-03, 0.69152637E-04,-0.17265367E-06, 0.12380927E-02,
     +  0.13322393E-03,-0.16959355E-04,-0.35077741E-03, 0.78206074E-04,
     + -0.62779407E-04, 0.10969207E-03,-0.13682789E-03,-0.13262259E-03,
     +  0.26856648E-03,-0.28027396E-03,
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
      l_g2_q1ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      l_g2_q1ex=l_g2_q1ex
       
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x22*x31        
     2  +coeff( 11)        *x33        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)    *x22*x33        
     6  +coeff( 15)        *x31        
     7  +coeff( 16)                *x51
     8  +coeff( 17)                *x52
      l_g2_q1ex=l_g2_q1ex
       
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)            *x43    
     2  +coeff( 20)    *x24            
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)    *x24    *x41    
c
      return
      end
      function x_g2_dent(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5096737E+01/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12237760E-02,-0.91166005E-01, 0.19846946E-01,-0.61752577E-02,
     +  0.51476941E-02, 0.14250478E-02, 0.19518358E-02,-0.28889384E-02,
     +  0.29269804E-02,-0.79325428E-02, 0.14178540E-02, 0.73242496E-03,
     + -0.23195827E-02,-0.46409314E-03,-0.43180017E-03,-0.10335656E-02,
     +  0.14340672E-03, 0.63504593E-03, 0.49298402E-03,-0.17256225E-02,
     +  0.27577585E-03, 0.36629278E-03,-0.59348351E-03, 0.13901960E-03,
     +  0.88559893E-04, 0.25570551E-02,-0.62616091E-04, 0.10626896E-03,
     + -0.76928780E-04,-0.67260367E-03, 0.59051299E-03, 0.10816435E-02,
     + -0.23666391E-03, 0.20895786E-03,-0.43812161E-04,-0.25478308E-03,
     + -0.52343804E-03,-0.91580063E-04,-0.37212097E-03,-0.81068327E-04,
     + -0.11458366E-03, 0.12052816E-03,-0.12258228E-03, 0.24044556E-03,
     +  0.19401683E-03, 0.50191622E-03, 0.37217480E-04, 0.33884597E-04,
     +  0.50474260E-04,-0.41239215E-04,-0.47437879E-04, 0.14109905E-03,
     + -0.34271160E-03,-0.49984512E-04, 0.13042981E-02, 0.68915397E-03,
     +  0.41357096E-04,-0.31022905E-03,-0.79552963E-03,-0.33645083E-04,
     +  0.28068600E-04,-0.11432755E-03, 0.11045118E-03, 0.55930013E-03,
     +  0.10074810E-03, 0.14839216E-04,-0.83000043E-04,-0.14405255E-03,
     +  0.28916734E-03,-0.18597321E-03, 0.14365466E-03,-0.87189343E-03,
     +  0.85875817E-03, 0.19626608E-04,-0.23708002E-04,-0.95800606E-04,
     +  0.30842355E-04, 0.18174222E-03, 0.33347293E-04,-0.17166416E-03,
     + -0.39138403E-03,-0.12789902E-03, 0.34820532E-04,-0.44818135E-03,
     +  0.35438588E-04, 0.73884126E-04,-0.39834798E-04,-0.18986483E-03,
     + -0.41906620E-04, 0.63896114E-04,
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
      x_g2_dent=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x23            
      x_g2_dent=x_g2_dent
  
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)*x11        *x41    
     3  +coeff( 12)    *x22            
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)        *x31        
      x_g2_dent=x_g2_dent
  
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)*x11*x22    *x41    
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)    *x24            
     5  +coeff( 23)    *x23*x31        
     6  +coeff( 24)*x11    *x31        
     7  +coeff( 25)    *x24*x31        
     8  +coeff( 26)    *x23    *x42    
      x_g2_dent=x_g2_dent
  
     9  +coeff( 27)                *x51
     1  +coeff( 28)            *x42    
     2  +coeff( 29)            *x41*x51
     3  +coeff( 30)    *x21    *x42*x51
     4  +coeff( 31)    *x24    *x41    
     5  +coeff( 32)    *x23    *x41*x51
     6  +coeff( 33)*x11*x21    *x41    
     7  +coeff( 34)*x11        *x42    
     8  +coeff( 35)    *x23        *x53
      x_g2_dent=x_g2_dent
  
     9  +coeff( 36)    *x22*x31        
     1  +coeff( 37)    *x21    *x41*x51
     2  +coeff( 38)    *x22    *x42    
     3  +coeff( 39)    *x23        *x51
     4  +coeff( 40)    *x21*x31*x41*x51
     5  +coeff( 41)    *x21    *x44*x51
     6  +coeff( 42)*x11*x23            
     7  +coeff( 43)*x11*x22*x31        
     8  +coeff( 44)*x11*x24            
      x_g2_dent=x_g2_dent
  
     9  +coeff( 45)*x11*x23    *x41    
     1  +coeff( 46)*x11*x24    *x41    
     2  +coeff( 47)        *x31*x41    
     3  +coeff( 48)    *x21*x32        
     4  +coeff( 49)            *x43    
     5  +coeff( 50)    *x21*x31    *x51
     6  +coeff( 51)    *x22*x31*x41    
     7  +coeff( 52)    *x21*x31*x42    
     8  +coeff( 53)    *x21    *x43    
      x_g2_dent=x_g2_dent
  
     9  +coeff( 54)    *x23*x31*x41    
     1  +coeff( 55)    *x21    *x44    
     2  +coeff( 56)    *x21    *x43*x51
     3  +coeff( 57)*x11    *x31*x41    
     4  +coeff( 58)    *x23*x31*x42    
     5  +coeff( 59)    *x23    *x43    
     6  +coeff( 60)*x11        *x41*x51
     7  +coeff( 61)*x11            *x52
     8  +coeff( 62)    *x23    *x41*x52
      x_g2_dent=x_g2_dent
  
     9  +coeff( 63)    *x24*x32*x41    
     1  +coeff( 64)    *x23*x31*x43    
     2  +coeff( 65)    *x23*x31*x42*x51
     3  +coeff( 66)    *x23    *x42*x52
     4  +coeff( 67)*x11*x21*x31*x42    
     5  +coeff( 68)    *x22*x33*x42*x51
     6  +coeff( 69)    *x24*x33*x42    
     7  +coeff( 70)    *x24*x31*x42*x52
     8  +coeff( 71)    *x21*x34*x42*x53
      x_g2_dent=x_g2_dent
  
     9  +coeff( 72)*x11*x22*x31*x44*x52
     1  +coeff( 73)*x11*x22*x33*x42*x54
     2  +coeff( 74)        *x31*x42    
     3  +coeff( 75)            *x42*x51
     4  +coeff( 76)    *x21    *x41*x52
     5  +coeff( 77)    *x22*x31*x41*x51
     6  +coeff( 78)    *x21    *x42*x52
     7  +coeff( 79)    *x21*x34*x41    
     8  +coeff( 80)    *x22    *x44    
      x_g2_dent=x_g2_dent
  
     9  +coeff( 81)    *x23    *x42*x51
     1  +coeff( 82)    *x22    *x42*x52
     2  +coeff( 83)    *x22        *x54
     3  +coeff( 84)    *x23    *x44    
     4  +coeff( 85)*x11    *x32    *x51
     5  +coeff( 86)    *x22*x32*x42*x51
     6  +coeff( 87)*x11*x22*x32        
     7  +coeff( 88)    *x24*x31*x43    
     8  +coeff( 89)    *x24*x33    *x51
      x_g2_dent=x_g2_dent
  
     9  +coeff( 90)*x11*x22    *x41*x51
c
      return
      end
      function t_g2_dent(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 65)
      data ncoeff/ 64/
      data avdat/  0.3591255E+01/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.30166306E-01, 0.54197001E+00,-0.13744401E-01, 0.64865779E-02,
     + -0.41489746E-01, 0.12055094E+00,-0.14327842E+00, 0.39569758E-01,
     +  0.44971790E-01,-0.14047395E-01,-0.13290204E-01,-0.49319498E-01,
     + -0.10811482E-02, 0.82294056E-02,-0.94535071E-02,-0.20718906E-01,
     +  0.13558418E-01, 0.50195385E-01, 0.31021112E-02,-0.48683495E-02,
     +  0.29558169E-02, 0.68955719E-02, 0.12912577E-01,-0.34282135E-02,
     + -0.40058414E-02,-0.61110901E-02, 0.78973053E-02, 0.20740772E-01,
     + -0.86518831E-03, 0.31854797E-02, 0.59151149E-03, 0.29944049E-02,
     +  0.35712083E-02, 0.12386102E-02,-0.57503083E-02, 0.77202549E-03,
     + -0.21567063E-01,-0.79742793E-04, 0.63189439E-03, 0.26471014E-02,
     +  0.95067348E-03,-0.14253287E-02, 0.65678502E-02,-0.12492167E-02,
     + -0.77589089E-02,-0.10984010E-01, 0.78175002E-03,-0.57494018E-03,
     +  0.38503085E-02, 0.36786363E-03, 0.94466860E-03, 0.74003957E-03,
     +  0.43792962E-02,-0.93101058E-02, 0.24672525E-02,-0.49088080E-03,
     + -0.39691911E-02,-0.11883003E-02, 0.43073285E-03, 0.98912208E-03,
     + -0.83538820E-03, 0.28450480E-02,-0.34860137E-02, 0.62773754E-02,
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_g2_dent=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_g2_dent=t_g2_dent
  
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)        *x31        
     5  +coeff( 14)            *x42    
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x22        *x51
      t_g2_dent=t_g2_dent
  
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)    *x24            
     5  +coeff( 23)*x11*x22    *x41    
     6  +coeff( 24)    *x22*x31        
     7  +coeff( 25)    *x21        *x52
     8  +coeff( 26)    *x22    *x42    
      t_g2_dent=t_g2_dent
  
     9  +coeff( 27)    *x23        *x51
     1  +coeff( 28)    *x24    *x41    
     2  +coeff( 29)*x11    *x31        
     3  +coeff( 30)            *x42*x51
     4  +coeff( 31)*x12                
     5  +coeff( 32)*x11*x22            
     6  +coeff( 33)    *x23*x31        
     7  +coeff( 34)*x11*x21    *x41    
     8  +coeff( 35)    *x22    *x41*x51
      t_g2_dent=t_g2_dent
  
     9  +coeff( 36)*x12*x21            
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)                *x52
     3  +coeff( 39)            *x43    
     4  +coeff( 40)    *x21    *x41*x51
     5  +coeff( 41)            *x41*x52
     6  +coeff( 42)*x11        *x42    
     7  +coeff( 43)    *x21    *x42*x51
     8  +coeff( 44)    *x22        *x52
      t_g2_dent=t_g2_dent
  
     9  +coeff( 45)    *x23    *x41*x51
     1  +coeff( 46)    *x24    *x42    
     2  +coeff( 47)*x12*x21    *x43    
     3  +coeff( 48)    *x22*x31*x41    
     4  +coeff( 49)    *x21    *x43    
     5  +coeff( 50)    *x21*x31*x41*x51
     6  +coeff( 51)    *x21    *x41*x52
     7  +coeff( 52)*x11*x22*x31        
     8  +coeff( 53)    *x22    *x43    
      t_g2_dent=t_g2_dent
  
     9  +coeff( 54)    *x21    *x44    
     1  +coeff( 55)    *x24        *x51
     2  +coeff( 56)*x11*x21    *x41*x51
     3  +coeff( 57)    *x21    *x43*x51
     4  +coeff( 58)*x11*x24            
     5  +coeff( 59)*x11*x23*x31        
     6  +coeff( 60)*x11*x23    *x41    
     7  +coeff( 61)    *x23*x31*x42    
     8  +coeff( 62)    *x23    *x43    
      t_g2_dent=t_g2_dent
  
     9  +coeff( 63)*x11*x24    *x41    
     1  +coeff( 64)    *x23    *x44    
c
      return
      end
      function y_g2_dent(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 32)
      data ncoeff/ 31/
      data avdat/ -0.1601798E-01/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.53396431E-03, 0.14910650E-02, 0.74911319E-01,-0.40905993E-02,
     +  0.12178769E-01, 0.46386705E-02, 0.32818567E-02, 0.45744391E-03,
     +  0.12615356E-01,-0.20530028E-01, 0.15645454E-02,-0.80270162E-02,
     + -0.29548118E-02,-0.21221708E-02, 0.72486338E-03, 0.33785603E-02,
     + -0.12440389E-03,-0.77157223E-03,-0.79411420E-03,-0.75312634E-03,
     + -0.85968454E-03,-0.85228268E-03, 0.12114133E-03, 0.30039620E-03,
     + -0.32696506E-03,-0.24215950E-03,-0.39247860E-03, 0.31737829E-03,
     + -0.13900037E-02,-0.26306428E-03, 0.90890675E-03,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_g2_dent=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31    *x51
      y_g2_dent=y_g2_dent
  
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x22            
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)    *x24            
     8  +coeff( 17)*x11                
      y_g2_dent=y_g2_dent
  
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)            *x41*x52
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x21*x31        
     6  +coeff( 24)    *x21        *x51
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)*x11        *x41    
      y_g2_dent=y_g2_dent
  
     9  +coeff( 27)    *x21    *x41*x51
     1  +coeff( 28)    *x23            
     2  +coeff( 29)    *x22    *x42    
     3  +coeff( 30)*x11*x21        *x51
     4  +coeff( 31)*x11*x23            
c
      return
      end
      function p_g2_dent(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/  0.2268176E-01/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12233705E-02, 0.62960829E-02,-0.46519842E-02,-0.76753736E-01,
     + -0.10761118E-01, 0.12679632E-01,-0.18090026E-01,-0.74196642E-03,
     + -0.27143310E-02, 0.12739912E-01, 0.23896082E-02, 0.11062548E-02,
     + -0.60107079E-02,-0.35842359E-02,-0.18624566E-02,-0.79696433E-03,
     +  0.52356103E-03, 0.29589578E-02, 0.11517716E-02, 0.37866300E-02,
     +  0.11114078E-02, 0.37959849E-02,-0.10335072E-02,-0.81959367E-03,
     + -0.13963006E-03,-0.13702587E-03, 0.50503062E-03, 0.20308203E-03,
     +  0.11517750E-02,-0.98794501E-03,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_g2_dent=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)            *x42    
      p_g2_dent=p_g2_dent
  
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x23    *x42    
     7  +coeff( 16)    *x21*x31        
     8  +coeff( 17)        *x31    *x51
      p_g2_dent=p_g2_dent
  
     9  +coeff( 18)    *x23            
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)            *x43    
     6  +coeff( 24)            *x41*x52
     7  +coeff( 25)*x11                
     8  +coeff( 26)        *x31*x41    
      p_g2_dent=p_g2_dent
  
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)*x11            *x51
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)    *x22    *x41*x51
c
      return
      end
      function l_g2_dent(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 17)
      data ncoeff/ 16/
      data avdat/ -0.4085935E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.54811980E-02, 0.17640470E+00,-0.40909559E-02,-0.95035229E-02,
     + -0.87489346E-02,-0.37135340E-01, 0.12515541E-01,-0.33246626E-02,
     +  0.48935828E-02,-0.44724578E-02,-0.13636403E-02, 0.12809893E-02,
     + -0.18866789E-02, 0.69965883E-02,-0.81427023E-02, 0.13733700E-01,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_g2_dent=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x21*x31        
      l_g2_dent=l_g2_dent
  
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)*x11        *x41    
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x23    *x41    
c
      return
      end
      function x_g2_dext(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.4784906E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12346556E-02, 0.24329734E+00,-0.34455706E-02, 0.11960944E+00,
     + -0.60741991E-01, 0.34219749E-01,-0.18456167E-01,-0.10252259E-02,
     + -0.60378797E-02, 0.98928753E-02, 0.10189338E-01,-0.27928366E-02,
     + -0.10704259E-01, 0.25534768E-01,-0.43537603E-02,-0.21917566E-02,
     +  0.27999929E-02,-0.24917754E-02, 0.30322059E-02,-0.39374229E-03,
     + -0.17471531E-02, 0.54318337E-02,-0.61270024E-03, 0.50363096E-03,
     + -0.76955301E-03, 0.19100563E-02, 0.21485738E-02,-0.42016848E-03,
     + -0.91964379E-02, 0.27826350E-03,-0.26452765E-03,-0.16892180E-03,
     +  0.81445987E-03,-0.57314662E-03, 0.38997806E-03,-0.73704415E-03,
     +  0.16325749E-03,-0.28286022E-03,-0.18291610E-03, 0.17068901E-02,
     +  0.20691479E-03, 0.18683173E-02,-0.28774986E-04,-0.13873095E-03,
     + -0.46799607E-02,-0.21371520E-02,-0.22476462E-02, 0.23166637E-02,
     +  0.14266276E-03, 0.19404017E-03, 0.43608720E-03,-0.53285476E-03,
     + -0.15775763E-02, 0.39210929E-04, 0.47566154E-03,-0.63105737E-03,
     + -0.81750681E-03,-0.16451959E-03,-0.54357538E-03, 0.61393199E-04,
     + -0.10056163E-03,-0.11179807E-03, 0.10038059E-02,-0.52003440E-03,
     + -0.12603245E-03,-0.30037574E-03, 0.24146819E-02,-0.55025239E-03,
     + -0.13977827E-02, 0.36794777E-03,-0.13592822E-02, 0.13533800E-02,
     + -0.12316639E-02,-0.71360490E-04, 0.55314107E-04, 0.35237501E-03,
     +  0.29073225E-03, 0.53338124E-04,-0.68876478E-04,-0.33648677E-04,
     +  0.21722438E-03, 0.13821578E-03, 0.34206506E-03,-0.34374616E-03,
     +  0.94232353E-04,-0.35730330E-03,-0.35457825E-03, 0.51111507E-03,
     +  0.34735122E-03, 0.23431241E-03,
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
      x_g2_dext=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x24            
      x_g2_dext=x_g2_dext
  
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x22            
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)            *x42    
     8  +coeff( 17)    *x22    *x41    
      x_g2_dext=x_g2_dext
  
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)        *x31        
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)            *x41*x51
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)    *x21        *x52
     8  +coeff( 26)    *x23*x31        
      x_g2_dext=x_g2_dext
  
     9  +coeff( 27)    *x23        *x51
     1  +coeff( 28)*x11    *x31        
     2  +coeff( 29)    *x23    *x42    
     3  +coeff( 30)*x11            *x51
     4  +coeff( 31)        *x31*x41    
     5  +coeff( 32)    *x21*x31    *x51
     6  +coeff( 33)*x11*x21    *x41    
     7  +coeff( 34)*x11        *x42    
     8  +coeff( 35)*x11*x22*x31        
      x_g2_dext=x_g2_dext
  
     9  +coeff( 36)*x11*x24            
     1  +coeff( 37)    *x22*x31        
     2  +coeff( 38)            *x43    
     3  +coeff( 39)    *x21*x31*x42    
     4  +coeff( 40)    *x21    *x43    
     5  +coeff( 41)    *x22    *x41*x51
     6  +coeff( 42)    *x21    *x42*x51
     7  +coeff( 43)*x11*x21            
     8  +coeff( 44)    *x21*x34        
      x_g2_dext=x_g2_dext
  
     9  +coeff( 45)    *x21    *x44    
     1  +coeff( 46)    *x23    *x41*x51
     2  +coeff( 47)    *x21    *x43*x51
     3  +coeff( 48)    *x23    *x43    
     4  +coeff( 49)*x11*x21        *x51
     5  +coeff( 50)    *x23    *x41*x52
     6  +coeff( 51)*x11*x22        *x51
     7  +coeff( 52)*x11*x23    *x41    
     8  +coeff( 53)*x11*x24    *x41    
      x_g2_dext=x_g2_dext
  
     9  +coeff( 54)        *x31*x41*x51
     1  +coeff( 55)    *x21    *x41*x52
     2  +coeff( 56)    *x23*x31*x41    
     3  +coeff( 57)    *x21*x31*x43    
     4  +coeff( 58)            *x44*x51
     5  +coeff( 59)    *x21    *x42*x52
     6  +coeff( 60)*x11*x21*x31        
     7  +coeff( 61)*x11    *x31*x41    
     8  +coeff( 62)*x11        *x41*x51
      x_g2_dext=x_g2_dext
  
     9  +coeff( 63)    *x23    *x42*x51
     1  +coeff( 64)    *x22        *x54
     2  +coeff( 65)*x11*x23            
     3  +coeff( 66)    *x24*x32*x41    
     4  +coeff( 67)    *x23    *x44    
     5  +coeff( 68)    *x21*x32*x42*x53
     6  +coeff( 69)*x11*x22*x32    *x51
     7  +coeff( 70)*x11*x22*x34    *x51
     8  +coeff( 71)    *x24*x34*x42*x51
      x_g2_dext=x_g2_dext
  
     9  +coeff( 72)*x11*x24*x34    *x51
     1  +coeff( 73)*x11*x22*x32*x43*x52
     2  +coeff( 74)        *x31*x42    
     3  +coeff( 75)                *x53
     4  +coeff( 76)    *x22*x31*x41    
     5  +coeff( 77)    *x22        *x52
     6  +coeff( 78)        *x31*x41*x52
     7  +coeff( 79)            *x41*x53
     8  +coeff( 80)    *x22*x33        
      x_g2_dext=x_g2_dext
  
     9  +coeff( 81)    *x22    *x43    
     1  +coeff( 82)    *x22*x32    *x51
     2  +coeff( 83)    *x22    *x42*x51
     3  +coeff( 84)    *x21*x31*x42*x51
     4  +coeff( 85)    *x22*x31    *x52
     5  +coeff( 86)    *x22*x33*x41    
     6  +coeff( 87)    *x24    *x42    
     7  +coeff( 88)    *x23*x31*x42    
     8  +coeff( 89)    *x22    *x44    
      x_g2_dext=x_g2_dext
  
     9  +coeff( 90)    *x23*x31*x41*x51
c
      return
      end
      function t_g2_dext(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 57)
      data ncoeff/ 56/
      data avdat/  0.5946798E+00/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.82595600E-03,-0.77001482E-01, 0.27706500E-01, 0.42968364E-02,
     +  0.16858380E-01,-0.66988054E-02, 0.14912280E-02, 0.16944426E-02,
     + -0.28722091E-02,-0.86190109E-03,-0.10594068E-02, 0.11595305E-02,
     +  0.25270991E-02,-0.63645663E-02,-0.25118873E-03,-0.38385432E-03,
     + -0.89806842E-03,-0.82317746E-03, 0.17559310E-03, 0.47591218E-03,
     + -0.25516632E-03, 0.44228783E-03, 0.41646860E-03,-0.14047907E-02,
     +  0.11814162E-03, 0.35644067E-03,-0.52722852E-03, 0.14537034E-02,
     + -0.78740675E-04,-0.56386376E-04,-0.97085984E-04, 0.63272899E-04,
     +  0.94878982E-04, 0.17396288E-03,-0.26385422E-03, 0.17910017E-03,
     + -0.34126345E-03, 0.10596903E-03,-0.26433967E-03,-0.79495461E-04,
     + -0.21176347E-04, 0.93279101E-04, 0.86625521E-04,-0.20134170E-03,
     + -0.90739479E-04, 0.83794926E-04,-0.98675817E-04, 0.20247228E-03,
     +  0.57733018E-03,-0.17481667E-03, 0.30297556E-03, 0.18041275E-03,
     +  0.17621883E-03, 0.16454696E-03,-0.54217211E-03, 0.40670997E-03,
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
      t_g2_dext=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)            *x41    
     8  +coeff(  8)    *x21*x31        
      t_g2_dext=t_g2_dext
  
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x22            
     2  +coeff( 11)                *x52
     3  +coeff( 12)*x11        *x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)            *x42    
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)    *x22    *x41    
      t_g2_dext=t_g2_dext
  
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)        *x31        
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)*x11*x22    *x41    
     7  +coeff( 25)*x11    *x31        
     8  +coeff( 26)    *x22        *x51
      t_g2_dext=t_g2_dext
  
     9  +coeff( 27)    *x23*x31        
     1  +coeff( 28)    *x23    *x42    
     2  +coeff( 29)        *x31*x41    
     3  +coeff( 30)*x11*x21            
     4  +coeff( 31)    *x22*x31        
     5  +coeff( 32)    *x21    *x41*x51
     6  +coeff( 33)        *x31*x41*x51
     7  +coeff( 34)    *x24            
     8  +coeff( 35)*x11*x21    *x41    
      t_g2_dext=t_g2_dext
  
     9  +coeff( 36)*x11        *x42    
     1  +coeff( 37)    *x23        *x51
     2  +coeff( 38)    *x22        *x52
     3  +coeff( 39)            *x42*x52
     4  +coeff( 40)*x12*x21            
     5  +coeff( 41)        *x31    *x51
     6  +coeff( 42)            *x41*x52
     7  +coeff( 43)                *x53
     8  +coeff( 44)    *x22    *x42    
      t_g2_dext=t_g2_dext
  
     9  +coeff( 45)            *x41*x53
     1  +coeff( 46)*x11*x23            
     2  +coeff( 47)*x11*x22*x31        
     3  +coeff( 48)    *x24    *x41    
     4  +coeff( 49)    *x21    *x44    
     5  +coeff( 50)    *x24        *x51
     6  +coeff( 51)    *x23    *x41*x51
     7  +coeff( 52)*x11*x24            
     8  +coeff( 53)*x11*x23    *x41    
      t_g2_dext=t_g2_dext
  
     9  +coeff( 54)    *x23*x31*x42    
     1  +coeff( 55)    *x23    *x43    
     2  +coeff( 56)*x11*x24    *x41    
c
      return
      end
      function y_g2_dext(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 56)
      data ncoeff/ 55/
      data avdat/  0.1680168E-01/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.20951123E-02,-0.63198674E-02,-0.25058087E-01, 0.46027191E-02,
     +  0.48247082E-03, 0.74997230E-03, 0.54112165E-02,-0.27061358E-01,
     +  0.16392533E-02, 0.41757047E-01,-0.50747750E-03,-0.11160618E-01,
     + -0.46072975E-02, 0.64945198E-02, 0.35963375E-02,-0.22161477E-03,
     +  0.14864587E-02,-0.16811809E-01,-0.30658303E-02,-0.29098694E-02,
     + -0.20667310E-02, 0.62732124E-02,-0.95887873E-02, 0.15857106E-04,
     + -0.16499730E-02, 0.53178789E-02, 0.16445885E-03,-0.86660427E-03,
     + -0.95230970E-03, 0.23122702E-03,-0.65019302E-03, 0.47415439E-02,
     + -0.85630233E-03,-0.31861715E-03,-0.29634584E-02, 0.44509687E-03,
     + -0.40354161E-03, 0.85712172E-03, 0.54808729E-03, 0.63162402E-03,
     + -0.25468771E-02, 0.74033625E-03,-0.15130045E-03, 0.26975348E-03,
     +  0.34341693E-03, 0.15361868E-03, 0.57970843E-03, 0.51932398E-03,
     + -0.15112117E-02,-0.98457257E-03,-0.17270279E-02,-0.44322677E-03,
     +  0.76940912E-03, 0.35322618E-03, 0.57187775E-03,
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
      x53 = x52*x5
c
c                  function
c
      y_g2_dext=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21    *x41    
      y_g2_dext=y_g2_dext
  
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x22            
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)*x11        *x41    
      y_g2_dext=y_g2_dext
  
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)    *x23            
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x21*x33        
     7  +coeff( 25)    *x22    *x41*x51
     8  +coeff( 26)    *x24            
      y_g2_dext=y_g2_dext
  
     9  +coeff( 27)            *x45    
     1  +coeff( 28)    *x21*x31        
     2  +coeff( 29)    *x22*x31        
     3  +coeff( 30)*x11            *x51
     4  +coeff( 31)                *x53
     5  +coeff( 32)    *x22    *x42    
     6  +coeff( 33)*x11*x21        *x51
     7  +coeff( 34)        *x31*x42    
     8  +coeff( 35)            *x43    
      y_g2_dext=y_g2_dext
  
     9  +coeff( 36)    *x21*x31*x41    
     1  +coeff( 37)    *x21        *x52
     2  +coeff( 38)    *x21    *x43    
     3  +coeff( 39)*x11*x22            
     4  +coeff( 40)    *x23        *x51
     5  +coeff( 41)    *x24    *x41    
     6  +coeff( 42)*x11*x23            
     7  +coeff( 43)        *x31    *x52
     8  +coeff( 44)    *x22*x31*x41    
      y_g2_dext=y_g2_dext
  
     9  +coeff( 45)*x11        *x42    
     1  +coeff( 46)*x11*x21    *x41    
     2  +coeff( 47)    *x22        *x52
     3  +coeff( 48)    *x21    *x44    
     4  +coeff( 49)    *x22    *x43    
     5  +coeff( 50)    *x23    *x42    
     6  +coeff( 51)    *x22    *x42*x51
     7  +coeff( 52)*x11*x22    *x41    
     8  +coeff( 53)    *x24        *x51
      y_g2_dext=y_g2_dext
  
     9  +coeff( 54)            *x45*x51
     1  +coeff( 55)    *x23    *x41*x52
c
      return
      end
      function p_g2_dext(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 40)
      data ncoeff/ 39/
      data avdat/  0.6917529E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.66459551E-03, 0.13265443E-02,-0.18815013E-01,-0.17686719E-02,
     + -0.56516863E-02, 0.63579855E-03,-0.11519421E-02, 0.28809131E-03,
     +  0.84610321E-02, 0.14060489E-02,-0.22312224E-04, 0.27807882E-05,
     + -0.26885266E-02,-0.20533290E-02, 0.75357349E-03,-0.33000084E-04,
     + -0.19043508E-02, 0.11938807E-02, 0.30829455E-03, 0.77621499E-03,
     + -0.80352061E-03, 0.38275539E-03,-0.59787999E-03, 0.10449941E-02,
     + -0.70107938E-03,-0.13878195E-03, 0.76483135E-04, 0.63188010E-04,
     + -0.19307973E-03, 0.14907411E-03,-0.19428738E-03, 0.33531908E-03,
     + -0.53460310E-04,-0.14965274E-03,-0.88066758E-04, 0.86339656E-04,
     +  0.28634179E-03,-0.45873344E-03, 0.12560659E-03,
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
      p_g2_dext=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)        *x31    *x51
      p_g2_dext=p_g2_dext
  
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)        *x33        
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x24            
     7  +coeff( 16)    *x22*x33        
     8  +coeff( 17)        *x31        
      p_g2_dext=p_g2_dext
  
     9  +coeff( 18)    *x23            
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)    *x22            
     5  +coeff( 23)            *x43    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)    *x22    *x41*x51
     8  +coeff( 26)    *x21*x31        
      p_g2_dext=p_g2_dext
  
     9  +coeff( 27)        *x31*x41    
     1  +coeff( 28)*x11            *x51
     2  +coeff( 29)            *x41*x52
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)*x11*x21        *x51
     5  +coeff( 32)    *x23        *x51
     6  +coeff( 33)*x11                
     7  +coeff( 34)    *x21        *x52
     8  +coeff( 35)                *x53
      p_g2_dext=p_g2_dext
  
     9  +coeff( 36)*x11*x22            
     1  +coeff( 37)    *x23    *x41    
     2  +coeff( 38)    *x24    *x41    
     3  +coeff( 39)    *x23*x31*x41    
c
      return
      end
      function l_g2_dext(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/ -0.6000567E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.48015579E-02,-0.34968928E+00, 0.44285110E-02,-0.88722140E-01,
     +  0.22835551E-01,-0.23215832E-01, 0.81579715E-01,-0.36787737E-01,
     +  0.73368005E-02,-0.13274551E-01, 0.38492084E-02, 0.18747026E-01,
     + -0.29701361E-01,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_g2_dext=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      l_g2_dext=l_g2_dext
  
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x23            
     2  +coeff( 11)*x11        *x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x23    *x41    
c
      return
      end
      function x_g2_q3en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.6406332E-03/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.20132267E-02, 0.15235345E+00,-0.19684969E-02, 0.12446510E+00,
     +  0.12673044E-01,-0.39962895E-01, 0.26231455E-01,-0.12710427E-01,
     + -0.47352444E-02, 0.69713872E-02, 0.12496504E-02,-0.19164070E-02,
     + -0.79649724E-02, 0.16923841E-01,-0.28663010E-02,-0.39846338E-02,
     +  0.20163788E-02,-0.26983719E-02, 0.18982055E-02,-0.25096245E-03,
     + -0.67151099E-03,-0.14211886E-02,-0.67676773E-03, 0.36533515E-02,
     +  0.18562636E-02,-0.27924692E-03,-0.53449296E-02,-0.19248319E-03,
     + -0.19518571E-03,-0.29223875E-03, 0.24082647E-03, 0.52427192E-03,
     + -0.36425394E-03,-0.26123988E-03,-0.39949303E-03, 0.12203925E-02,
     +  0.15960161E-02,-0.24943589E-02, 0.11810172E-04,-0.13925608E-02,
     +  0.18627113E-02, 0.26080685E-03,-0.12003812E-02, 0.27637056E-03,
     + -0.15801797E-03,-0.48935448E-03,-0.18746212E-03,-0.11047482E-02,
     +  0.72050031E-03,-0.71803435E-04,-0.87419315E-03,-0.47703586E-04,
     +  0.15377993E-03, 0.22964802E-03, 0.23363797E-04, 0.14075630E-02,
     + -0.10195444E-03,-0.98339678E-03,-0.30679192E-03, 0.31948250E-04,
     + -0.79415833E-04,-0.11055073E-02,-0.75498472E-04, 0.99833369E-05,
     + -0.30731040E-03,-0.16469350E-03,-0.20274749E-04,-0.22869950E-03,
     + -0.56871137E-03,-0.48647160E-03, 0.13053223E-02, 0.20438579E-02,
     + -0.20087587E-04,-0.60908358E-04, 0.90105335E-04,-0.23874374E-03,
     +  0.13258024E-03, 0.27670790E-03,-0.52356940E-04, 0.69976086E-04,
     + -0.67580484E-04, 0.78765315E-03, 0.36830588E-04, 0.21424961E-03,
     + -0.17000080E-03,-0.71066563E-04,-0.26109090E-03,-0.37750628E-03,
     +  0.49225957E-03,-0.17983170E-03,
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
      x_g2_q3en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)*x11                
      x_g2_q3en=x_g2_q3en
  
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x23*x31        
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x21*x31        
     8  +coeff( 17)    *x22        *x51
      x_g2_q3en=x_g2_q3en
  
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)        *x31        
     3  +coeff( 21)            *x41*x51
     4  +coeff( 22)    *x21*x31*x41    
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)*x11*x22    *x41    
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)*x11    *x31        
      x_g2_q3en=x_g2_q3en
  
     9  +coeff( 27)    *x23    *x42    
     1  +coeff( 28)        *x31*x41    
     2  +coeff( 29)    *x21*x31    *x51
     3  +coeff( 30)    *x21        *x52
     4  +coeff( 31)                *x53
     5  +coeff( 32)*x11*x21    *x41    
     6  +coeff( 33)*x11        *x42    
     7  +coeff( 34)            *x43    
     8  +coeff( 35)    *x24            
      x_g2_q3en=x_g2_q3en
  
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)    *x21    *x42*x51
     2  +coeff( 38)    *x21    *x44    
     3  +coeff( 39)*x11            *x51
     4  +coeff( 40)    *x21    *x43*x51
     5  +coeff( 41)    *x23    *x43    
     6  +coeff( 42)*x11*x22*x31        
     7  +coeff( 43)    *x23*x31*x43    
     8  +coeff( 44)*x11*x22        *x51
      x_g2_q3en=x_g2_q3en
  
     9  +coeff( 45)    *x21*x34    *x52
     1  +coeff( 46)*x11*x24            
     2  +coeff( 47)*x11        *x43*x51
     3  +coeff( 48)*x11*x24    *x41    
     4  +coeff( 49)    *x21*x32*x43*x54
     5  +coeff( 50)    *x22*x31        
     6  +coeff( 51)    *x22    *x41    
     7  +coeff( 52)        *x31*x42    
     8  +coeff( 53)    *x22*x31*x41    
      x_g2_q3en=x_g2_q3en
  
     9  +coeff( 54)    *x21*x32*x41    
     1  +coeff( 55)    *x22    *x42    
     2  +coeff( 56)    *x24    *x41    
     3  +coeff( 57)    *x22*x32*x41    
     4  +coeff( 58)    *x23    *x41*x51
     5  +coeff( 59)    *x21    *x42*x52
     6  +coeff( 60)*x11*x21*x31        
     7  +coeff( 61)*x11    *x31*x41    
     8  +coeff( 62)    *x24    *x42    
      x_g2_q3en=x_g2_q3en
  
     9  +coeff( 63)    *x24    *x41*x51
     1  +coeff( 64)    *x23    *x41*x52
     2  +coeff( 65)    *x21    *x42*x53
     3  +coeff( 66)    *x24    *x43    
     4  +coeff( 67)*x11    *x32    *x51
     5  +coeff( 68)*x11*x23    *x41    
     6  +coeff( 69)    *x23*x34    *x53
     7  +coeff( 70)*x11*x22*x32    *x53
     8  +coeff( 71)*x11*x24*x31*x43*x52
      x_g2_q3en=x_g2_q3en
  
     9  +coeff( 72)*x11*x23*x33*x44*x52
     1  +coeff( 73)        *x31    *x51
     2  +coeff( 74)    *x21*x32        
     3  +coeff( 75)    *x21*x32    *x51
     4  +coeff( 76)    *x22    *x41*x51
     5  +coeff( 77)    *x21*x31*x41*x51
     6  +coeff( 78)    *x21    *x41*x52
     7  +coeff( 79)            *x42*x52
     8  +coeff( 80)    *x21        *x53
      x_g2_q3en=x_g2_q3en
  
     9  +coeff( 81)                *x54
     1  +coeff( 82)    *x22    *x43    
     2  +coeff( 83)        *x33*x41*x51
     3  +coeff( 84)    *x22    *x42*x51
     4  +coeff( 85)    *x21*x31*x42*x51
     5  +coeff( 86)            *x44*x51
     6  +coeff( 87)    *x22*x33*x41    
     7  +coeff( 88)    *x21*x34*x41    
     8  +coeff( 89)    *x23*x31*x42    
      x_g2_q3en=x_g2_q3en
  
     9  +coeff( 90)    *x21*x33*x42    
c
      return
      end
      function t_g2_q3en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 57)
      data ncoeff/ 56/
      data avdat/  0.5616265E-03/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13853516E-03,-0.57295226E-01, 0.20468943E-01, 0.32094033E-02,
     + -0.31712255E-02, 0.12366900E-01,-0.42252922E-02, 0.91133284E-03,
     +  0.12450112E-02,-0.11349046E-02,-0.20866771E-02, 0.87690604E-03,
     +  0.24848629E-02,-0.49375058E-02, 0.23085689E-03,-0.62120665E-03,
     +  0.78198420E-04, 0.43341477E-03,-0.20083238E-03,-0.10917373E-03,
     +  0.41289962E-03,-0.10077258E-02, 0.15952419E-04, 0.85279375E-04,
     +  0.16236863E-03,-0.37675843E-03, 0.95600006E-03, 0.40034673E-04,
     + -0.11656649E-04,-0.45102934E-04, 0.85994994E-04,-0.11519941E-03,
     +  0.11494390E-03, 0.24468981E-03,-0.41265349E-03,-0.40657764E-04,
     +  0.49957173E-03, 0.72940551E-04, 0.39362905E-04,-0.92461662E-04,
     + -0.15759126E-04,-0.11417341E-03,-0.48738197E-03, 0.51744999E-04,
     + -0.73663105E-04, 0.95884621E-04,-0.57394012E-04, 0.78191835E-04,
     + -0.67449160E-04,-0.16221008E-03,-0.47243171E-04,-0.12363672E-03,
     +  0.12992301E-03, 0.27551726E-03, 0.24192188E-03,-0.11197652E-03,
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
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_g2_q3en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)            *x41    
      t_g2_q3en=t_g2_q3en
  
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x23            
     3  +coeff( 12)*x11        *x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)        *x31        
      t_g2_q3en=t_g2_q3en
  
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)            *x42    
     6  +coeff( 24)*x11    *x31        
     7  +coeff( 25)    *x22    *x41    
     8  +coeff( 26)    *x23*x31        
      t_g2_q3en=t_g2_q3en
  
     9  +coeff( 27)    *x23    *x42    
     1  +coeff( 28)        *x31*x41    
     2  +coeff( 29)            *x43    
     3  +coeff( 30)    *x21    *x41*x51
     4  +coeff( 31)                *x53
     5  +coeff( 32)*x11*x21    *x41    
     6  +coeff( 33)*x11        *x42    
     7  +coeff( 34)    *x22    *x42    
     8  +coeff( 35)    *x23        *x51
      t_g2_q3en=t_g2_q3en
  
     9  +coeff( 36)*x12*x21            
     1  +coeff( 37)    *x23    *x41*x51
     2  +coeff( 38)    *x23    *x42*x51
     3  +coeff( 39)    *x22*x31        
     4  +coeff( 40)            *x42*x51
     5  +coeff( 41)*x12                
     6  +coeff( 42)    *x24            
     7  +coeff( 43)    *x21    *x42*x51
     8  +coeff( 44)            *x43*x51
      t_g2_q3en=t_g2_q3en
  
     9  +coeff( 45)    *x22        *x52
     1  +coeff( 46)            *x42*x52
     2  +coeff( 47)    *x21        *x53
     3  +coeff( 48)*x11*x23            
     4  +coeff( 49)*x11*x22*x31        
     5  +coeff( 50)    *x24        *x51
     6  +coeff( 51)        *x33*x41*x51
     7  +coeff( 52)    *x23        *x52
     8  +coeff( 53)*x11*x24            
      t_g2_q3en=t_g2_q3en
  
     9  +coeff( 54)    *x21    *x44*x51
     1  +coeff( 55)*x11*x24    *x41    
     2  +coeff( 56)*x12*x21    *x43    
c
      return
      end
      function y_g2_q3en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 52)
      data ncoeff/ 51/
      data avdat/  0.2350425E-01/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28556697E-02,-0.82641067E-02,-0.42674217E-01, 0.50821416E-02,
     + -0.15081099E-02, 0.89479907E-03, 0.61560981E-02,-0.32230701E-01,
     +  0.20462053E-02, 0.52197583E-01,-0.11522265E-01,-0.58681928E-02,
     +  0.80964230E-02, 0.47419514E-02,-0.20440957E-03, 0.17329210E-02,
     + -0.18885517E-01,-0.40978398E-02,-0.37024494E-02,-0.21221333E-02,
     +  0.78523690E-02,-0.11870306E-01,-0.22352384E-02, 0.61446433E-02,
     +  0.41383336E-03,-0.81831875E-03,-0.53439982E-03,-0.96029392E-03,
     +  0.56832191E-02,-0.10906468E-02,-0.37318745E-02, 0.28524463E-03,
     + -0.85454888E-03, 0.64732652E-03, 0.11208754E-02,-0.31446391E-02,
     +  0.76325011E-03,-0.75073325E-03,-0.30167575E-03, 0.44568777E-03,
     + -0.17170837E-03, 0.17259569E-02, 0.37881936E-03, 0.21419881E-03,
     +  0.18530573E-02, 0.81849180E-03,-0.20155334E-02,-0.26316585E-02,
     + -0.18769072E-02,-0.51467604E-03, 0.96328062E-03,
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
      x31 = x3
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
      y_g2_q3en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21    *x41    
      y_g2_q3en=y_g2_q3en
  
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)    *x22            
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)            *x42*x51
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x22    *x41    
      y_g2_q3en=y_g2_q3en
  
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)    *x24            
     7  +coeff( 25)            *x45    
     8  +coeff( 26)    *x21*x31        
      y_g2_q3en=y_g2_q3en
  
     9  +coeff( 27)*x11                
     1  +coeff( 28)    *x22*x31        
     2  +coeff( 29)    *x22    *x42    
     3  +coeff( 30)*x11*x21        *x51
     4  +coeff( 31)            *x43    
     5  +coeff( 32)*x11            *x51
     6  +coeff( 33)                *x53
     7  +coeff( 34)*x11*x22            
     8  +coeff( 35)    *x23        *x51
      y_g2_q3en=y_g2_q3en
  
     9  +coeff( 36)    *x24    *x41    
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)    *x23        *x52
     3  +coeff( 39)        *x31*x42    
     4  +coeff( 40)    *x21*x31*x41    
     5  +coeff( 41)        *x31    *x52
     6  +coeff( 42)    *x21    *x43    
     7  +coeff( 43)*x11        *x42    
     8  +coeff( 44)*x11*x21    *x41    
      y_g2_q3en=y_g2_q3en
  
     9  +coeff( 45)    *x23    *x41    
     1  +coeff( 46)    *x22        *x52
     2  +coeff( 47)    *x22    *x43    
     3  +coeff( 48)    *x23    *x42    
     4  +coeff( 49)    *x22    *x42*x51
     5  +coeff( 50)*x11*x22    *x41    
     6  +coeff( 51)    *x24        *x51
c
      return
      end
      function p_g2_q3en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 45)
      data ncoeff/ 44/
      data avdat/  0.6855232E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.64845412E-03, 0.15905603E-02,-0.17579088E-01,-0.14800041E-02,
     + -0.69582285E-02, 0.84107078E-03,-0.12386935E-02, 0.30802118E-03,
     +  0.94645349E-02, 0.15626206E-02,-0.94537463E-04, 0.41431349E-06,
     + -0.35445546E-02,-0.24681434E-02,-0.14552261E-04, 0.83084259E-03,
     + -0.73241508E-05,-0.30458363E-03, 0.12982435E-05,-0.19652836E-02,
     +  0.14815801E-02, 0.37451697E-03,-0.24560184E-03, 0.11250444E-02,
     + -0.65449585E-03,-0.66294946E-03, 0.10817681E-02,-0.22290983E-03,
     + -0.10168316E-02,-0.10712700E-03, 0.10977562E-03, 0.58519745E-04,
     + -0.19221539E-03,-0.98588163E-04, 0.12841013E-03,-0.79523576E-04,
     +  0.97470729E-04, 0.17829399E-03, 0.72903560E-04, 0.27840625E-03,
     +  0.27213609E-03,-0.16756743E-03,-0.43087182E-03, 0.31512923E-03,
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
      p_g2_q3en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)        *x31    *x51
      p_g2_q3en=p_g2_q3en
  
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)        *x33        
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)    *x24            
     8  +coeff( 17)    *x22*x33        
      p_g2_q3en=p_g2_q3en
  
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)        *x33    *x52
     2  +coeff( 20)        *x31        
     3  +coeff( 21)    *x23            
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x21*x31        
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)            *x43    
     8  +coeff( 26)    *x21    *x41*x51
      p_g2_q3en=p_g2_q3en
  
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)*x11*x21        *x51
     2  +coeff( 29)    *x22    *x41*x51
     3  +coeff( 30)*x11                
     4  +coeff( 31)        *x31*x41    
     5  +coeff( 32)*x11            *x51
     6  +coeff( 33)            *x41*x52
     7  +coeff( 34)*x11*x21            
     8  +coeff( 35)    *x21*x31*x41    
      p_g2_q3en=p_g2_q3en
  
     9  +coeff( 36)                *x53
     1  +coeff( 37)*x11*x22            
     2  +coeff( 38)*x11*x21    *x41    
     3  +coeff( 39)*x11        *x42    
     4  +coeff( 40)    *x21    *x43    
     5  +coeff( 41)    *x23        *x51
     6  +coeff( 42)    *x21    *x41*x52
     7  +coeff( 43)    *x24    *x41    
     8  +coeff( 44)    *x24        *x51
c
      return
      end
      function l_g2_q3en(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/ -0.5018328E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44915E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.45277937E-02,-0.22816305E+00,-0.28644461E-01, 0.13980801E-01,
     + -0.20141087E-01, 0.52533623E-01,-0.18050563E-01, 0.37397270E-02,
     +  0.46357475E-02,-0.82039358E-02,-0.40952875E-02,-0.21622083E-02,
     +  0.26010384E-02, 0.11819396E-01,-0.19623082E-01,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_g2_q3en=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)            *x41    
      l_g2_q3en=l_g2_q3en
  
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x23            
     2  +coeff( 11)            *x42    
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)*x11        *x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x23    *x41    
c
      return
      end
      function x_g2_q3ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.8109767E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44149E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.75525539E-02, 0.51040728E-01,-0.14516902E-03,-0.58925947E-04,
     +  0.31067163E+00,-0.23548942E-01, 0.28026540E-01,-0.19422010E-01,
     + -0.10606864E-01, 0.10769982E-01,-0.34792845E-02, 0.42757285E-02,
     + -0.23348131E-02,-0.68500121E-02,-0.47516814E-02, 0.11639501E-01,
     + -0.11165810E-02, 0.22890647E-02, 0.21663390E-02,-0.17720745E-02,
     +  0.95408480E-03, 0.46248181E-03,-0.10718156E-02,-0.13814511E-02,
     +  0.12538406E-02,-0.44152254E-03, 0.23118812E-02,-0.24666506E-03,
     + -0.36069608E-03,-0.22233720E-03, 0.73401199E-03,-0.43787225E-02,
     + -0.34830297E-03,-0.47006336E-03, 0.18042398E-02, 0.10712836E-02,
     + -0.36308466E-03,-0.31567959E-03,-0.17724789E-03, 0.85331971E-03,
     + -0.67363674E-03, 0.46861163E-03, 0.32854351E-03, 0.36965925E-03,
     + -0.24775465E-04,-0.11300758E-03,-0.20941076E-03, 0.53682399E-03,
     + -0.21230362E-03,-0.21223288E-02, 0.51828602E-03,-0.10893322E-02,
     + -0.53010939E-03,-0.10046158E-03,-0.12940752E-03,-0.14358305E-03,
     +  0.54199155E-03, 0.28062833E-03, 0.19943836E-03, 0.16077828E-02,
     +  0.50958950E-03,-0.33415965E-03,-0.22202205E-03,-0.24651206E-03,
     + -0.70532702E-03,-0.91186193E-04, 0.42862562E-03, 0.89230474E-04,
     +  0.22574503E-03,-0.46299730E-04,-0.17619648E-03,-0.51877869E-03,
     + -0.26311912E-03,-0.30060558E-03,-0.47189929E-04,-0.24509465E-03,
     +  0.65989196E-04,-0.27995944E-03, 0.21353866E-03, 0.39015737E-03,
     + -0.66407584E-03,-0.20007258E-03, 0.55441207E-04, 0.60289465E-04,
     +  0.66242843E-04, 0.26118485E-04, 0.17921155E-03, 0.59786486E-04,
     + -0.96427095E-04, 0.57195688E-04,
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
      x_g2_q3ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_g2_q3ex=x_g2_q3ex
  
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x22            
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)            *x41*x51
      x_g2_q3ex=x_g2_q3ex
  
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)                *x53
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)    *x22    *x41    
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)    *x24            
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)*x11            *x51
      x_g2_q3ex=x_g2_q3ex
  
     9  +coeff( 27)*x11*x22    *x41    
     1  +coeff( 28)        *x31*x41    
     2  +coeff( 29)    *x21*x31    *x51
     3  +coeff( 30)    *x21        *x52
     4  +coeff( 31)    *x23*x31        
     5  +coeff( 32)    *x23    *x42    
     6  +coeff( 33)            *x43    
     7  +coeff( 34)            *x42*x51
     8  +coeff( 35)    *x21    *x43    
      x_g2_q3ex=x_g2_q3ex
  
     9  +coeff( 36)    *x21    *x42*x51
     1  +coeff( 37)    *x22        *x52
     2  +coeff( 38)*x11*x21            
     3  +coeff( 39)*x11    *x31        
     4  +coeff( 40)    *x24    *x41    
     5  +coeff( 41)    *x24        *x51
     6  +coeff( 42)    *x23    *x41*x51
     7  +coeff( 43)*x11*x21    *x41    
     8  +coeff( 44)*x11*x22        *x51
      x_g2_q3ex=x_g2_q3ex
  
     9  +coeff( 45)        *x31    *x51
     1  +coeff( 46)        *x31*x41*x51
     2  +coeff( 47)            *x41*x52
     3  +coeff( 48)            *x42*x52
     4  +coeff( 49)                *x54
     5  +coeff( 50)    *x21    *x44    
     6  +coeff( 51)    *x22    *x42*x51
     7  +coeff( 52)    *x21    *x43*x51
     8  +coeff( 53)    *x23        *x52
      x_g2_q3ex=x_g2_q3ex
  
     9  +coeff( 54)    *x21*x32    *x52
     1  +coeff( 55)*x11        *x42    
     2  +coeff( 56)*x11        *x41*x51
     3  +coeff( 57)    *x23    *x41*x52
     4  +coeff( 58)*x11*x23            
     5  +coeff( 59)*x11*x22*x31        
     6  +coeff( 60)    *x23    *x44    
     7  +coeff( 61)    *x24    *x41*x52
     8  +coeff( 62)    *x23    *x42*x52
      x_g2_q3ex=x_g2_q3ex
  
     9  +coeff( 63)*x11*x24            
     1  +coeff( 64)*x11*x23    *x41    
     2  +coeff( 65)*x11*x24    *x41    
     3  +coeff( 66)        *x31*x42    
     4  +coeff( 67)    *x21    *x41*x52
     5  +coeff( 68)        *x31*x41*x52
     6  +coeff( 69)            *x41*x53
     7  +coeff( 70)    *x21*x34        
     8  +coeff( 71)    *x23*x31*x41    
      x_g2_q3ex=x_g2_q3ex
  
     9  +coeff( 72)    *x21    *x42*x52
     1  +coeff( 73)            *x42*x53
     2  +coeff( 74)    *x22        *x54
     3  +coeff( 75)*x11    *x33*x41    
     4  +coeff( 76)    *x24*x32*x42    
     5  +coeff( 77)*x11*x23        *x51
     6  +coeff( 78)*x11*x23    *x41*x51
     7  +coeff( 79)*x11*x21*x33*x41*x51
     8  +coeff( 80)*x11*x23    *x42*x51
      x_g2_q3ex=x_g2_q3ex
  
     9  +coeff( 81)    *x24*x34*x42*x51
     1  +coeff( 82)*x11*x22*x33    *x52
     2  +coeff( 83)    *x22    *x42    
     3  +coeff( 84)    *x21*x31*x41*x51
     4  +coeff( 85)    *x21*x31    *x52
     5  +coeff( 86)    *x24*x31        
     6  +coeff( 87)    *x22    *x43    
     7  +coeff( 88)        *x32*x42*x51
     8  +coeff( 89)    *x21*x31*x41*x52
      x_g2_q3ex=x_g2_q3ex
  
     9  +coeff( 90)        *x31*x42*x52
c
      return
      end
      function t_g2_q3ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 44)
      data ncoeff/ 43/
      data avdat/  0.3036403E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44149E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23791976E-02,-0.14279668E-01,-0.14572630E-05, 0.41488218E-03,
     +  0.11189245E+00, 0.57157180E-02,-0.96265534E-02,-0.17212562E-02,
     +  0.17815313E-02,-0.76301297E-03,-0.11170853E-02, 0.99182501E-03,
     +  0.18561390E-03, 0.40193068E-03,-0.14127358E-02,-0.12777594E-02,
     + -0.63287560E-03,-0.26459608E-03, 0.40464271E-04, 0.49081963E-03,
     + -0.50958304E-03, 0.10190512E-02,-0.60772603E-04,-0.74470518E-04,
     + -0.15035339E-03,-0.91818802E-04,-0.14639217E-03,-0.21259929E-03,
     + -0.64116932E-04,-0.24343246E-03, 0.41144979E-03,-0.83245075E-04,
     + -0.75185664E-04, 0.44417466E-04, 0.32285854E-03, 0.52555195E-04,
     +  0.19139348E-03, 0.10105217E-03, 0.25233932E-03,-0.43600606E-03,
     + -0.31116986E-03, 0.37370759E-03,-0.32849197E-03,
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
      t_g2_q3ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      t_g2_q3ex=t_g2_q3ex
  
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)                *x53
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x21        *x52
      t_g2_q3ex=t_g2_q3ex
  
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x24            
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x21*x31        
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)    *x21*x31*x41    
     8  +coeff( 26)*x11            *x51
      t_g2_q3ex=t_g2_q3ex
  
     9  +coeff( 27)            *x42*x51
     1  +coeff( 28)            *x41*x52
     2  +coeff( 29)        *x33*x41    
     3  +coeff( 30)    *x22        *x52
     4  +coeff( 31)            *x42*x52
     5  +coeff( 32)    *x21*x31    *x51
     6  +coeff( 33)        *x31*x41*x51
     7  +coeff( 34)*x11        *x42    
     8  +coeff( 35)    *x21    *x43    
      t_g2_q3ex=t_g2_q3ex
  
     9  +coeff( 36)    *x23        *x51
     1  +coeff( 37)            *x41*x53
     2  +coeff( 38)*x11*x23            
     3  +coeff( 39)    *x24    *x41    
     4  +coeff( 40)    *x23    *x42    
     5  +coeff( 41)    *x24        *x51
     6  +coeff( 42)    *x23    *x41*x51
     7  +coeff( 43)            *x42*x53
c
      return
      end
      function y_g2_q3ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 44)
      data ncoeff/ 43/
      data avdat/  0.2150922E-01/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44149E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23963295E-02,-0.66321655E-02,-0.49519196E-01, 0.50223442E-02,
     + -0.23948189E-02, 0.36090300E-02,-0.25059134E-01, 0.97451342E-03,
     +  0.33925001E-01,-0.39958153E-02,-0.42790123E-02, 0.50453721E-02,
     +  0.38000096E-02, 0.13038784E-02,-0.13218240E-01,-0.35215127E-02,
     +  0.53409780E-02,-0.91970162E-02,-0.33800085E-02, 0.27481679E-03,
     + -0.73696743E-03,-0.74952282E-03, 0.37140984E-02, 0.37885937E-02,
     + -0.82936644E-03, 0.50204428E-03,-0.30216109E-03,-0.28047876E-02,
     + -0.45552416E-03, 0.10500326E-02, 0.23270314E-03,-0.66484150E-03,
     +  0.43191595E-03, 0.12284821E-02, 0.38077738E-03, 0.15511081E-02,
     +  0.14617748E-02,-0.90617046E-03, 0.35189110E-03,-0.26685116E-03,
     + -0.21130822E-02,-0.17810366E-02, 0.10535924E-02,
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
      x31 = x3
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
      y_g2_q3ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31    *x51
      y_g2_q3ex=y_g2_q3ex
  
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x23            
      y_g2_q3ex=y_g2_q3ex
  
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x22    *x41*x51
     2  +coeff( 20)            *x45    
     3  +coeff( 21)    *x21*x31        
     4  +coeff( 22)*x11*x21            
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)    *x24            
     7  +coeff( 25)*x11*x21        *x51
     8  +coeff( 26)        *x31*x41    
      y_g2_q3ex=y_g2_q3ex
  
     9  +coeff( 27)*x11                
     1  +coeff( 28)            *x43    
     2  +coeff( 29)    *x22*x31        
     3  +coeff( 30)            *x41*x52
     4  +coeff( 31)*x11            *x51
     5  +coeff( 32)    *x21        *x52
     6  +coeff( 33)*x11*x21    *x41    
     7  +coeff( 34)    *x23        *x51
     8  +coeff( 35)    *x21*x31*x41    
      y_g2_q3ex=y_g2_q3ex
  
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)    *x23    *x41    
     2  +coeff( 38)            *x41*x53
     3  +coeff( 39)*x11*x22            
     4  +coeff( 40)        *x31*x44    
     5  +coeff( 41)    *x23    *x42    
     6  +coeff( 42)    *x24    *x41    
     7  +coeff( 43)    *x24        *x51
c
      return
      end
      function p_g2_q3ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 60)
      data ncoeff/ 59/
      data avdat/ -0.8157184E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44149E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.92995056E-03,-0.17696018E-02, 0.29566993E-02, 0.12796267E-01,
     +  0.68586558E-03, 0.49838410E-02, 0.11819820E-01,-0.24157062E-02,
     +  0.21647224E-02,-0.10099518E-02,-0.20418661E-01,-0.31601414E-02,
     +  0.90933550E-03,-0.28270732E-02, 0.73778140E-02,-0.17064959E-02,
     +  0.42570401E-02, 0.31495602E-02,-0.33625634E-03,-0.62257051E-03,
     +  0.13619812E-02, 0.11665968E-02, 0.74680592E-03,-0.23948080E-02,
     + -0.22547031E-02, 0.29247592E-03, 0.42067471E-03, 0.37615624E-03,
     +  0.36595028E-03, 0.36581210E-03,-0.11328419E-02, 0.20035241E-03,
     + -0.78977711E-04, 0.22202150E-03,-0.24574288E-03,-0.69551832E-04,
     +  0.12444268E-02,-0.16056813E-03, 0.16435231E-03,-0.65842672E-03,
     + -0.13791365E-03,-0.14837783E-03,-0.66360657E-03,-0.38063648E-03,
     + -0.60544442E-03, 0.20305582E-03,-0.25222576E-03,-0.45530294E-03,
     + -0.33259718E-03, 0.19687046E-03, 0.89942553E-03, 0.71942090E-03,
     + -0.34128258E-03, 0.12900589E-03, 0.74431108E-03,-0.18928938E-03,
     +  0.61738089E-03,-0.41050045E-03,-0.27307015E-03,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      p_g2_q3ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)            *x42    
      p_g2_q3ex=p_g2_q3ex
  
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x22        *x51
      p_g2_q3ex=p_g2_q3ex
  
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)            *x43    
     4  +coeff( 22)    *x21    *x41*x51
     5  +coeff( 23)                *x53
     6  +coeff( 24)    *x24            
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)    *x21*x31        
      p_g2_q3ex=p_g2_q3ex
  
     9  +coeff( 27)    *x22*x31        
     1  +coeff( 28)            *x42*x51
     2  +coeff( 29)*x11*x21        *x51
     3  +coeff( 30)    *x22    *x41*x51
     4  +coeff( 31)            *x41*x53
     5  +coeff( 32)*x11                
     6  +coeff( 33)*x11            *x51
     7  +coeff( 34)        *x31    *x52
     8  +coeff( 35)*x11*x22            
      p_g2_q3ex=p_g2_q3ex
  
     9  +coeff( 36)*x11*x21    *x41    
     1  +coeff( 37)    *x24    *x41    
     2  +coeff( 38)    *x21*x31*x41    
     3  +coeff( 39)        *x31*x42    
     4  +coeff( 40)    *x23    *x41    
     5  +coeff( 41)    *x22*x31*x41    
     6  +coeff( 42)*x11        *x42    
     7  +coeff( 43)    *x21    *x43    
     8  +coeff( 44)            *x43*x51
      p_g2_q3ex=p_g2_q3ex
  
     9  +coeff( 45)    *x22        *x52
     1  +coeff( 46)    *x21    *x41*x52
     2  +coeff( 47)            *x42*x52
     3  +coeff( 48)                *x54
     4  +coeff( 49)*x11*x23            
     5  +coeff( 50)*x11*x22    *x41    
     6  +coeff( 51)    *x23    *x42    
     7  +coeff( 52)    *x22    *x43    
     8  +coeff( 53)    *x24        *x51
      p_g2_q3ex=p_g2_q3ex
  
     9  +coeff( 54)*x11*x21    *x41*x51
     1  +coeff( 55)    *x22    *x42*x51
     2  +coeff( 56)    *x23        *x52
     3  +coeff( 57)            *x41*x54
     4  +coeff( 58)    *x23    *x42*x51
     5  +coeff( 59)    *x23        *x53
c
      return
      end
      function l_g2_q3ex(x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 18)
      data ncoeff/ 17/
      data avdat/ -0.5632350E-02/
      data xmin/
     1 -0.49974E-02,-0.47995E-01,-0.29970E-02,-0.26782E-01,-0.44149E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49966E-02, 0.47890E-01, 0.29992E-02, 0.14503E-01, 0.44952E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.51977853E-02,-0.22862345E+00,-0.27678281E-01, 0.13941159E-01,
     + -0.23077806E-01, 0.52787311E-01,-0.13060484E-01, 0.31855190E-02,
     +  0.47115963E-02,-0.80636283E-02,-0.82897237E-02,-0.49476214E-02,
     + -0.25135593E-02, 0.24816373E-02, 0.42414987E-02, 0.11573529E-01,
     + -0.19374723E-01,
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
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_g2_q3ex=avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)            *x41    
      l_g2_q3ex=l_g2_q3ex
  
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x23            
     3  +coeff( 12)            *x42    
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x23    *x41    
c
      return
      end
