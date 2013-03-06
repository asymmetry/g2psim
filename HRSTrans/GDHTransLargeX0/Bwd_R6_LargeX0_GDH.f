C Right HRS with septum reverse functions. 
C Generated with larger x0 range for Vince (-1.6cm < x0 < 1.6 cm)
C Usage is the same as before.
C                                                   -JJL 2/6/07
      function r6_largex0_txfit       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.2441093E-02/
      data xmin/
     1 -0.52553E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.46408E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.30284817E-02, 0.85081957E-01,-0.83913520E-03,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
c
c                  function
c
      r6_largex0_txfit       =avdat
     1  +coeff(  1)    
     2  +coeff(  2)*x11
     3  +coeff(  3)*x12
c
      return
      end
      function r6_largex0_delta       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1584944E-02/
      data xmin/
     1 -0.67583E+00,-0.25849E-01,-0.39793E-01,-0.36128E-01,-0.15997E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.61323E+00, 0.22501E-01, 0.40101E-01, 0.32451E-01, 0.15975E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.39418694E-02, 0.48422143E-01, 0.38988586E-02, 0.41510155E-02,
     +  0.26872354E-02,-0.16418886E-02, 0.20470600E-03, 0.11169014E-02,
     + -0.79277536E-03, 0.59718965E-03, 0.33109690E-03, 0.14010859E-03,
     +  0.63675630E-03, 0.22563360E-03,-0.28825470E-02, 0.39774007E-02,
     + -0.19968199E-02, 0.71180123E-03,-0.21079915E-03, 0.11509127E-03,
     +  0.79375721E-04,-0.22415987E-02,-0.31205101E-03, 0.85190943E-04,
     +  0.10034697E-02,-0.82441598E-04, 0.72684888E-04, 0.25267771E-03,
     +  0.60045865E-03, 0.19291657E-03,-0.21618135E-03, 0.27629867E-03,
     + -0.12091933E-02,-0.14038467E-03, 0.31200091E-02,-0.77243960E-02,
     + -0.13367868E-02, 0.79845255E-02,-0.30970140E-02, 0.49037044E-03,
     + -0.89551730E-03, 0.40157171E-03,-0.18245415E-04,-0.13404409E-02,
     +  0.22180653E-02,-0.71895763E-03, 0.35983376E-03, 0.15900227E-04,
     +  0.17612035E-03,-0.23222662E-03,-0.90202047E-04,-0.63233339E-04,
     +  0.16500799E-03,-0.38627510E-04, 0.23533944E-02, 0.67616842E-03,
     + -0.39440705E-03, 0.45787846E-03, 0.13267828E-03,-0.16090751E-03,
     +  0.40421681E-04, 0.30342115E-02,-0.41376832E-02,-0.86545460E-02,
     +  0.88005383E-02,-0.34609304E-02,-0.64105709E-03, 0.36864675E-03,
     +  0.30903955E-03, 0.42533630E-03,-0.30864446E-03,-0.58864674E-03,
     +  0.36419203E-03, 0.11500115E-03, 0.25219740E-02, 0.92386355E-03,
     + -0.13475715E-03, 0.15916453E-02, 0.18928262E-02,-0.67517598E-03,
     + -0.38137094E-02,-0.12300305E-03, 0.27383477E-03,-0.13242220E-03,
     +  0.16005991E-03, 0.16982279E-03,-0.50015259E-03, 0.43285138E-03,
     + -0.33144854E-04,-0.45986733E-04,
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
      x13 = x12*x1
      x14 = x13*x1
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
      r6_largex0_delta       =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)*x12                
     4  +coeff(  4)*x11*x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21            
     8  +coeff(  8)        *x32        
      r6_largex0_delta       =r6_largex0_delta       
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)*x11*x21*x31        
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x21    *x42    
      r6_largex0_delta       =r6_largex0_delta       
     9  +coeff( 18)    *x24            
     1  +coeff( 19)    *x21*x31        
     2  +coeff( 20)*x12*x21            
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)    *x22*x32        
     5  +coeff( 23)    *x24*x31        
     6  +coeff( 24)*x12*x22    *x41    
     7  +coeff( 25)*x12        *x44    
     8  +coeff( 26)*x11        *x41    
      r6_largex0_delta       =r6_largex0_delta       
     9  +coeff( 27)    *x22*x31        
     1  +coeff( 28)    *x22    *x41    
     2  +coeff( 29)*x11    *x31*x41    
     3  +coeff( 30)*x11        *x42    
     4  +coeff( 31)        *x31    *x51
     5  +coeff( 32)            *x41*x51
     6  +coeff( 33)*x11*x21*x32        
     7  +coeff( 34)    *x23    *x41    
     8  +coeff( 35)    *x22*x31*x41    
      r6_largex0_delta       =r6_largex0_delta       
     9  +coeff( 36)    *x21*x32*x41    
     1  +coeff( 37)    *x22    *x42    
     2  +coeff( 38)    *x21*x31*x42    
     3  +coeff( 39)    *x21    *x43    
     4  +coeff( 40)        *x32    *x51
     5  +coeff( 41)        *x31*x41*x51
     6  +coeff( 42)            *x42*x51
     7  +coeff( 43)*x12*x23            
     8  +coeff( 44)*x11*x22*x32        
      r6_largex0_delta       =r6_largex0_delta       
     9  +coeff( 45)*x11*x21*x33        
     1  +coeff( 46)*x11        *x44    
     2  +coeff( 47)    *x21*x34*x41    
     3  +coeff( 48)            *x41    
     4  +coeff( 49)*x13                
     5  +coeff( 50)    *x23            
     6  +coeff( 51)*x12        *x41    
     7  +coeff( 52)*x11*x21    *x41    
     8  +coeff( 53)*x12*x22            
      r6_largex0_delta       =r6_largex0_delta       
     9  +coeff( 54)    *x23*x31        
     1  +coeff( 55)    *x21*x33        
     2  +coeff( 56)*x11*x21*x31*x41    
     3  +coeff( 57)*x12        *x42    
     4  +coeff( 58)*x11        *x43    
     5  +coeff( 59)    *x22        *x51
     6  +coeff( 60)*x11    *x31    *x51
     7  +coeff( 61)*x11        *x41*x51
     8  +coeff( 62)    *x22*x33        
      r6_largex0_delta       =r6_largex0_delta       
     9  +coeff( 63)*x11*x21*x32*x41    
     1  +coeff( 64)    *x22*x32*x41    
     2  +coeff( 65)    *x22*x31*x42    
     3  +coeff( 66)    *x22    *x43    
     4  +coeff( 67)        *x31*x44    
     5  +coeff( 68)    *x22*x31    *x51
     6  +coeff( 69)*x11    *x32    *x51
     7  +coeff( 70)    *x21*x32    *x51
     8  +coeff( 71)    *x22    *x41*x51
      r6_largex0_delta       =r6_largex0_delta       
     9  +coeff( 72)    *x21*x31*x41*x51
     1  +coeff( 73)    *x21    *x42*x51
     2  +coeff( 74)*x13*x22*x31        
     3  +coeff( 75)*x11*x22*x33        
     4  +coeff( 76)*x11*x22    *x43    
     5  +coeff( 77)    *x21    *x41*x52
     6  +coeff( 78)*x11*x21*x31*x44    
     7  +coeff( 79)*x14*x21*x33        
     8  +coeff( 80)*x14*x21*x32*x41    
      r6_largex0_delta       =r6_largex0_delta       
     9  +coeff( 81)*x13*x22*x32*x41    
     1  +coeff( 82)*x11    *x31        
     2  +coeff( 83)            *x43    
     3  +coeff( 84)*x11*x23            
     4  +coeff( 85)*x13    *x31        
     5  +coeff( 86)*x11*x21    *x42    
     6  +coeff( 87)        *x31*x43    
     7  +coeff( 88)            *x44    
     8  +coeff( 89)*x11*x21        *x51
      r6_largex0_delta       =r6_largex0_delta       
     9  +coeff( 90)    *x21*x31    *x51
c
      return
      end
      function r6_largex0_theta       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.7813832E-03/
      data xmin/
     1 -0.67583E+00,-0.25849E-01,-0.39793E-01,-0.36128E-01,-0.15997E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.61323E+00, 0.22501E-01, 0.40101E-01, 0.32451E-01, 0.15975E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29442890E-02,-0.39248238E-03,-0.53575411E-01, 0.68268436E-02,
     + -0.36589282E-02, 0.17584240E-02,-0.11705750E-02, 0.17643815E-02,
     +  0.14471002E-01, 0.57887179E-02,-0.12876722E-02, 0.28066247E-03,
     + -0.11906703E-01, 0.15193662E-01,-0.66189799E-02,-0.11138616E-02,
     +  0.14003222E-02, 0.51322929E-03, 0.54982636E-03, 0.30202705E-02,
     + -0.87028258E-02, 0.15372867E-01,-0.43759312E-01,-0.57280441E-02,
     +  0.46955019E-01,-0.18120684E-01, 0.34041372E-02,-0.12823894E-02,
     +  0.13383236E-02,-0.83381608E-02,-0.25489336E-03, 0.76658610E-03,
     + -0.10438929E-02, 0.43656095E-02,-0.44719209E-02, 0.38536622E-02,
     +  0.44841850E-02, 0.17598674E-01,-0.49540669E-01, 0.53106762E-01,
     + -0.20073844E-01,-0.83916392E-02, 0.21603245E-02,-0.45879898E-02,
     +  0.22057667E-02,-0.25541467E-04,-0.75669358E-02,-0.82798957E-04,
     + -0.24953883E-03,-0.41395423E-03,-0.41757367E-03,-0.70772064E-02,
     +  0.52821526E-03,-0.17499829E-02, 0.54816650E-02,-0.15556345E-02,
     + -0.25252223E-02,-0.27159199E-02, 0.26025830E-02, 0.98187728E-02,
     + -0.29310831E-03, 0.22110785E-02, 0.18014500E-02,-0.40802178E-02,
     + -0.40235128E-02,-0.22923941E-01, 0.17686784E-01,-0.63165776E-02,
     +  0.19678457E-02,-0.30395530E-01, 0.59628394E-02, 0.18978676E-01,
     +  0.44783915E-03, 0.44472683E-02, 0.21529607E-02,-0.49846079E-02,
     + -0.22923654E-03,-0.31180857E-03, 0.22205631E-02, 0.89652492E-02,
     +  0.62206894E-03,-0.44877934E-02, 0.13611293E-02, 0.13189689E-02,
     + -0.94545609E-03, 0.13613647E-02, 0.89714061E-02, 0.15721451E-01,
     + -0.65173996E-02,-0.14217433E-02,
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
      x13 = x12*x1
      x14 = x13*x1
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
c
c                  function
c
      r6_largex0_theta       =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)*x11*x21            
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)*x11*x21*x31        
      r6_largex0_theta       =r6_largex0_theta       
     9  +coeff(  9)    *x21*x33        
     1  +coeff( 10)                *x51
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)*x11*x22            
      r6_largex0_theta       =r6_largex0_theta       
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x22    *x41    
     2  +coeff( 20)            *x43    
     3  +coeff( 21)    *x22*x32        
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)    *x21*x32*x41    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)    *x21*x31*x42    
     8  +coeff( 26)    *x21    *x43    
      r6_largex0_theta       =r6_largex0_theta       
     9  +coeff( 27)    *x23    *x42    
     1  +coeff( 28)        *x31    *x51
     2  +coeff( 29)            *x41*x51
     3  +coeff( 30)*x13*x21*x32        
     4  +coeff( 31)            *x41    
     5  +coeff( 32)        *x32        
     6  +coeff( 33)        *x31*x41    
     7  +coeff( 34)        *x32*x41    
     8  +coeff( 35)*x11*x21*x32        
      r6_largex0_theta       =r6_largex0_theta       
     9  +coeff( 36)*x11*x21*x31*x41    
     1  +coeff( 37)            *x44    
     2  +coeff( 38)    *x22*x33        
     3  +coeff( 39)    *x22*x32*x41    
     4  +coeff( 40)    *x22*x31*x42    
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)    *x21    *x44    
     7  +coeff( 43)        *x32    *x51
     8  +coeff( 44)        *x31*x41*x51
      r6_largex0_theta       =r6_largex0_theta       
     9  +coeff( 45)            *x42*x51
     1  +coeff( 46)*x12*x21        *x51
     2  +coeff( 47)*x14*x22*x32        
     3  +coeff( 48)        *x31        
     4  +coeff( 49)*x13                
     5  +coeff( 50)*x12        *x41    
     6  +coeff( 51)*x11*x21    *x41    
     7  +coeff( 52)        *x31*x42    
     8  +coeff( 53)*x11*x23            
      r6_largex0_theta       =r6_largex0_theta       
     9  +coeff( 54)    *x23    *x41    
     1  +coeff( 55)*x11        *x43    
     2  +coeff( 56)        *x31*x43    
     3  +coeff( 57)*x11*x23*x31        
     4  +coeff( 58)    *x24*x31        
     5  +coeff( 59)    *x23*x32        
     6  +coeff( 60)*x11*x21*x33        
     7  +coeff( 61)*x14        *x41    
     8  +coeff( 62)*x11*x23    *x41    
      r6_largex0_theta       =r6_largex0_theta       
     9  +coeff( 63)    *x24    *x41    
     1  +coeff( 64)*x11*x22*x31*x41    
     2  +coeff( 65)    *x23*x31*x41    
     3  +coeff( 66)*x11*x21*x32*x41    
     4  +coeff( 67)    *x21*x33*x41    
     5  +coeff( 68)        *x34*x41    
     6  +coeff( 69)*x11*x22    *x42    
     7  +coeff( 70)    *x21*x32*x42    
     8  +coeff( 71)        *x33*x42    
      r6_largex0_theta       =r6_largex0_theta       
     9  +coeff( 72)    *x21*x31*x43    
     1  +coeff( 73)    *x21        *x51
     2  +coeff( 74)*x11*x22*x33        
     3  +coeff( 75)*x14    *x31*x41    
     4  +coeff( 76)*x13        *x43    
     5  +coeff( 77)*x11*x21        *x51
     6  +coeff( 78)*x11    *x31    *x51
     7  +coeff( 79)*x14*x22*x31        
     8  +coeff( 80)*x13*x21*x33        
      r6_largex0_theta       =r6_largex0_theta       
     9  +coeff( 81)*x14    *x31*x42    
     1  +coeff( 82)*x13        *x44    
     2  +coeff( 83)    *x22*x31    *x51
     3  +coeff( 84)*x11    *x32    *x51
     4  +coeff( 85)    *x22    *x41*x51
     5  +coeff( 86)        *x32*x41*x51
     6  +coeff( 87)*x13*x23*x32        
     7  +coeff( 88)*x13*x21*x32*x42    
     8  +coeff( 89)*x12*x22    *x44    
      r6_largex0_theta       =r6_largex0_theta       
     9  +coeff( 90)*x11    *x33    *x51
c
      return
      end
      function r6_largex0_phi         (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.3705922E-03/
      data xmin/
     1 -0.67583E+00,-0.25849E-01,-0.39793E-01,-0.36128E-01,-0.15997E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.61323E+00, 0.22501E-01, 0.40101E-01, 0.32451E-01, 0.15975E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.76069019E-03,-0.51631555E-02,-0.46790289E-03,-0.29668087E-01,
     + -0.25104263E-02,-0.64609991E-03, 0.19060638E-02, 0.99874260E-02,
     + -0.12733574E-01, 0.27220577E-02, 0.11339004E-02, 0.17595605E-03,
     + -0.68490757E-02, 0.52845310E-02, 0.87090116E-02,-0.32428117E-04,
     +  0.16961724E-01, 0.17685500E-02,-0.29332731E-01,-0.25652202E-02,
     +  0.17493308E-02,-0.35082102E-02, 0.37291129E-02,-0.39431760E-02,
     +  0.29883678E-02, 0.68213638E-04,-0.77262992E-03, 0.23515155E-02,
     +  0.31587521E-02,-0.28231451E-02,-0.64326200E-03, 0.86394614E-02,
     + -0.17171795E-01, 0.37094292E-02,-0.46626772E-02, 0.14536219E-01,
     + -0.29045525E-02,-0.90061678E-02,-0.13358745E-01, 0.14179836E-01,
     + -0.27194943E-02, 0.40093567E-02, 0.12948631E-01, 0.38812710E-02,
     + -0.29681096E-01, 0.24463027E-02,-0.12765319E-01, 0.10111972E-01,
     +  0.15333164E-01,-0.10848352E-01, 0.84593808E-02,-0.33329104E-02,
     +  0.24820291E-02,-0.11248270E-01,-0.18012095E-02,-0.35706304E-01,
     +  0.17164035E-01,-0.87075685E-04,-0.98603331E-02,-0.16321668E-01,
     + -0.65745194E-02,-0.18361930E-01, 0.17629107E-01, 0.15877314E-01,
     +  0.29246191E-01,-0.16727062E-01,-0.13947943E-02, 0.33474367E-01,
     + -0.80247801E-02,-0.91881603E-02, 0.41195785E-03,-0.84544308E-02,
     +  0.63194670E-02, 0.15258566E-01,-0.72519188E-02, 0.10354561E-01,
     +  0.35561200E-01,-0.18219687E-01,-0.27683924E-02,-0.45940671E-01,
     +  0.55880673E-01, 0.12045313E-01, 0.10797638E-01, 0.71332366E-02,
     + -0.18926198E-01,-0.48727289E-03, 0.46758596E-02, 0.73359185E-02,
     + -0.15035718E-01,-0.97713741E-02,
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
      x13 = x12*x1
      x14 = x13*x1
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
c
c                  function
c
      r6_largex0_phi         =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)        *x31        
     5  +coeff(  5)            *x41    
     6  +coeff(  6)*x12                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)*x11    *x31        
      r6_largex0_phi         =r6_largex0_phi         
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)*x12        *x41    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)*x11    *x31*x41    
     8  +coeff( 17)            *x43    
      r6_largex0_phi         =r6_largex0_phi         
     9  +coeff( 18)*x12    *x32        
     1  +coeff( 19)*x11    *x31*x42    
     2  +coeff( 20)    *x21        *x51
     3  +coeff( 21)    *x22*x34        
     4  +coeff( 22)*x11*x23    *x42    
     5  +coeff( 23)    *x21*x31    *x51
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)    *x21*x31        
     8  +coeff( 26)        *x32        
      r6_largex0_phi         =r6_largex0_phi         
     9  +coeff( 27)*x11*x22            
     1  +coeff( 28)*x11    *x32        
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)    *x21*x31*x41    
     5  +coeff( 32)*x11        *x42    
     6  +coeff( 33)        *x31*x42    
     7  +coeff( 34)*x12*x21*x31        
     8  +coeff( 35)*x11*x21*x32        
      r6_largex0_phi         =r6_largex0_phi         
     9  +coeff( 36)    *x22*x32        
     1  +coeff( 37)    *x21*x33        
     2  +coeff( 38)*x11*x22    *x41    
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)*x11    *x32*x41    
     5  +coeff( 41)*x14    *x31        
     6  +coeff( 42)*x14        *x41    
     7  +coeff( 43)    *x23    *x42    
     8  +coeff( 44)*x11*x24*x31        
      r6_largex0_phi         =r6_largex0_phi         
     9  +coeff( 45)    *x22*x33*x41    
     1  +coeff( 46)    *x23        *x51
     2  +coeff( 47)    *x23*x32    *x51
     3  +coeff( 48)*x12*x21*x31*x41*x51
     4  +coeff( 49)    *x23*x31*x41*x51
     5  +coeff( 50)*x12*x21    *x42*x51
     6  +coeff( 51)        *x32*x41    
     7  +coeff( 52)*x11    *x33        
     8  +coeff( 53)*x12*x21    *x41    
      r6_largex0_phi         =r6_largex0_phi         
     9  +coeff( 54)*x12    *x31*x41    
     1  +coeff( 55)    *x21*x31*x42    
     2  +coeff( 56)        *x31*x43    
     3  +coeff( 57)            *x44    
     4  +coeff( 58)                *x51
     5  +coeff( 59)*x12    *x33        
     6  +coeff( 60)    *x22*x33        
     7  +coeff( 61)*x13    *x31*x41    
     8  +coeff( 62)    *x23*x31*x41    
      r6_largex0_phi         =r6_largex0_phi         
     9  +coeff( 63)*x12    *x32*x41    
     1  +coeff( 64)*x11*x21*x32*x41    
     2  +coeff( 65)    *x22*x32*x41    
     3  +coeff( 66)*x11    *x33*x41    
     4  +coeff( 67)    *x21*x33*x41    
     5  +coeff( 68)*x11    *x32*x42    
     6  +coeff( 69)*x11*x21    *x43    
     7  +coeff( 70)*x11    *x31*x43    
     8  +coeff( 71)        *x31    *x51
      r6_largex0_phi         =r6_largex0_phi         
     9  +coeff( 72)*x14*x21*x31        
     1  +coeff( 73)*x14    *x32        
     2  +coeff( 74)*x11*x23*x32        
     3  +coeff( 75)    *x23*x33        
     4  +coeff( 76)*x11*x21*x34        
     5  +coeff( 77)    *x23*x32*x41    
     6  +coeff( 78)*x11*x21*x33*x41    
     7  +coeff( 79)    *x21*x34*x41    
     8  +coeff( 80)    *x23*x31*x42    
      r6_largex0_phi         =r6_largex0_phi         
     9  +coeff( 81)    *x22*x32*x42    
     1  +coeff( 82)        *x34*x42    
     2  +coeff( 83)*x11*x22    *x43    
     3  +coeff( 84)    *x23    *x43    
     4  +coeff( 85)    *x22*x31*x43    
     5  +coeff( 86)    *x22        *x51
     6  +coeff( 87)*x14*x21*x32        
     7  +coeff( 88)*x13*x22*x32        
     8  +coeff( 89)*x11*x23*x33        
      r6_largex0_phi         =r6_largex0_phi         
     9  +coeff( 90)*x11*x24    *x42    
c
      return
      end
      function r6_largex0_y00         (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.3046953E-03/
      data xmin/
     1 -0.67583E+00,-0.25849E-01,-0.39793E-01,-0.36128E-01,-0.15997E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.61323E+00, 0.22501E-01, 0.40101E-01, 0.32451E-01, 0.15975E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.17881116E-02,-0.11179538E-03,-0.53051598E-02, 0.50579596E-01,
     +  0.19572682E-02,-0.38595523E-02,-0.32353792E-01,-0.28470013E-01,
     + -0.17547189E-02, 0.13679782E-01, 0.61742850E-02, 0.47940826E-02,
     +  0.31043391E-03,-0.46841131E-03, 0.17417255E-02, 0.15455346E-02,
     +  0.91192726E-03, 0.25456723E-01,-0.42737382E-02, 0.11291203E-01,
     + -0.13873124E-01,-0.81193727E-02,-0.13591440E-01, 0.12249426E-01,
     +  0.96788993E-02,-0.26435347E-01,-0.23252158E-02, 0.78431086E-03,
     +  0.56084967E-03,-0.24873433E-02,-0.92601515E-02,-0.15897330E-01,
     + -0.21008400E-01, 0.32961417E-01,-0.10062514E-02, 0.71843946E-03,
     +  0.93204323E-02,-0.13470071E-02,-0.31179123E-03, 0.65567857E-02,
     +  0.17166083E-02, 0.69060316E-02, 0.28991058E-01, 0.36894374E-01,
     +  0.28510159E-02,-0.13707466E-02, 0.49681813E-02, 0.11784635E-01,
     +  0.13844735E-01,-0.33836376E-02,-0.45498811E-01, 0.65573208E-01,
     + -0.75437434E-01,-0.60916143E-02, 0.61282851E-02, 0.37983733E-04,
     + -0.28803141E-02, 0.98821977E-02,-0.21561824E-01, 0.12054563E-01,
     +  0.76116691E-02,-0.11577561E-01, 0.55087539E-02, 0.83448794E-02,
     +  0.12648228E-02,-0.17805016E-02,-0.32764949E-01,-0.16904676E-01,
     + -0.48565292E-02, 0.34910452E-01, 0.49152398E-02, 0.90219313E-02,
     +  0.11818107E-01, 0.12982429E-01, 0.25803769E-01,-0.34892663E-01,
     + -0.62211566E-02, 0.36536250E-01,-0.34062527E-01,-0.63390732E-02,
     + -0.17223727E-01, 0.36215037E-01,-0.15733983E-01, 0.17696408E-03,
     + -0.17664217E-02,-0.78673055E-02,-0.10070235E-01,-0.14986739E-01,
     +  0.12438956E-01,-0.15010736E-01,
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
      x13 = x12*x1
      x14 = x13*x1
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
c
c                  function
c
      r6_largex0_y00         =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x12                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11    *x31        
     8  +coeff(  8)    *x21*x31        
      r6_largex0_y00         =r6_largex0_y00         
     9  +coeff(  9)        *x32        
     1  +coeff( 10)*x11        *x41    
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)*x13                
     6  +coeff( 15)*x12*x21            
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)*x11*x21*x31        
      r6_largex0_y00         =r6_largex0_y00         
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)*x11    *x32        
     2  +coeff( 20)        *x33        
     3  +coeff( 21)*x12        *x41    
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x22    *x41    
     6  +coeff( 24)*x11    *x31*x41    
     7  +coeff( 25)    *x21*x31*x41    
     8  +coeff( 26)            *x43    
      r6_largex0_y00         =r6_largex0_y00         
     9  +coeff( 27)*x14                
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)    *x24            
     3  +coeff( 30)*x12    *x32        
     4  +coeff( 31)*x11*x21*x32        
     5  +coeff( 32)    *x22*x32        
     6  +coeff( 33)*x11    *x32*x41    
     7  +coeff( 34)*x11    *x31*x42    
     8  +coeff( 35)*x14*x21            
      r6_largex0_y00         =r6_largex0_y00         
     9  +coeff( 36)*x13*x22            
     1  +coeff( 37)*x14    *x31        
     2  +coeff( 38)*x13    *x32        
     3  +coeff( 39)*x14        *x41    
     4  +coeff( 40)*x13    *x31*x41    
     5  +coeff( 41)*x13        *x42    
     6  +coeff( 42)*x12    *x31*x42    
     7  +coeff( 43)*x11        *x44    
     8  +coeff( 44)        *x31*x44    
      r6_largex0_y00         =r6_largex0_y00         
     9  +coeff( 45)    *x21        *x51
     1  +coeff( 46)        *x31    *x51
     2  +coeff( 47)*x12*x21*x33        
     3  +coeff( 48)*x13*x21*x31*x41    
     4  +coeff( 49)*x11*x21*x33*x41    
     5  +coeff( 50)    *x21*x34*x41    
     6  +coeff( 51)*x11    *x32*x43    
     7  +coeff( 52)        *x33*x43    
     8  +coeff( 53)        *x32*x44    
      r6_largex0_y00         =r6_largex0_y00         
     9  +coeff( 54)    *x21*x31    *x51
     1  +coeff( 55)    *x21    *x41*x51
     2  +coeff( 56)*x11    *x32*x44    
     3  +coeff( 57)    *x23        *x51
     4  +coeff( 58)    *x21*x32    *x51
     5  +coeff( 59)    *x21*x31*x41*x51
     6  +coeff( 60)    *x21    *x42*x51
     7  +coeff( 61)    *x22*x32    *x51
     8  +coeff( 62)    *x21*x31*x42*x51
      r6_largex0_y00         =r6_largex0_y00         
     9  +coeff( 63)        *x31*x44*x51
     1  +coeff( 64)*x11                
     2  +coeff( 65)*x11*x21            
     3  +coeff( 66)    *x21*x32        
     4  +coeff( 67)        *x32*x41    
     5  +coeff( 68)*x11        *x42    
     6  +coeff( 69)    *x21    *x42    
     7  +coeff( 70)        *x31*x42    
     8  +coeff( 71)    *x23*x31        
      r6_largex0_y00         =r6_largex0_y00         
     9  +coeff( 72)*x11    *x33        
     1  +coeff( 73)    *x21*x33        
     2  +coeff( 74)*x11*x22    *x41    
     3  +coeff( 75)    *x22*x31*x41    
     4  +coeff( 76)    *x21*x32*x41    
     5  +coeff( 77)    *x22    *x42    
     6  +coeff( 78)    *x21*x31*x42    
     7  +coeff( 79)        *x32*x42    
     8  +coeff( 80)*x11        *x43    
      r6_largex0_y00         =r6_largex0_y00         
     9  +coeff( 81)    *x21    *x43    
     1  +coeff( 82)        *x31*x43    
     2  +coeff( 83)            *x44    
     3  +coeff( 84)                *x51
     4  +coeff( 85)*x12*x23            
     5  +coeff( 86)*x12*x22*x31        
     6  +coeff( 87)*x11*x22*x32        
     7  +coeff( 88)    *x23*x32        
     8  +coeff( 89)    *x23*x31*x41    
      r6_largex0_y00         =r6_largex0_y00         
     9  +coeff( 90)*x12*x21    *x42    
c
      return
      end
