      function delta_r12p5 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.8675047E-03/
      data xmin/
     1 -0.69647E+00,-0.25366E-01,-0.94322E-01,-0.34681E-01,-0.99237E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.61730E+00, 0.38656E-01, 0.94322E-01, 0.34681E-01, 0.99987E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.57736752E-05, 0.46993420E-01, 0.10767586E-01, 0.15088027E-01,
     +  0.57990588E-02, 0.26363484E-02, 0.44245147E-02, 0.12310279E-02,
     +  0.16267250E-02,-0.17092453E-02,-0.13197486E-02, 0.54015401E-02,
     +  0.73287118E-03,-0.43749347E-03, 0.10190976E-02,-0.30467438E-02,
     + -0.18863851E-02, 0.37169273E-02, 0.44369372E-03, 0.90896647E-03,
     + -0.51580206E-03, 0.16229262E-03, 0.70138974E-03, 0.30491862E-02,
     + -0.70686440E-03, 0.23058206E-02,-0.39973045E-02,
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
      delta_r12p5 =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11*x21            
     6  +coeff(  6)*x12                
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)*x11            *x51
      delta_r12p5 =delta_r12p5 
     9  +coeff(  9)    *x21    *x44*x51
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x21*x32*x42    
     3  +coeff( 12)    *x24*x32*x42*x51
     4  +coeff( 13)    *x22            
     5  +coeff( 14)        *x32        
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)    *x21    *x42    
      delta_r12p5 =delta_r12p5 
     9  +coeff( 18)    *x22    *x42*x51
     1  +coeff( 19)            *x42    
     2  +coeff( 20)*x11    *x32        
     3  +coeff( 21)        *x31*x41*x51
     4  +coeff( 22)    *x22*x32        
     5  +coeff( 23)*x11*x21    *x42    
     6  +coeff( 24)    *x23*x32        
     7  +coeff( 25)    *x22    *x44    
     8  +coeff( 26)    *x23    *x42*x51
      delta_r12p5 =delta_r12p5 
     9  +coeff( 27)*x14*x21*x31*x41    
c
      return
      end
      function theta_r12p5 (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 54)
      data ncoeff/ 53/
      data avdat/  0.1139717E-02/
      data xmin/
     1 -0.69647E+00,-0.25366E-01,-0.94322E-01,-0.34681E-01,-0.99237E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.61730E+00, 0.38656E-01, 0.94322E-01, 0.34681E-01, 0.99987E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.13018986E-01,-0.17989747E-02,-0.67721188E-01, 0.30841276E-01,
     +  0.90162903E-02,-0.23660799E-04,-0.45266105E-02, 0.13382224E-02,
     + -0.14516278E-02, 0.11692584E-02,-0.46065749E-03, 0.27639263E-02,
     + -0.52839547E-03,-0.20396975E-02,-0.65947877E-03, 0.22226393E-02,
     +  0.25520648E-02,-0.21277766E-02, 0.65209856E-03,-0.13294400E-03,
     + -0.81532943E-03,-0.10638921E-03, 0.52157114E-02, 0.52906678E-03,
     +  0.48528495E-03,-0.47454936E-03,-0.11189295E-02,-0.86779119E-02,
     + -0.10897228E-01,-0.58184280E-02,-0.20042718E-02, 0.59534814E-02,
     +  0.29081251E-01, 0.79838256E-03, 0.12553356E-02,-0.11175508E-02,
     +  0.55435393E-03,-0.11137167E-01,-0.57937768E-02, 0.26353071E-01,
     + -0.14629376E-01, 0.86127373E-03,-0.69670222E-03,-0.16595932E-02,
     + -0.88354503E-03,-0.12106807E-02,-0.66828611E-03, 0.12143422E-01,
     +  0.18520562E-01,-0.10636737E+00, 0.48597801E-01,-0.12820365E+00,
     +  0.66467211E-01,
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
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      theta_r12p5 =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)    *x21            
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x12                
     6  +coeff(  6)            *x42    
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)    *x21        *x51
      theta_r12p5 =theta_r12p5 
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11*x22            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)*x11        *x42    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)*x11*x21        *x51
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)            *x42*x51
      theta_r12p5 =theta_r12p5 
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x22*x31*x41    
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)    *x21    *x42*x51
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)    *x23        *x51
     7  +coeff( 25)    *x21*x31*x41*x51
     8  +coeff( 26)    *x22        *x52
      theta_r12p5 =theta_r12p5 
     9  +coeff( 27)    *x23    *x42    
     1  +coeff( 28)    *x22    *x42*x51
     2  +coeff( 29)    *x23    *x42*x51
     3  +coeff( 30)    *x22            
     4  +coeff( 31)*x13                
     5  +coeff( 32)    *x21    *x42*x52
     6  +coeff( 33)    *x22    *x42*x52
     7  +coeff( 34)        *x32        
     8  +coeff( 35)*x12*x21            
      theta_r12p5 =theta_r12p5 
     9  +coeff( 36)*x11    *x31*x41    
     1  +coeff( 37)        *x32    *x51
     2  +coeff( 38)    *x21    *x42*x53
     3  +coeff( 39)*x13*x22    *x42    
     4  +coeff( 40)    *x23    *x42*x52
     5  +coeff( 41)    *x22    *x42*x53
     6  +coeff( 42)        *x31*x41    
     7  +coeff( 43)*x11    *x32        
     8  +coeff( 44)        *x31*x41*x51
      theta_r12p5 =theta_r12p5 
     9  +coeff( 45)*x11    *x31*x41*x51
     1  +coeff( 46)*x11        *x42*x51
     2  +coeff( 47)            *x42*x53
     3  +coeff( 48)    *x22*x31*x43*x51
     4  +coeff( 49)    *x23*x31*x43*x51
     5  +coeff( 50)    *x22*x31*x43*x52
     6  +coeff( 51)    *x21*x31*x43*x53
     7  +coeff( 52)    *x23*x31*x43*x52
     8  +coeff( 53)    *x22*x31*x43*x53
c
      return
      end
      function phi_r12p5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 33)
      data ncoeff/ 32/
      data avdat/ -0.1441005E-03/
      data xmin/
     1 -0.69647E+00,-0.25366E-01,-0.94322E-01,-0.33747E-01,-0.99237E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.61730E+00, 0.38656E-01, 0.91823E-01, 0.34681E-01, 0.99987E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.54512545E-03, 0.13134193E-01,-0.32253928E-01,-0.15507543E-01,
     + -0.14199584E-02,-0.86875809E-02, 0.11640772E-01, 0.54584993E-02,
     +  0.49686572E-02, 0.67184470E-02, 0.44761444E-02,-0.40542046E-02,
     + -0.12952300E-02, 0.36826448E-02,-0.20279777E-02,-0.17131567E-02,
     + -0.29164685E-01,-0.31024863E-02, 0.14796787E-02,-0.28886616E-02,
     + -0.67139226E-02, 0.31597917E-02,-0.31096593E-02, 0.15055757E-02,
     +  0.22910985E-02,-0.16289565E-02,-0.13834921E-03,-0.68066741E-03,
     + -0.10443472E-02, 0.27242964E-02, 0.59991929E-03,-0.44662245E-02,
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
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      phi_r12p5   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)        *x31        
     4  +coeff(  4)*x11        *x41    
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)*x11    *x31        
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)*x12        *x41    
      phi_r12p5   =phi_r12p5   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)*x12    *x31        
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)        *x32*x41    
     6  +coeff( 15)    *x21    *x43    
     7  +coeff( 16)    *x22    *x43    
     8  +coeff( 17)*x12*x21*x32*x43    
      phi_r12p5   =phi_r12p5   
     9  +coeff( 18)*x11*x21*x31        
     1  +coeff( 19)        *x33        
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)*x11*x22*x31        
     4  +coeff( 22)        *x31*x44    
     5  +coeff( 23)*x13*x21    *x41    
     6  +coeff( 24)*x11    *x31    *x51
     7  +coeff( 25)    *x22    *x41*x51
     8  +coeff( 26)    *x22*x31    *x51
      phi_r12p5   =phi_r12p5   
     9  +coeff( 27)    *x21            
     1  +coeff( 28)*x13        *x41    
     2  +coeff( 29)*x13    *x31        
     3  +coeff( 30)*x12*x22    *x41    
     4  +coeff( 31)*x11        *x41*x51
     5  +coeff( 32)*x11*x23*x31        
c
      return
      end
      function y00_r12p5   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 81)
      data ncoeff/ 80/
      data avdat/ -0.1543536E-02/
      data xmin/
     1 -0.69647E+00,-0.25366E-01,-0.94322E-01,-0.33747E-01,-0.99237E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.61730E+00, 0.38656E-01, 0.91823E-01, 0.34681E-01, 0.99987E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11014852E-02, 0.28395487E-01, 0.49620355E-03,-0.21119084E-04,
     +  0.89608967E-01,-0.42480439E-01,-0.15857782E-01,-0.12434438E-03,
     +  0.58140140E-02, 0.35755064E-01, 0.18268752E-03, 0.60576349E-04,
     + -0.70810947E-02,-0.21063462E-01, 0.51112130E-01,-0.22051633E-03,
     +  0.75370688E-02,-0.80817612E-02, 0.16385472E-02,-0.11838055E-01,
     +  0.29290779E-02, 0.54737912E-02, 0.46160235E-02, 0.24955336E-01,
     + -0.11128540E-01, 0.10230500E-02, 0.77975444E-02,-0.34559972E-02,
     + -0.11674209E-01, 0.66808020E-02,-0.61562792E-02,-0.42555798E-02,
     +  0.12602223E-01,-0.71947174E-02,-0.59552938E-02,-0.13140393E-01,
     + -0.40968922E-02, 0.54361154E-02,-0.29526188E-02, 0.47562667E-02,
     +  0.15584632E-02,-0.25094401E-01, 0.55150143E-02,-0.20912884E-01,
     + -0.89343387E-03, 0.33375998E-02, 0.28162424E-02,-0.95712729E-02,
     +  0.49631624E-02, 0.48779980E-02,-0.93566012E-02, 0.18782029E-02,
     +  0.16377709E-02, 0.48007304E-03, 0.13612343E-03, 0.26129487E-02,
     +  0.11158281E-01,-0.15605728E-01, 0.34412337E-02,-0.56789936E-02,
     +  0.29238684E-02, 0.21202663E-01,-0.15846625E-01, 0.13622635E-02,
     +  0.19074230E-02, 0.26321040E-01, 0.24950968E-02, 0.38518905E-03,
     +  0.14737046E-03, 0.17412951E-02, 0.38872431E-02,-0.70210737E-02,
     +  0.14970715E-02,-0.14919582E-02,-0.10284946E-02,-0.40725188E-03,
     +  0.18648041E-02,-0.21296130E-02,-0.84341420E-02,-0.39403578E-02,
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
      y00_r12p5   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x41    
     6  +coeff(  6)*x11    *x31        
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x31*x41    
      y00_r12p5   =y00_r12p5   
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)*x12    *x31        
     5  +coeff( 14)*x11*x21*x31        
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)        *x31*x42    
      y00_r12p5   =y00_r12p5   
     9  +coeff( 18)        *x31    *x51
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)            *x41*x51
     5  +coeff( 23)    *x21*x32*x41    
     6  +coeff( 24)    *x23*x31        
     7  +coeff( 25)    *x21*x31    *x51
     8  +coeff( 26)    *x24            
      y00_r12p5   =y00_r12p5   
     9  +coeff( 27)*x11*x22    *x41    
     1  +coeff( 28)*x13*x21    *x41    
     2  +coeff( 29)*x12        *x41    
     3  +coeff( 30)*x13    *x31        
     4  +coeff( 31)*x11    *x31    *x51
     5  +coeff( 32)*x11        *x41*x51
     6  +coeff( 33)    *x21    *x41*x51
     7  +coeff( 34)*x11*x21    *x41*x51
     8  +coeff( 35)            *x41*x52
      y00_r12p5   =y00_r12p5   
     9  +coeff( 36)*x11*x21    *x41    
     1  +coeff( 37)            *x43    
     2  +coeff( 38)*x13        *x41    
     3  +coeff( 39)    *x22    *x42    
     4  +coeff( 40)*x11        *x43    
     5  +coeff( 41)    *x24*x31        
     6  +coeff( 42)    *x23    *x42    
     7  +coeff( 43)    *x22    *x41*x51
     8  +coeff( 44)    *x24    *x42    
      y00_r12p5   =y00_r12p5   
     9  +coeff( 45)        *x33    *x52
     1  +coeff( 46)        *x33        
     2  +coeff( 47)*x11    *x33        
     3  +coeff( 48)*x12*x21*x31        
     4  +coeff( 49)*x11    *x31*x42    
     5  +coeff( 50)*x12*x21    *x41    
     6  +coeff( 51)*x12*x22*x31        
     7  +coeff( 52)        *x31    *x52
     8  +coeff( 53)*x12    *x33    *x51
      y00_r12p5   =y00_r12p5   
     9  +coeff( 54)*x11*x21            
     1  +coeff( 55)                *x51
     2  +coeff( 56)    *x21*x33        
     3  +coeff( 57)*x11*x22*x31        
     4  +coeff( 58)*x11*x21*x32*x41    
     5  +coeff( 59)*x14    *x31        
     6  +coeff( 60)*x11*x21*x31    *x51
     7  +coeff( 61)*x14        *x41    
     8  +coeff( 62)    *x21*x33*x42    
      y00_r12p5   =y00_r12p5   
     9  +coeff( 63)*x11*x22*x32*x41    
     1  +coeff( 64)*x11    *x31    *x52
     2  +coeff( 65)*x11        *x41*x52
     3  +coeff( 66)    *x22*x33*x42    
     4  +coeff( 67)    *x22    *x41*x52
     5  +coeff( 68)        *x32*x41    
     6  +coeff( 69)    *x21        *x51
     7  +coeff( 70)*x12    *x33        
     8  +coeff( 71)*x12    *x32*x41    
      y00_r12p5   =y00_r12p5   
     9  +coeff( 72)*x13*x21*x31        
     1  +coeff( 73)*x12    *x31    *x51
     2  +coeff( 74)    *x24    *x41    
     3  +coeff( 75)*x12        *x41*x51
     4  +coeff( 76)        *x34    *x51
     5  +coeff( 77)    *x21*x33    *x51
     6  +coeff( 78)*x13*x21*x32        
     7  +coeff( 79)*x13*x22*x31        
     8  +coeff( 80)*x11*x24*x31        
c
      return
      end
