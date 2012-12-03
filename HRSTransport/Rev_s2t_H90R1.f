c  Transport function fitted on trajectories from SNAKE. Here the trajectories combine from 1, 2, 4GeV, with 5T g2p target field. 
c  (x,m) is xf, theatf, yf, phif, p, note this is different from normal case. p is the momentum of the particle.  

      function s2t_sr5p69_2gev_txfit                        (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  4)
      data ncoeff/  3/
      data avdat/ -0.2055424E-02/
      data xmin/
     1 -0.60683E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.57697E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29396296E-02, 0.47811743E-01, 0.16085932E-03,
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
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
c
c                  function
c
      s2t_sr5p69_2gev_txfit                        =avdat
     1  +coeff(  1)    
     2  +coeff(  2)*x11
     3  +coeff(  3)*x12
c
      return
      end
      function s2t_sr5p69_2gev_pyfit                        (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  3)
      data ncoeff/  2/
      data avdat/  0.1112069E+00/
      data xmin/
     1  0.94655E-01, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15321E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.31114595E-02, 0.24103777E-01,
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
c          set up monomials   functions
      x11 = x1
c
c                  function
c
      s2t_sr5p69_2gev_pyfit                        =avdat
     1  +coeff(  1)    
     2  +coeff(  2)*x11
c
      return
      end
   

      function s2t_sr5p69_theta                     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1266531E+00/
      data xmin/
     1 -0.62411E-01,-0.18267E-01,-0.59419E-01,-0.11117E+00, 0.95543E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15554E+00, 0.39464E-01, 0.15700E+00, 0.92662E-02, 0.41993E+01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23894185E-01,-0.37870798E+01,-0.46646260E-01,-0.40849119E+00,
     +  0.40399462E+00, 0.41269126E+01, 0.27491945E+00,-0.22369785E+01,
     +  0.66102070E+00, 0.28629756E+01, 0.33461137E-02, 0.61817592E+00,
     + -0.17541131E+00, 0.13850525E+00, 0.75335629E-01,-0.21910113E+00,
     + -0.18929407E+00, 0.40399574E-01,-0.70259482E-01, 0.10015142E+00,
     +  0.19589485E+00, 0.11737118E+01,-0.18363287E+00,-0.12538783E+02,
     +  0.65080633E+01,-0.12424981E+02, 0.14241733E+02,-0.12961671E+01,
     + -0.54516110E+01,-0.93346238E+00, 0.15955850E+01, 0.72892475E+01,
     +  0.72293550E-01,-0.20286281E-01,-0.81209309E-01, 0.12778061E+02,
     +  0.96893275E+00,-0.93060148E+00,-0.62582213E-02, 0.57337868E+00,
     +  0.22097619E+00,-0.80936015E+00,-0.27538414E+01, 0.12162517E+01,
     +  0.29791622E+01, 0.10291787E+02, 0.13519272E+01,-0.33758788E+01,
     +  0.37683630E+01,-0.72251196E+01,-0.19234230E+00, 0.20134732E+00,
     +  0.12415381E+00,-0.35109486E-01,-0.74394554E-01, 0.15084922E+00,
     + -0.12243018E+02, 0.65039791E-01, 0.45359921E+00,-0.14984195E+01,
     + -0.42374745E+00,-0.54812393E+01,-0.41083574E+01,-0.14658412E+00,
     + -0.68840594E+01,-0.55978066E+00, 0.21097827E+01,-0.58427796E+01,
     + -0.39053507E-01, 0.16226881E+01, 0.12902580E+01,-0.24011185E+01,
     +  0.15092519E+00,-0.19721635E+00,-0.63920755E+01, 0.52710290E+01,
     +  0.65528497E-01, 0.55205064E+01, 0.98444023E+01,-0.13241970E-01,
     +  0.44322562E+00, 0.24767270E+01, 0.33296292E+01,-0.17606243E+01,
     + -0.26935276E-01, 0.19425501E+00, 0.43221036E+00, 0.18972932E+00,
     +  0.87853782E-02,-0.72694910E-02,
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
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      s2t_sr5p69_theta                     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)*x11                
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x12                
     5  +coeff(  5)        *x31        
     6  +coeff(  6)                *x51
     7  +coeff(  7)*x13                
     8  +coeff(  8)    *x23            
      s2t_sr5p69_theta                     =s2t_sr5p69_theta                     
     9  +coeff(  9)*x11    *x31        
     1  +coeff( 10)*x11            *x51
     2  +coeff( 11)*x14                
     3  +coeff( 12)*x12    *x31        
     4  +coeff( 13)        *x32        
     5  +coeff( 14)*x12            *x51
     6  +coeff( 15)        *x31    *x51
     7  +coeff( 16)                *x52
     8  +coeff( 17)*x13    *x31        
      s2t_sr5p69_theta                     =s2t_sr5p69_theta                     
     9  +coeff( 18)*x11            *x52
     1  +coeff( 19)*x13            *x51
     2  +coeff( 20)*x14    *x31        
     3  +coeff( 21)*x12    *x32        
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)*x14    *x32        
     6  +coeff( 24)*x11*x22            
     7  +coeff( 25)*x11*x21        *x51
     8  +coeff( 26)    *x21    *x41*x51
      s2t_sr5p69_theta                     =s2t_sr5p69_theta                     
     9  +coeff( 27)*x11*x21    *x42    
     1  +coeff( 28)    *x21        *x52
     2  +coeff( 29)*x12*x22    *x41    
     3  +coeff( 30)*x12    *x31*x41    
     4  +coeff( 31)*x11*x23        *x51
     5  +coeff( 32)    *x24        *x51
     6  +coeff( 33)*x12    *x31    *x51
     7  +coeff( 34)        *x32    *x51
     8  +coeff( 35)        *x31    *x52
      s2t_sr5p69_theta                     =s2t_sr5p69_theta                     
     9  +coeff( 36)    *x23    *x41*x51
     1  +coeff( 37)    *x21    *x41*x52
     2  +coeff( 38)*x13*x22        *x51
     3  +coeff( 39)*x13            *x52
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)*x14*x22        *x51
     6  +coeff( 42)*x13*x23        *x51
     7  +coeff( 43)*x12*x24        *x51
     8  +coeff( 44)*x12*x22        *x52
      s2t_sr5p69_theta                     =s2t_sr5p69_theta                     
     9  +coeff( 45)*x11        *x41    
     1  +coeff( 46)    *x21        *x51
     2  +coeff( 47)*x13*x21            
     3  +coeff( 48)*x11        *x41*x51
     4  +coeff( 49)*x13*x22            
     5  +coeff( 50)*x11*x23    *x41    
     6  +coeff( 51)        *x31*x41*x51
     7  +coeff( 52)*x11    *x31*x41*x51
     8  +coeff( 53)*x11    *x31    *x52
      s2t_sr5p69_theta                     =s2t_sr5p69_theta                     
     9  +coeff( 54)*x12        *x41*x52
     1  +coeff( 55)*x14    *x31    *x51
     2  +coeff( 56)*x12    *x32    *x51
     3  +coeff( 57)*x11*x21            
     4  +coeff( 58)            *x42    
     5  +coeff( 59)*x12*x21            
     6  +coeff( 60)*x12        *x41    
     7  +coeff( 61)        *x31*x41    
     8  +coeff( 62)            *x41*x51
      s2t_sr5p69_theta                     =s2t_sr5p69_theta                     
     9  +coeff( 63)    *x24            
     1  +coeff( 64)*x13        *x41    
     2  +coeff( 65)*x12*x21    *x41    
     3  +coeff( 66)*x11    *x31*x41    
     4  +coeff( 67)*x12        *x42    
     5  +coeff( 68)*x12*x23            
     6  +coeff( 69)*x11    *x32        
     7  +coeff( 70)*x11        *x43    
     8  +coeff( 71)*x11*x22        *x51
      s2t_sr5p69_theta                     =s2t_sr5p69_theta                     
     9  +coeff( 72)    *x23        *x51
     1  +coeff( 73)        *x32*x41    
     2  +coeff( 74)*x12        *x41*x51
     3  +coeff( 75)*x11*x21    *x41*x51
     4  +coeff( 76)    *x22    *x41*x51
     5  +coeff( 77)*x13        *x42    
     6  +coeff( 78)*x12*x21    *x42    
     7  +coeff( 79)*x11*x22    *x42    
     8  +coeff( 80)    *x21*x31*x42    
      s2t_sr5p69_theta                     =s2t_sr5p69_theta                     
     9  +coeff( 81)*x11        *x42*x51
     1  +coeff( 82)    *x21    *x42*x51
     2  +coeff( 83)*x13*x23            
     3  +coeff( 84)*x12*x24            
     4  +coeff( 85)*x13*x21*x31        
     5  +coeff( 86)*x12*x22*x31        
     6  +coeff( 87)*x11*x23*x31        
     7  +coeff( 88)    *x24*x31        
     8  +coeff( 89)*x11*x21*x32        
      s2t_sr5p69_theta                     =s2t_sr5p69_theta                     
     9  +coeff( 90)    *x22*x32        
c
      return
      end
      function s2t_sr5p69_phi                       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 76)
      data ncoeff/ 75/
      data avdat/ -0.9588001E-03/
      data xmin/
     1 -0.62411E-01,-0.18267E-01,-0.59419E-01,-0.11117E+00, 0.95543E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15554E+00, 0.39464E-01, 0.15700E+00, 0.92662E-02, 0.41993E+01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21962149E+01,-0.38783202E+01,-0.48051329E+01, 0.31108800E+00,
     +  0.25762427E+01, 0.11257217E+02, 0.54453883E+01, 0.24435802E+01,
     + -0.18620904E+00,-0.76644903E+00,-0.55179262E+00,-0.45624261E+01,
     + -0.20621407E+02,-0.42789366E-01,-0.48589492E+00,-0.18619319E+01,
     + -0.12432082E+02, 0.45184299E+00,-0.60865195E-02, 0.19494985E+00,
     +  0.17531042E+01, 0.86988807E+00, 0.20006139E+01, 0.22968011E+00,
     + -0.19177313E+01,-0.29976417E-02,-0.15591797E+00,-0.94863717E-02,
     + -0.49013373E-01, 0.48813593E-01, 0.17802741E+00, 0.26373247E-01,
     + -0.29682571E-02, 0.22931237E+01, 0.35125401E+02, 0.90475864E+01,
     + -0.45222473E+01,-0.15414256E+00,-0.26051058E+02, 0.22800982E+01,
     +  0.12079324E+01, 0.17930396E+00,-0.16525537E-01, 0.15274862E+01,
     + -0.76223058E+00, 0.21809661E+00, 0.91049933E+00, 0.77390009E+00,
     +  0.31818995E+00, 0.42038538E-01, 0.52685934E+00, 0.12156507E+00,
     + -0.32062411E+00,-0.15039532E+00,-0.10449097E+00, 0.97089045E-01,
     +  0.28345326E-01,-0.17734368E+00, 0.24734035E-01,-0.25522583E-02,
     +  0.38652144E-01,-0.40646534E-01,-0.24220709E-01,-0.15945429E+00,
     +  0.45034323E-01, 0.21324148E-01,-0.60722616E-01, 0.46743415E-02,
     + -0.37379045E-01,-0.67474484E-01,-0.37250713E-01,-0.11170021E-01,
     +  0.20563926E-02, 0.56950646E-02,-0.88832555E-02,
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
      s2t_sr5p69_phi                       =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      s2t_sr5p69_phi                       =s2t_sr5p69_phi                       
     9  +coeff(  9)*x11    *x31        
     1  +coeff( 10)        *x33        
     2  +coeff( 11)*x11        *x41    
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)        *x32*x41    
     5  +coeff( 14)*x12                
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11    *x32        
     8  +coeff( 17)        *x34        
      s2t_sr5p69_phi                       =s2t_sr5p69_phi                       
     9  +coeff( 18)*x11    *x31*x41    
     1  +coeff( 19)                *x51
     2  +coeff( 20)*x12    *x31        
     3  +coeff( 21)*x11    *x33        
     4  +coeff( 22)        *x33*x42    
     5  +coeff( 23)    *x21    *x43    
     6  +coeff( 24)*x12    *x32        
     7  +coeff( 25)*x11    *x33*x41    
     8  +coeff( 26)*x11            *x51
      s2t_sr5p69_phi                       =s2t_sr5p69_phi                       
     9  +coeff( 27)*x12    *x34        
     1  +coeff( 28)*x12            *x51
     2  +coeff( 29)        *x31    *x51
     3  +coeff( 30)            *x41*x51
     4  +coeff( 31)*x12*x21            
     5  +coeff( 32)        *x33    *x51
     6  +coeff( 33)*x11    *x32    *x51
     7  +coeff( 34)*x11    *x32*x41    
     8  +coeff( 35)        *x34*x41    
      s2t_sr5p69_phi                       =s2t_sr5p69_phi                       
     9  +coeff( 36)        *x32*x43    
     1  +coeff( 37)        *x31*x44    
     2  +coeff( 38)*x11    *x34        
     3  +coeff( 39)        *x34*x42    
     4  +coeff( 40)        *x34*x43    
     5  +coeff( 41)        *x33*x44    
     6  +coeff( 42)        *x34*x44    
     7  +coeff( 43)        *x34    *x51
     8  +coeff( 44)        *x31*x42    
      s2t_sr5p69_phi                       =s2t_sr5p69_phi                       
     9  +coeff( 45)        *x33*x41    
     1  +coeff( 46)    *x21    *x42    
     2  +coeff( 47)        *x32*x42    
     3  +coeff( 48)        *x31*x43    
     4  +coeff( 49)            *x44    
     5  +coeff( 50)*x12        *x41    
     6  +coeff( 51)*x11*x21    *x41    
     7  +coeff( 52)*x11        *x43    
     8  +coeff( 53)*x12    *x31*x41    
      s2t_sr5p69_phi                       =s2t_sr5p69_phi                       
     9  +coeff( 54)*x11    *x32*x42    
     1  +coeff( 55)*x11    *x31*x43    
     2  +coeff( 56)        *x33*x43    
     3  +coeff( 57)*x11        *x44    
     4  +coeff( 58)    *x21    *x44    
     5  +coeff( 59)        *x32*x44    
     6  +coeff( 60)        *x32    *x51
     7  +coeff( 61)        *x31*x41*x51
     8  +coeff( 62)            *x42*x51
      s2t_sr5p69_phi                       =s2t_sr5p69_phi                       
     9  +coeff( 63)*x12    *x33        
     1  +coeff( 64)*x12*x21    *x41    
     2  +coeff( 65)*x11    *x34*x41    
     3  +coeff( 66)*x12    *x31*x42    
     4  +coeff( 67)*x11    *x33*x42    
     5  +coeff( 68)*x12        *x43    
     6  +coeff( 69)*x11*x21    *x43    
     7  +coeff( 70)*x11    *x32*x43    
     8  +coeff( 71)*x11    *x31*x44    
      s2t_sr5p69_phi                       =s2t_sr5p69_phi                       
     9  +coeff( 72)    *x21*x31*x44    
     1  +coeff( 73)*x11    *x31    *x51
     2  +coeff( 74)    *x21*x31    *x51
     3  +coeff( 75)*x11        *x41*x51
c
      return
      end
      function s2t_sr5p69_y00                       (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 37)
      data ncoeff/ 36/
      data avdat/ -0.1458034E-03/
      data xmin/
     1 -0.62411E-01,-0.18267E-01,-0.59419E-01,-0.11117E+00, 0.95543E+00,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.15554E+00, 0.39464E-01, 0.15700E+00, 0.92662E-02, 0.41993E+01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.34981528E-02, 0.15370658E+00,-0.93168035E-01, 0.12996192E-02,
     + -0.44973017E-03, 0.16209338E-01, 0.81733987E-01, 0.17380847E-01,
     + -0.18419795E+00, 0.11451633E-02,-0.14943263E+00, 0.89332864E-01,
     +  0.39067608E-02, 0.71394013E-03, 0.16621275E-01, 0.16465823E-02,
     +  0.46573453E-01, 0.10169813E-01, 0.80173546E-02,-0.96865378E-01,
     +  0.17356069E+00,-0.33322803E-02,-0.14231360E-01, 0.34248035E-02,
     +  0.49868613E-02, 0.60624443E-02,-0.21534087E-01,-0.36354342E-02,
     +  0.39165621E-03, 0.16561626E-01,-0.86316653E-02, 0.74928678E-02,
     +  0.32582316E-02, 0.42175623E-02,-0.81290705E-02, 0.60676676E-02,
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
c
c                  function
c
      s2t_sr5p69_y00                       =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x12                
     6  +coeff(  6)*x11*x21            
     7  +coeff(  7)    *x22            
     8  +coeff(  8)*x11    *x31        
      s2t_sr5p69_y00                       =s2t_sr5p69_y00                       
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)        *x32        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)*x11*x21*x31        
      s2t_sr5p69_y00                       =s2t_sr5p69_y00                       
     9  +coeff( 18)*x11    *x32        
     1  +coeff( 19)        *x33        
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)*x11    *x31    *x51
     5  +coeff( 23)*x12*x21*x31        
     6  +coeff( 24)*x13    *x34        
     7  +coeff( 25)                *x51
     8  +coeff( 26)        *x31    *x51
      s2t_sr5p69_y00                       =s2t_sr5p69_y00                       
     9  +coeff( 27)    *x21    *x42    
     1  +coeff( 28)*x14                
     2  +coeff( 29)*x11    *x33        
     3  +coeff( 30)    *x21*x31*x41*x51
     4  +coeff( 31)        *x32*x43    
     5  +coeff( 32)        *x34    *x51
     6  +coeff( 33)        *x33    *x52
     7  +coeff( 34)*x13*x21*x31    *x51
     8  +coeff( 35)*x11*x21    *x41*x53
      s2t_sr5p69_y00                       =s2t_sr5p69_y00                       
     9  +coeff( 36)*x12    *x33*x41*x51
c
      return
      end
