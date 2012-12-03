      subroutine hlmhltz(x,y,z,r,amp,bx,by,bz)
c************************************************************************
c hlmhltz - calculates field for a pair of helmholtz coils
c           x,y,z = field point coords, origin exactly between coils (meters)
c               r = coil radius and separation
c             amp = current in coils
c              bi = x,y,z components of field (tesla)
c 7/8/87
c*************************************************************************
      real x,y,z,r,amp,bx,by,bz,rho,br1,br2,by1,by2,theta,y1,y2
c      write(6,*)"in hlmhltz",r, amp
      y1=(r/2)+y
      y2=y-(r/2)
      rho=sqrt(x**2+z**2)
      if (rho.eq.0.)then
      theta=0.
      go to 1
      endif
      theta=atan2(z,x)
 1    call loop(r,amp,rho,y1,br1,by1)
      call loop(r,amp,rho,y2,br2,by2)
      by=by1+by2
      br=br1+br2
      bx=br*cos(theta)
      bz=br*sin(theta)
      return
      end

      subroutine loop(a,amp,rho,y,br,by)
c************************************************************************
c  loop - calculates the magnetic field for a current loop
c         a = radius of loop
c       amp = current in loop
c       rho = distance from loop axis (meters)
c         y = distance from loop along axis (meters)
c        br = radial component of field (tesla)
c        by = axial component of field
c
c  reference: garrett, journal of applied physics, 34 sept 1963, p. 2567
c  7/8/87
c******************************************************************
      real a,amp,rho,y,br,by,r1sq,r2sq,ksq,kpsq
c	real vk,ve
      real*8 dksq,dvk,dve
c evaluate parameters
      r1sq=(a+rho)**2+y**2
      r2sq=(a-rho)**2+y**2
      ksq=4.*a*rho/r1sq
      kpsq=1-ksq
c evaluate elliptical integrals
      dksq=ksq
      call fb01ad(dksq,dvk,dve)
c      ve=dve
c      vk=dvk
c      write(6,*)ksq,vk,ve
c calculate fields
      if(r2sq.eq.0.)then
c      write(6,*)' bull''s-eye!!! you just hit a coil.'
      by=1.69e+38
      go to 15
      endif
      by=(amp*2.0e-07/sqrt(r1sq))*((dvk-dve)+(a-rho)*2*a*dve/r2sq)
 15   if(rho.eq.0.)then
      br=0.
      return
      endif
      if(kpsq.eq.0.)kpsq=0.3e-36
      br=-(amp*1.e-07*y)/(rho*sqrt(r1sq)*kpsq)*((2-ksq)*(dvk-dve)
     &     -ksq*dvk)
      return
      end      

      subroutine fb01ad(c,  vk,ve)
      implicit real*8(a-h,o-z)
c*ibm real*8 xlg/  z7fffffffffffffff /
      real * 8 xlg/'7fffffffffffffff'x/
      d=1d0-c
      if(d .gt. 0d0) e=-log(d)
c**** harwell version of fb01ad
      if(c .ge. 1d0)go to 2
           ve=e*((((((((((
     a     3.18591956555015718d-5*d  +.989833284622538479d-3)*d
     b    +.643214658643830177d-2)*d +.16804023346363385d-1)*d
     c    +.261450147003138789d-1)*d +.334789436657616262d-1)*d
     d    +.427178905473830956d-1)*d +.585936612555314917d-1)*d
     e    +.937499997212031407d-1)*d +.249999999999901772d0)*d)
     f    +(((((((((
     g     .149466217571813268d-3*d  +.246850333046072273d-2)*d
     h    +.863844217360407443d-2)*d+.107706350398664555d-1)*d
     i    +.782040406095955417d-2)*d +.759509342255943228d-2)*d
     j    +.115695957452954022d-1)*d +.218318116761304816d-1)*d
     k    +.568051945675591566d-1)*d +.443147180560889526d0)*d
     l    +1d0
c****
c**** routine modified to calculate vk and ve always
c****
c****
           vk=e*((((((((((
     a     .297002809665556121d-4*d   +.921554634963249846d-3)*d
     b    +.597390429915542916d-2)*d  +.155309416319772039d-1)*d
     c    +.239319133231107901d-1)*d  +.301248490128989303d-1)*d
     d    +.373777397586236041d-1)*d  +.48828041906862398d-1)*d
     e    +.703124997390383521d-1)*d  +.124999999999908081d0)*d
     f    +.5d0)+(((((((((
     g     .139308785700664673d-3*d   +.229663489839695869d-2)*d
     h    +.800300398064998537d-2)*d  +.984892932217689377d-2)*d
     i    +.684790928262450512d-2)*d  +.617962744605331761d-2)*d
     j    +.878980187455506468d-2)*d  +.149380135326871652d-1)*d
     k    +.308851462713051899d-1)*d  +.965735902808562554d-1)*d
     l    +1.38629436111989062d0
      return
    2 ve=1d0
      vk=xlg
      return
      end
