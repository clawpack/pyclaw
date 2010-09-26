c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  ilr,asdq,bmasdq,bpasdq)
c     =====================================================
      use geoclaw_module

      implicit none
c
c     # Riemann solver in the transverse direction using an einfeldt
c     Jacobian.

c-----------------------last modified 1/10/05----------------------

      integer ixy,maxm,meqn,mwaves,mbc,mx,ilr

      double precision  ql(1-mbc:maxm+mbc, meqn)
      double precision  qr(1-mbc:maxm+mbc, meqn)
      double precision  asdq(1-mbc:maxm+mbc, meqn)
      double precision  bmasdq(1-mbc:maxm+mbc, meqn)
      double precision  bpasdq(1-mbc:maxm+mbc, meqn)
      double precision  aux1(1-mbc:maxm+mbc,*)
      double precision  aux2(1-mbc:maxm+mbc,*)
      double precision  aux3(1-mbc:maxm+mbc,*)

      double precision  s(3)
      double precision  r(3,3)
      double precision  beta(3)
      double precision  g,tol,abs_tol
      double precision  hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br
      double precision  uhat,vhat,hhat,roe1,roe3,s1,s2,s3,s1l,s3r
      double precision  delf1,delf2,delf3,dxdcd,dxdcu
      double precision  dxdcm,dxdcp,topo1,topo3,eta

      integer i,m,mw,mu,mv

      g=grav
      tol=drytolerance
      abs_tol=drytolerance

      if (ixy.eq.1) then
	     mu = 2
	     mv = 3
      else
	     mu = 3
	     mv = 2
      endif


      do i=2-mbc,mx+mbc

         hl=qr(i-1,1)
         hr=ql(i,1)
         hul=qr(i-1,mu)
         hur=ql(i,mu)
         hvl=qr(i-1,mv)
         hvr=ql(i,mv)

c===========determine velocity from momentum===========================
       if (hl.lt.abs_tol) then
          hl=0.d0
          ul=0.d0
          vl=0.d0
       else
          ul=hul/hl
          vl=hvl/hl
       endif

       if (hr.lt.abs_tol) then
          hr=0.d0
          ur=0.d0
          vr=0.d0
       else
          ur=hur/hr
          vr=hvr/hr
       endif

       do mw=1,mwaves
          s(mw)=0.d0
          beta(mw)=0.d0
          do m=1,meqn
             r(m,mw)=0.d0
          enddo
       enddo
      dxdcp = 1.d0
      dxdcm = 1.d0

       if (hl.le.drytolerance.and.hr.le.drytolerance) go to 90

*      !check and see if cell that transverse waves are going in is high and dry
       if (ilr.eq.1) then
            eta = qr(i-1,1) + aux2(i-1,1)
            topo1 = aux1(i-1,1)
            topo3 = aux3(i-1,1)
c            s1 = vl-sqrt(g*hl)
c            s3 = vl+sqrt(g*hl)
c            s2 = 0.5d0*(s1+s3)
       else
            eta = ql(i,1) + aux2(i,1)
            topo1 = aux1(i,1)
            topo3 = aux3(i,1)
c            s1 = vr-sqrt(g*hr)
c            s3 = vr+sqrt(g*hr)
c            s2 = 0.5d0*(s1+s3)
       endif
       if (eta.lt.max(topo1,topo3)) go to 90

      if (icoordsys.eq.2) then
         if (ixy.eq.2) then
            dxdcp=(Rearth*pi/180.d0)
            dxdcm = dxdcp
         else
            if (ilr.eq.1) then
               dxdcp = Rearth*pi*cos(aux3(i-1,3))/180.d0
               dxdcm = Rearth*pi*cos(aux1(i-1,3))/180.d0
            else
               dxdcp = Rearth*pi*cos(aux3(i,3))/180.d0
               dxdcm = Rearth*pi*cos(aux1(i,3))/180.d0
            endif
         endif
      endif

c=====Determine some speeds necessary for the Jacobian=================
            vhat=(vr*dsqrt(hr))/(dsqrt(hr)+dsqrt(hl)) +
     &        (vl*dsqrt(hl))/(dsqrt(hr)+dsqrt(hl))

            uhat=(ur*dsqrt(hr))/(dsqrt(hr)+dsqrt(hl)) +
     &        (ul*dsqrt(hl))/(dsqrt(hr)+dsqrt(hl))
            hhat=(hr+hl)/2.d0

            roe1=vhat-dsqrt(g*hhat)
            roe3=vhat+dsqrt(g*hhat)

            s1l=vl-dsqrt(g*hl)
            s3r=vr+dsqrt(g*hr)

            s1=dmin1(roe1,s1l)
            s3=dmax1(roe3,s3r)

            s2=0.5d0*(s1+s3)

           s(1)=s1
           s(2)=s2
           s(3)=s3
c=======================Determine asdq decomposition (beta)============
         delf1=asdq(i,1)
         delf2=asdq(i,mu)
         delf3=asdq(i,mv)

         beta(1) = (s3*delf1/(s3-s1))-(delf3/(s3-s1))
         beta(2) = -s2*delf1 + delf2
         beta(3) = (delf3/(s3-s1))-(s1*delf1/(s3-s1))
c======================End =================================================

c=====================Set-up eigenvectors===================================
         r(1,1) = 1.d0
         r(2,1) = s2
         r(3,1) = s1

         r(1,2) = 0.d0
         r(2,2) = 1.d0
         r(3,2) = 0.d0

         r(1,3) = 1.d0
         r(2,3) = s2
         r(3,3) = s3
c============================================================================
90      continue
c============= compute fluctuations==========================================

            do  m=1,meqn
               bmasdq(i,m)=0.0d0
               bpasdq(i,m)=0.0d0
            enddo
            do  mw=1,3
               if (s(mw).lt.0.d0) then
                 bmasdq(i,1) =bmasdq(i,1) + dxdcm*s(mw)*beta(mw)*r(1,mw)
                 bmasdq(i,mu)=bmasdq(i,mu)+ dxdcm*s(mw)*beta(mw)*r(2,mw)
                 bmasdq(i,mv)=bmasdq(i,mv)+ dxdcm*s(mw)*beta(mw)*r(3,mw)
               elseif (s(mw).gt.0.d0) then
                 bpasdq(i,1) =bpasdq(i,1) + dxdcp*s(mw)*beta(mw)*r(1,mw)
                 bpasdq(i,mu)=bpasdq(i,mu)+ dxdcp*s(mw)*beta(mw)*r(2,mw)
                 bpasdq(i,mv)=bpasdq(i,mv)+ dxdcp*s(mw)*beta(mw)*r(3,mw)
               endif
            enddo
c========================================================================
         enddo
c

c
      return
      end
