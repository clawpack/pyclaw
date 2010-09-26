

c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c      # Set initial sea level flat unless iqinit = 1, in which case
c      # an initial perturbation of the q(i,j,1) is specified and has
c      # been strored in qinitwork.


       use geoclaw_module

       implicit double precision (a-h,o-z)
       dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
       dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

       include 'qinit.i'


       do i=1-mbc,mx+mbc
          x = xlower + (i-0.5d0)*dx
          do j=1-mbc,my+mbc
             y = ylower + (j-0.5d0)*dy
             q(i,j,1)=dmax1(0.d0,sealevel-aux(i,j,1))
             q(i,j,2)=0.d0
             q(i,j,3)=0.d0
             enddo
          enddo

       if (iqinit .gt. 0) then

         do i=1-mbc,mx+mbc
           x = xlower + (i-0.5d0)*dx
           xim = x - 0.5d0*dx
           xip = x + 0.5d0*dx
           do j=1-mbc,my+mbc
             y = ylower + (j-0.5d0)*dy
             yjm = y - 0.5d0*dy
             yjp = y + 0.5d0*dy


             if (xip.gt.xlowqinit.and.xim.lt.xhiqinit
     &          .and.yjp.gt.ylowqinit.and.yjm.lt.yhiqinit) then

                   xipc=min(xip,xhiqinit)
                   ximc=max(xim,xlowqinit)
                   xc=0.5d0*(xipc+ximc)

                   yjpc=min(yjp,yhiqinit)
                   yjmc=max(yjm,ylowqinit)
                   yc=0.5d0*(yjmc+yjpc)

                   dq = topointegral(ximc,xc,xipc,yjmc,yc,yjpc,
     &                    xlowqinit,ylowqinit,dxqinit,dyqinit,
     &                    mxqinit,myqinit,qinitwork,1)
                   dq=dq/((xipc-ximc)*(yjpc-yjmc)*aux(i,j,2))

                   if (iqinit.lt.4) then 
                      if (aux(i,j,1).le.0.d0) then
                          q(i,j,iqinit) = q(i,j,iqinit) + dq
                      endif
                   elseif (iqinit.eq.4) then
                      q(i,j,1) = max(dq-aux(i,j,1),0.d0)
                   endif
             endif
            enddo
           enddo
         endif

       return
       end
