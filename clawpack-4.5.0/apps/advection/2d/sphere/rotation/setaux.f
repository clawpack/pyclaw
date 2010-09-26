c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dxc,dyc,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays for advection on a curvilinear grid
c
c     # on input, (xc(i),yc(j)) gives uniformly spaced computational grid.
c     # on output, 
c     #   aux(i,j,1) is edge velocity at "left" boundary of grid point (i,j)
c     #   aux(i,j,2) is edge velocity at "bottom" boundary of grid point (i,j)
c     #   aux(i,j,3) = kappa  is ratio of cell area to dxc*dyc
c     
      implicit double precision (a-h,o-z)
      parameter (maxm3 = 1005)
      dimension xc(-3:maxm3), yc(-3:maxm3)
      dimension xp(-3:maxm3,-3:maxm3), yp(-3:maxm3,-3:maxm3)
      dimension zp(-3:maxm3,-3:maxm3)
      dimension theta(-3:maxm3,-3:maxm3), phi(-3:maxm3,-3:maxm3)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, 3)
      common /comsphere/ Rsphere

      Rsphere = 1.d0  ! try this
 
      pi = 4.d0*datan(1.d0)
c     write(6,*) 'in setaux',mx,my,xlower,ylower,dxc,dyc
c
      if (mbc .gt. 4) then
	 write(6,*)'***  increase size of local arrays in setaux ***'
	 stop
	 endif
      if (mx+mbc+1.gt.maxm3 .or. my+mbc+1.gt.maxm3) then
	 write(6,*)'***  increase size of maxm3 in setaux ***'
	 stop
	 endif
c
c     # Set xc and yc so that (xc(i),yc(j)) is the 
c     # lower left corner of (i,j) cell in computational space:
c
      do 10 i=1-mbc,mx+mbc+1
         xc(i) = xlower + (i-1.d0) * dxc
   10    continue
c
      do 12 j=1-mbc,my+mbc+1
         yc(j) = ylower + (j-1.d0) * dyc
   12    continue

c     # compute cell corners on sphere and angles phi, theta
c     # related to latitude and longitude
c
      do 15 j=1-mbc,my+mbc+1
         do 15 i=1-mbc,mx+mbc+1

c           # map computational point to (xp,yp,zp) on sphere:
            call mapc2m(xc(i),yc(j),xp(i,j),yp(i,j),zp(i,j))

c           # compute longitude theta from positive x axis:
            r = dsqrt(xp(i,j)**2 + yp(i,j)**2)
            if (r .gt. 1.d-2*dxc) then
                theta(i,j) = dacos(xp(i,j)/r)
              elseif (yp(i,j).gt.0.d0) then
                theta(i,j) = 0.d0
              endif
            if (yp(i,j) .lt. 0.d0) then
                theta(i,j) = -theta(i,j)
                endif

c           # compute phi, angle down from north pole:
            if (zp(i,j) .gt. 0.d0) then
                phi(i,j) = pi/2.d0 - dacos(r/Rsphere)  
              else
                phi(i,j) = pi/2.d0 + dacos(r/Rsphere)  
              endif
   15	    continue

c
      areamax = 0.d0
      areamin = 1.d3
      velmax = 0.d0
      do 20 j=1-mbc,my+mbc
         do 20 i=1-mbc,mx+mbc
c
c           # compute edge velocities by differencing stream function:
c
	    aux(i,j,1) = (stream(xp(i,j+1),yp(i,j+1),zp(i,j+1))
     &			  - stream(xp(i,j),yp(i,j),zp(i,j))) / dyc
c
	    aux(i,j,2) = -(stream(xp(i+1,j),yp(i+1,j),zp(i+1,j)) 
     & 			   - stream(xp(i,j),yp(i,j),zp(i,j))) / dxc
            if (dmax1(dabs(aux(i,j,1)), dabs(aux(i,j,2))) .gt. velmax)
     &        then
                velmax = dmax1(dabs(aux(i,j,1)), dabs(aux(i,j,2)))
                ivmax = i
                jvmax = j
              endif

c           # compute area of physical cell from four corners:
c           # find area on the sphere of two spherical triangles obtained
c           # by subdividing rectangle.  See
c           #    http://mathforum.org/library/drmath/view/65316.html

c           # corners are labeled  
c           #      1: (i,j)      2: (i+1,j)    3: (i,j+1)    4: (i+1,j+1)

            beta12 = dsin(phi(i,j))*dsin(phi(i+1,j))*
     &                  dcos(theta(i,j) - theta(i+1,j))
     &                + dcos(phi(i,j))*dcos(phi(i+1,j))
            beta23 = dsin(phi(i,j+1))*dsin(phi(i+1,j))*
     &                  dcos(theta(i,j+1) - theta(i+1,j))
     &                 + dcos(phi(i,j+1))*dcos(phi(i+1,j))
            beta13 = dsin(phi(i,j+1))*dsin(phi(i,j))*
     &                  dcos(theta(i,j+1) - theta(i,j))
     &                 + dcos(phi(i,j+1))*dcos(phi(i,j))
            beta24 = dsin(phi(i+1,j+1))*dsin(phi(i+1,j))*
     &                  dcos(theta(i+1,j+1) - theta(i+1,j))
     &                 + dcos(phi(i+1,j+1))*dcos(phi(i+1,j))
            beta34 = dsin(phi(i+1,j+1))*dsin(phi(i,j+1))*
     &                  dcos(theta(i+1,j+1) - theta(i,j+1))
     &                 + dcos(phi(i+1,j+1))*dcos(phi(i,j+1))

c           # great circles distances between corners:
            d12 = Rsphere * dacos(beta12)
            d23 = Rsphere * dacos(beta23)
            d13 = Rsphere * dacos(beta13)
            d24 = Rsphere * dacos(beta24)
            d34 = Rsphere * dacos(beta34)

            s123 = 0.5d0 * (d12 + d23 + d13)
            s234 = 0.5d0 * (d23 + d34 + d24)

c           # spherical excess for each triangle:
            t123 = dtan(s123/2.d0)*dtan((s123-d12)/2.d0)
     &             *dtan((s123-d23)/2.d0)*dtan((s123-d13)/2.d0)
            t123 = dmax1(t123,0.d0)
            E123 = 4.d0*datan(sqrt(t123))

            t234 = dtan(s234/2.d0)*dtan((s234-d23)/2.d0)
     &             *dtan((s234-d34)/2.d0)*dtan((s234-d24)/2.d0)
            t234 = dmax1(t234,0.d0)
            E234 = 4.d0*datan(sqrt(t234))

            area = Rsphere**2 * (E123 + E234)
c
c           # capacity kappa:
	    aux(i,j,3) = area / (dxc*dyc)
            if (i.gt.0 .and. i.le.mx .and. j.gt.0 .and. j.le.my) then
               areamax = dmax1(areamax,area)
               areamin = dmin1(areamin,area)
               endif
c           write(27,*) i,j,area,aux(i,j,3)
            if (aux(i,j,3) .lt. 0.5d0 .or. aux(i,j,3).gt.3.d0) then
               write(27,*) '   '
               write(27,*) i,j,area,aux(i,j,3)
               write(27,*) '   ',mx,my,xlower,ylower,dxc,dyc
               write(27,*) '   ',xc(i),yc(j),theta(i,j),phi(i,j)
               write(27,*) '   ',xp(i,j),yp(i,j),zp(i,j)
               endif

            write(21,210) i,j,xc(i),yc(j),xp(i,j),yp(i,j),zp(i,j),
     &                    theta(i,j),phi(i,j),t123,t234,
     &                    aux(i,j,1),aux(i,j,2),aux(i,j,3)
  210       format(2i3,5e14.5,/,'      ',4e14.5,/,'      ',3e14.5,/)
c
   20       continue
c
c      write(27,*) areamax,areamin
c      write(6,*) 'area ratio = ',areamax/areamin
c      write(6,*) 'max velocity = ',velmax, '  at ', ivmax,jvmax
       return
       end
