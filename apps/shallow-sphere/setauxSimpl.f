c     ============================================
      subroutine setaux(mx,my,xlower,ylower,dxc,dyc,
     &                  maux,aux,Rsphere)
c     ============================================
c
c     # set auxiliary arrays for shallow water equations on the sphere
c
c     # on input, (xc(i),yc(j)) gives uniformly spaced computational grid.
c     # on output, 
c     The aux array has the following elements:
c         1  kappa = ratio of cell area to dxc*dyc
c         2  enx = x-component of normal vector to left edge in tangent plane
c         3  eny = y-component of normal vector to left edge in tangent plane
c         4  enz = z-component of normal vector to left edge in tangent plane
c         5  etx = x-component of tangent vector to left edge in tangent plane
c         6  ety = y-component of tangent vector to left edge in tangent plane
c         7  etz = z-component of tangent vector to left edge in tangent plane
c         8  enx = x-component of normal vector to bottom edge in tangent plane
c         9  eny = y-component of normal vector to bottom edge in tangent plane
c        10  enz = z-component of normal vector to bottom edge in tangent plane
c        11  etx = x-component of tangent vector to bottom edge in tangent plane
c        12  ety = y-component of tangent vector to bottom edge in tangent plane
c        13  etz = z-component of tangent vector to bottom edge in tangent plane
c        14  erx = x-component of unit vector in radial direction at cell ctr
c        15  ery = y-component of unit vector in radial direction at cell ctr
c        16  erz = z-component of unit vector in radial direction at cell ctr
c
c     
      implicit double precision (a-h,o-z)
      dimension xc(1:mx+1), yc(1:my+1)
      dimension xp(1:mx+1,1:my+1), yp(1:mx+1,1:my+1)
      dimension zp(1:mx+1,1:my+1)
      dimension theta(1:mx+1,1:my+1), phi(1:mx+1,1:my+1)
      dimension aux(maux,1:mx,1:my)
cf2py integer optional,intent(in) mx
cf2py integer optional,intent(in) my
cf2py double precision intent(in) xlower
cf2py double precision intent(in) ylower
cf2py double precision intent(in) dx
cf2py double precision intent(in) dy
cf2py integer optional, intent(in)  maux
cf2py intent(in,out) aux
cf2py double precision intent(in) Rsphere

      pi = 4.d0*datan(1.d0)
c
c      if (mbc .gt. 4) then
c	 write(6,*)'***  increase size of local arrays in setaux ***'
c	 stop
c	 endif
c      if (mx+mbc+1.gt.1005 .or. my+mbc+1.gt.1005) then
c	 write(6,*)'***  increase size of 1005 in setaux ***'
c	 stop
c	 endif
c
c     # Set xc and yc so that (xc(i),yc(j)) is the 
c     # lower left corner of (i,j) cell in computational space:
c
      do 10 i=1,mx+1
         xc(i) = xlower + (i-1.d0) * dxc
   10    continue
c
      do 12 j=1,my+1
         yc(j) = ylower + (j-1.d0) * dyc
   12    continue

c     # compute cell corners on sphere and angles phi, theta
c     # related to latitude and longitude
c
      do 15 j=1,my+1
         do 15 i=1,mx+1

c           # map computational point to (xp,yp,zp) on sphere:
            call mapc2m(xc(i),yc(j),xp(i,j),yp(i,j),zp(i,j),Rsphere)

c           # compute longitude theta from positive x axis:
            r = dsqrt(xp(i,j)**2 + yp(i,j)**2)

c            theta(i,j) = atan2(xp(i,j),yp(i,j))

             if (r.gt.1.d-4) then
               theta(i,j) = dacos(xp(i,j)/r)
             elseif (yp(i,j).gt.0.d0) then
               theta(i,j) = 0.d0
             endif
             if (yp(i,j).lt.0.d0) then
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
      do 20 j=1,my
         do 20 i=1,mx
c
c           # compute normal and tangent vectors to left edge (in tangent plane)
c
c           # tangent vector is edge vector:
            etx = xp(i,j+1) - xp(i,j)
            ety = yp(i,j+1) - yp(i,j)
            etz = zp(i,j+1) - zp(i,j)
            aux(5,i,j) = etx
            aux(6,i,j) = ety
            aux(7,i,j) = etz

c           # normal to sphere in radial direction at midpoint of edge:
            erx = 0.5d0*(xp(i,j) + xp(i,j+1))
            ery = 0.5d0*(yp(i,j) + yp(i,j+1))
            erz = 0.5d0*(zp(i,j) + zp(i,j+1))
c
c           # normal to edge in tangent plane is cross product of et and er:
            enx = ety*erz - etz*ery
            eny = etz*erx - etx*erz
            enz = etx*ery - ety*erx
            ennorm = dsqrt(enx**2 + eny**2 + enz**2)
            aux(2,i,j) = enx/ennorm
            aux(3,i,j) = eny/ennorm
            aux(4,i,j) = enz/ennorm


c           # compute normal and tangent vectors to bottom edge (in tang. pl.)
c
c           # tangent vector is edge vector:
            etx = xp(i+1,j) - xp(i,j)
            ety = yp(i+1,j) - yp(i,j)
            etz = zp(i+1,j) - zp(i,j)
            aux(11,i,j) = etx
            aux(12,i,j) = ety
            aux(13,i,j) = etz

c           # normal to sphere in radial direction at midpoint of edge:
            erx = 0.5d0*(xp(i,j) + xp(i+1,j))
            ery = 0.5d0*(yp(i,j) + yp(i+1,j))
            erz = 0.5d0*(zp(i,j) + zp(i+1,j))
c
c           # normal to edge in tangent plane is cross product of er and et:
            enx = ery*etz - erz*ety
            eny = erz*etx - erx*etz
            enz = erx*ety - ery*etx
            ennorm = dsqrt(enx**2 + eny**2 + enz**2)
            aux(8,i,j) = enx/ennorm
            aux(9,i,j) = eny/ennorm
            aux(10,i,j) = enz/ennorm

c           # normal to sphere in radial direction at cell center:
            xcm = xlower+(i-0.5)*dxc
            ycm = ylower+(j-0.5)*dyc
            call mapc2m(xcm,ycm,xpm,ypm,zpm,Rsphere)

            aux(14,i,j) = xpm
            aux(15,i,j) = ypm
            aux(16,i,j) = zpm

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

c           area = Rsphere**2 * (E123 + E234)
            area = (E123 + E234)
c
c           # capacity kappa:
            aux(1,i,j) = area / (dxc*dyc)
           
   20       continue


       return


       end

