c
c
c
c
c     =================================================
      subroutine cellave(xlow,ylow,dx,dy,wl)
c     =================================================
      implicit double precision (a-h,o-z)
      external fss
      logical fl(5),alll,allr
      dimension x(10),y(10),xx(5),yy(5)
      common/fsscorn/ xc0,yc0,xc1,yc1
c   
c     # compute wl, fraction of cell that lies in left state.
c     # For initial data with two states ql and qr separated by a 
c     # discontinuity. The curve along which the discontinuity lies is
c     # specified by the function fdisc, which should return a value that
c     # is negative on the side where ql lies and positive on the qr side.
c
c     # xlow,ylow is the coordinate of the lower left corner of the cell.
c     # dx, dy are grid spacing in x and y.
c
      xx(1) = xlow
      xx(2) = xlow
      xx(3) = xlow+dx
      xx(4) = xlow+dx
      xx(5) = xx(1)
      yy(1) = ylow
      yy(2) = ylow+dy
      yy(3) = ylow+dy
      yy(4) = ylow
      yy(5) = yy(1)
      alll = .true.
      allr = .true.
c
      do 20 i=1,4
         fl(i) = fdisc(xx(i),yy(i)) .lt. 0.d0
         alll = alll .and. fl(i)
         allr = allr .and. (.not. fl(i))
   20    continue
      fl(5) = fl(1)
c
      if (alll) then
         wl = 1.d0
         return
         endif
      if (allr) then
         wl = 0.d0
         return
         endif
c
      iv = 0
      do 40 i=1,4
          if (fl(i)) then
               iv = iv+1
               x(iv) = xx(i)
               y(iv) = yy(i)
               endif
          if (fl(i).neqv.fl(i+1)) then
               iv = iv+1
               xc0 = xx(i)
               yc0 = yy(i)
               xc1 = xx(i+1)
               yc1 = yy(i+1)
               ss = zeroin(0.d0, 1.d0, fss, 1d-8)
c              write(27,*) 'xc,yc,ss:',xc0,yc0,xc1,yc1,ss
               x(iv) = xx(i) + ss*(xx(i+1)-xx(i))
               y(iv) = yy(i) + ss*(yy(i+1)-yy(i))
               endif
   40     continue
c
c     # compute area:
c
      if (iv.eq.0) then
         wl = 0.d0
         return
         endif
c
      x(iv+1) = x(1)
      y(iv+1) = y(1)
      area = 0.d0
      do 50 i=1,iv
         area = area + .5d0*(y(i)+y(i+1))*(x(i+1)-x(i))
c        write(27,*) '  x,y:',x(i),y(i)
   50    continue
c
      wl = area / (dx*dy)
c     write(27,*) 'area,wl:',area,wl
c
      return
      end
