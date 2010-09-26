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
c
c
c
c
c
c     =================================================
      function fss(s)
c     =================================================
      implicit double precision (a-h,o-z)
      common/fsscorn/ xc0,yc0,xc1,yc1
c   
c     # compute fdisc at distance s between corners (xc0,yc0) and (xc1,yc1)
c
      x = xc0 + s*(xc1-xc0)
      y = yc0 + s*(yc1-yc0)
      fss = fdisc(x,y)
      return
      end
c
c
c
c     =================================================
      function zeroin(ax,bx,f,tol)                                         
c     =================================================
      implicit double precision (a-h,o-z)
      external f
c         
c      a zero of the function  f(x)  is computed in the interval ax,bx .        
c      (Standard routine from netlib)
c         
c  input..
c         
c  ax     left endpoint of initial interval       
c  bx     right endpoint of initial interval      
c  f      function subprogram which evaluates f(x) for any x in       
c         the interval  ax,bx 
c  tol    desired length of the interval of uncertainty of the        
c         final result ( .ge. 0.d0)      
c         
c         
c  output..         
c         
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx   
c         
c         
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval 
c  ax,bx  to within a tolerance  4*macheps*dabs(x) + tol, where macheps
c  is the relative machine precision.   
c
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).    
c         
c         
c         
c  compute eps, the relative machine precision    
c         
      eps = 1.d0     
   10 eps = eps/2.d0 
      tol1 = 1.d0 + eps        
      if (tol1 .gt. 1.d0) go to 10       
c         
c initialization    
c         
      a = ax        
      b = bx        
      fa = f(a)     
      fb = f(b)     
c         
c begin step        
c         
   20 c = a         
      fc = fa       
      d = b - a     
      e = d         
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b         
      b = c         
      c = a         
      fa = fb       
      fb = fc       
      fc = fa       
c         
c convergence test  
c         
   40 tol1 = 2.d0*eps*dabs(b) + 0.5*tol   
      xm = .5*(c - b)         
      if (dabs(xm) .le. tol1) go to 90   
      if (fb .eq. 0.d0) go to 90         
c         
c is bisection necessary      
c         
      if (dabs(e) .lt. tol1) go to 70    
      if (dabs(fa) .le. dabs(fb)) go to 70
c         
c is quadratic interpolation possible   
c         
      if (a .ne. c) go to 50  
c         
c linear interpolation        
c         
      s = fb/fa     
      p = 2.d0*xm*s  
      q = 1.d0 - s   
      go to 60      
c         
c inverse quadratic interpolation       
c         
   50 q = fa/fc     
      r = fb/fc     
      s = fb/fa     
      p = s*(2.d0*xm*q*(q - r) - (b - a)*(r - 1.d0))
      q = (q - 1.d0)*(r - 1.d0)*(s - 1.d0) 
c         
c adjust signs      
c         
   60 if (p .gt. 0.d0) q = -q  
      p = dabs(p)    
c         
c is interpolation acceptable 
c         
      if ((2.d0*p) .ge. (3.d0*xm*q - dabs(tol1*q))) go to 70   
      if (p .ge. dabs(0.5*e*q)) go to 70 
      e = d         
      d = p/q       
      go to 80      
c         
c bisection         
c         
   70 d = xm        
      e = d         
c         
c complete step     
c         
   80 a = b         
      fa = fb       
      if (dabs(d) .gt. tol1) b = b + d   
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)     
      if ((fb*(fc/dabs(fc))) .gt. 0.d0) go to 20    
      go to 30      
c         
c done    
c         
   90 zeroin = b    
      return        
      end 
