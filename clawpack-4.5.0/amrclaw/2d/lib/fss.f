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
