c
c
c =========================================================
      subroutine copyq1(maxmx,meqn,mbc,mx,q1,q2)
c =========================================================
c
c     # copy the contents of q1 into q2
c
      implicit double precision (a-h,o-z)
      dimension q1(1-mbc:maxmx+mbc, meqn)
      dimension q2(1-mbc:maxmx+mbc, meqn)
c
      do 10 i = 1-mbc, mx+mbc
          do 10 m=1,meqn
             q2(i,m) = q1(i,m)
   10        continue
      return
      end
