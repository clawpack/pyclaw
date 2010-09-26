c
c
c =========================================================
      subroutine copyq2(maxmx,maxmy,meqn,mbc,mx,my,q1,q2)
c =========================================================
c
c     # copy the contents of q1 into q2
c
      implicit double precision (a-h,o-z)
      dimension q1(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension q2(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
c
c
      do 10 j = 1-mbc, my+mbc
        do 10 i = 1-mbc, mx+mbc
          do 10 m=1,meqn
             q2(i,j,m) = q1(i,j,m)
   10        continue
      return
      end
