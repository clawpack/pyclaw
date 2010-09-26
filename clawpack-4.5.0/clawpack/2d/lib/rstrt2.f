c
c
c =========================================================
      subroutine rstrt2(maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &           q,t,dx,dy,dtv,cflv,nv,method,mthlim,info)
c =========================================================
c
c     # Read in  the data needed for a restart
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension dtv(5),cflv(4),nv(2),method(7),mthlim(mwaves)
c
      read(99,*) mx,my
      read(99,*) dx,dy
      read(99,*) t
      read(99,*) dtv
      read(99,*) cflv
      read(99,*) nv
      read(99,*) method
      read(99,*) mthlim
      read(99,*) info
c
      do 12 m=1,meqn
         do 11 j=1,my
            do 10 i=1,mx
               read(99,*) q(i,j,m)
   10       continue
   11    continue
   12 continue
      return
      end
