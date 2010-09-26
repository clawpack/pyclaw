c
c
c
c
c =========================================================
      subroutine dump2(maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &           q,t,dx,dy,dtv,cflv,nv,method,mthlim,info)
c =========================================================
c
c     # Dump the data needed for a possible restart
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension dtv(5),cflv(4),nv(2),method(7),mthlim(mwaves)
c
      write(99,*) mx,my
      write(99,*) dx,dy
      write(99,*) t
      write(99,*) dtv
      write(99,*) cflv
      write(99,*) nv
      write(99,*) method
      write(99,*) mthlim
      write(99,*) info
c
      do 12 m=1,meqn
         do 11 j=1,my
            do 10 i=1,mx
               write(99,*) q(i,j,m)
   10       continue
   11    continue
   12 continue
      return
      end
