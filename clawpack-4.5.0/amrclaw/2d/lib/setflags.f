c
c --------------------------------------------------------------
c
      subroutine setflags(iflags,isize,jsize,
     1                    rctold,idim3,mitot,mjtot,mptr)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      dimension rctold(mitot,mjtot,idim3)
      integer*1 iflags(0:isize+1,0:jsize+1)

c :::::::::::::::::::::: SETFLAGS ::::::::::::::::::::::::::::::::::
c transfer flagged arrays into 1 large array of entire domain
c makes buffering, projecting, etc. easier without searching 
c through all kinds of grids
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      ibeg = node(ndilo,mptr) - nghost
      jbeg = node(ndjlo,mptr) - nghost

      do 10 j = nghost+1, mjtot-nghost
      do 10 i = nghost+1, mitot-nghost
        iflags(ibeg+i,jbeg+j) = iflags(ibeg+i,jbeg+j) + rctold(i,j,1)
 10   continue
c
 99   return
      end
