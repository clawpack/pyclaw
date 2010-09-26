c
c ----------------------------------------------------
c
      subroutine domcopy(iflags2,iflags,isize,jsize)

      implicit double precision (a-h, o-z)

      include  "call.i"

      integer*1  iflags2(0:isize+1,0:jsize+1)
      integer*1  iflags (0:isize+1,0:jsize+1)

c
c ::::::::::::::::::::::::::: DOMCOPY :::::::::::::::::::::
c 
c  domain flags are in iflags. copy into iflags2.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::


      do 10 j = 0, jsize+1
      do 10 i = 0, isize+1
         iflags2(i,j) = iflags(i,j)
 10   continue
c
c  take care of periodicity again
c
      if (xperdom) then
         do 35 j = 0, jsize+1
           iflags2(0,j)       = iflags2(isize,j)
           iflags2(isize+1,j) = iflags2(1,j)
 35      continue
       else
         do 55 j = 1, jsize
           if (iflags2(1,j) .eq. 1) iflags2(0,j) = 1
           if (iflags2(isize,j) .eq. 1) iflags2(isize+1,j) = 1
 55      continue
      endif
      if (yperdom) then
         do 45 i = 0, isize+1
           iflags2(i,0)       = iflags2(i,jsize)
           iflags2(i,jsize+1) = iflags2(i,1)
 45      continue
       else if (spheredom) then
         do 46 i = 0, isize+1
         iflags2(i,0)       = iflags2(isize+1-i,1)
         iflags2(i,jsize+1) = iflags2(isize+1-i,jsize)
 46      continue
       else
         do 65 i = 1, isize
           if (iflags2(i,1) .eq. 1) iflags2(i,0) = 1
           if (iflags2(i,jsize) .eq. 1) iflags2(i,jsize+1) = 1
 65      continue
      endif

      if (dprint) then
         write(outunit,*)" from domcopy: domflags "
         do 40 jj = 1, jsize
         j = jsize + 1 - jj
         write(outunit,100)(iflags2(i,j),i=1,isize)
 100     format(80i1)
 40      continue
      endif

 
      return
      end
