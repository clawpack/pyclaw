c
c ----------------------------------------------------
c
      subroutine domshrink(iflags2,iflags,idim,jdim)

      implicit double precision (a-h, o-z)

      include  "call.i"

      integer*1  iflags2(0:idim+1,0:jdim+1)
      integer*1  iflags (0:idim+1,0:jdim+1)

c
c :::::::::::::::::::::::::  DOMSHRINK ::::::::::::::::::::::::::::
c
c  shrink domain flags one cell for allowable properly nested domain
c  This is needed even for lcheck = lbase. More shrinking needed
c  for finer levels.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (dprint) then
        write(outunit,*)" from domshrink: on entry, iflags2"
        do 10 jj = 1, jdim
        j = jdim + 1 - jj
        write(outunit,100)(iflags2(i,j),i=1,idim)
 100    format(80i1)
 10     continue
      endif

      do 40 j = 1, jdim
      do 40 i = 1, idim
         iflags(i,j) = iflags2(i,j)
         if (iflags2(i  ,j  ) .le. 0 .or.
     1       iflags2(i+1,j  ) .le. 0 .or. iflags2(i-1,j  ) .le. 0 .or. 
     2       iflags2(i+1,j+1) .le. 0 .or. iflags2(i-1,j+1) .le. 0 .or. 
     3       iflags2(i  ,j-1) .le. 0 .or. iflags2(i  ,j+1) .le. 0 .or.
     4       iflags2(i+1,j-1) .le. 0 .or. iflags2(i-1,j-1) .le. 0) then
                 iflags(i,j) = 0
          endif
 40   continue
c
c if border of domain touches a physical boundary then set domain in
c ghost cell as well
c
C WHY DOESNT THIS HAVE TO HAVE PERIODIC CLAUSE(AND spheredom)
       if (.not. xperdom) then
         do 55 j = 1, jdim
           if (iflags(1,j) .eq. 1) iflags(0,j) = 1
           if (iflags(idim,j) .eq. 1) iflags(idim+1,j) = 1
 55      continue
       endif
       if (.not. yperdom) then
         do 65 i = 1, idim
           if (iflags(i,1) .eq. 1) iflags(i,0) = 1
           if (iflags(i,jdim) .eq. 1) iflags(i,jdim+1) = 1
 65      continue
       endif

      if (dprint) then
        write(outunit,*)" from domshrink: on exit, iflags"
        do 70 jj = 1, jdim
        j = jdim + 1 - jj
        write(outunit,100)(iflags(i,j),i=1,idim)
 70     continue
      endif

      return
      end
