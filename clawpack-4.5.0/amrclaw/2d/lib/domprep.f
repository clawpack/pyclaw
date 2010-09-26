c
c ----------------------------------------------------
c
      subroutine domprep(domflags,lbase,ibase,jbase)

      implicit double precision (a-h, o-z)

      include  "call.i"

      integer*1  domflags(0:ibase+1,0:jbase+1)

c
c ::::::::::::::::::::::::::: PREPDOM :::::::::::::::::::::
c 
c  prepare 2 dimensional array of domain for proper nesting
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::


      do 10 j = 0, jbase+1
      do 10 i = 0, ibase+1
         domflags(i,j) = 0
 10   continue

      mptr = lstart(lbase)
 15   continue
c     domain flags appears to be 1 based indexing, so 0 a ghost cell.
c     should change it to be 0 based, like grids, so border is at -1.
      do 20 j = node(ndjlo,mptr) + 1, node(ndjhi,mptr) + 1
      do 20 i = node(ndilo,mptr) + 1, node(ndihi,mptr) + 1
         domflags(i,j) = 1
 20   continue
      mptr = node(levelptr, mptr)
      if (mptr .ne. 0) go to 15

c
c take care of periodic domains or if border of domain touches a 
c  physical boundary then set domain in ghost cell as well
c
      if (xperdom) then
         do 25 j = 0, jbase+1
           domflags(0,j)       = domflags(ibase,j)
           domflags(ibase+1,j) = domflags(1,j)
 25      continue
       else
         do 65 j = 1, jbase
           domflags(0,j)       = domflags(1,j)
           domflags(ibase+1,j) = domflags(ibase,j)
 65      continue
      endif
      if (yperdom) then
         do 35 i = 0, ibase+1
           domflags(i,0)       = domflags(i,jbase)
           domflags(i,jbase+1) = domflags(i,1)
 35      continue
       else if (spheredom) then
         do 36 i = 0, ibase+1
           domflags(i,0)       = domflags(ibase+1-i,1)
           domflags(i,jbase+1) = domflags(ibase+1-i,jbase)
 36     continue
       else
         do 55 i = 1, ibase
           domflags(i,0)       = domflags(i,1)
           domflags(i,jbase+1) = domflags(i,jbase)
 55      continue
      endif
c
c the 4 corners
c
        if (domflags(0,1)+domflags(1,0) .eq. 2) domflags(0,0)=1
        if (domflags(ibase,0)+domflags(ibase+1,1) .eq. 2) 
     .          domflags(ibase+1,0)=1
        if (domflags(ibase,jbase+1)+domflags(ibase+1,jbase) .eq. 2) 
     .          domflags(ibase+1,jbase+1)=1
        if (domflags(0,jbase)+domflags(1,jbase+1) .eq. 2) 
     .          domflags(0,jbase+1)=1

      if (dprint) then
         write(outunit,*)" from domprep: domflags at level  ", lbase
         do 40 jj = 1, jbase
         j = jbase + 1 - jj
         write(outunit,100)(domflags(i,j),i=1,ibase)
 100     format(80i1)
 40      continue
      endif

      return
      end
