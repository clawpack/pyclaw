c
c -----------------------------------------------------------
c
      subroutine colate (badpts, len, lcheck, 
     1                   iflags,domflags,isize,jsize,npts)
c
      implicit  double precision (a-h,o-z)
      dimension badpts(2,len)
      integer*1   iflags  (0:isize+1,0:jsize+1)
      integer*1   domflags(0:isize+1,0:jsize+1)

      include  "call.i"
c
c
c *************************************************************
c
c colate = takes the error plane with flagged pts at level lcheck
c          and puts their (i,j) cell centered
c          indices into the badpts array.
c          To insure proper nesting,  get rid of flagged point
c          that don't fit into properly nested domain (in iflags2)
c
c *************************************************************
c
c     # if pt. flagged but domain not flagged, turn it off
c     # note that this results in flags of 1,  not 2 of 3.


      if (dprint) then
         write(outunit,*)" from colate: iflags"
         do 48 jj = 1, jsize
           j = jsize + 1 - jj
           write(outunit,101)(iflags(i,j),i=1,isize)
 48      continue
         write(outunit,*)" from colate: domflags"
         do 49 jj = 1, jsize
           j = jsize + 1 - jj
           write(outunit,101)(domflags(i,j),i=1,isize)
 101       format(80i1)
 49      continue
      endif

      do 10 j = 1, jsize 
      do 10 i = 1, isize 
        iflags(i,j) = min(iflags(i,j),domflags(i,j))
 10   continue


      index  = 0
c
c     give points the indices from integer region space.
      do 20 j   = 1, jsize
      do 20 i   = 1, isize
        if (iflags(i,j) .ne. goodpt) then
          index = index + 1
          badpts(1,index) = dble(i)-.5
          badpts(2,index) = dble(j)-.5
        endif
 20   continue
c
 99   npts = index 
      if (gprint) then
        write(outunit,100) npts, lcheck
 100    format( i5,' flagged points colated on level ',i4)
      endif

      return
      end
