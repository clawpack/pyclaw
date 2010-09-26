c
c ---------------------------------------------------------
c
      subroutine outmsh(mptr,outgrd,nvar,naux)
      implicit double precision (a-h,o-z)
      logical  outgrd

      include  "call.i"

c
c ::::::::::::::::::::: OUTMSH :::::::::::::::::::::::::::::::::::::
c
c outmsh - output the grid descriptor and optionally the values on
c          the grid (for a single grid - see "outtre" for outputing
c          a subtree)
c input parameters:
c    mptr   - ptr to grid descriptor
c    outgrd - if true, output value on grid
c special case
c     if grid has level < 1, nothing is printed. (forest pointer
c has level = 0).  this simplifies logic of outtre; also, any grid
c with level <= 0 is not properly set (the forest pointer
c is used only to provide pointers into the tree of coarsest grids)
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
100   format(1x,47h+----------------------------------------------,
     *30h-----------------------------+)
c
      lev = node(nestlevel,mptr)
c
      write(outunit,100)
      write(outunit,101) mptr
101   format(1x,10h! grid no:,i4,62x,1h!)
      write(outunit,102) node(nestlevel,mptr),rnode(timemult,mptr),
     .             node(levelptr,mptr)
102   format(1x,1h!,11h nestlevel=,i3,12h, time mult=,f8.5,
     1       13h, level ptr =,i4,24x,1h!)
      write(outunit,103) node(store1,mptr),node(store2,mptr),
     1                   node(cfluxptr,mptr),node(ffluxptr,mptr)
 103  format(1x,'! storage locs =',2i8,'  bndry locs =',2i8,14x,1h!)
      write(outunit,104)
      write(outunit,111) rnode(cornxlo,mptr),rnode(cornyhi,mptr),
     1             rnode(cornxhi,mptr),rnode(cornyhi,mptr)
      write(outunit,111) rnode(cornxlo,mptr),rnode(cornylo,mptr),
     1             rnode(cornxhi,mptr),rnode(cornylo,mptr)
      write(outunit,112)
      write(outunit,113) node(ndilo,mptr),node(ndjhi,mptr),
     1                   node(ndihi,mptr),node(ndjhi,mptr)
      write(outunit,113) node(ndilo,mptr),node(ndjlo,mptr),
     1                   node(ndihi,mptr),node(ndjlo,mptr)
112   format(1x,23h! integer index space :,53x,1h!)
113   format(1x,2h! ,18x,2(1h(,i8,2h, ,i8,1h)),16x,1h!)
104   format(1x,23h! corners of rectangle:,53x,1h!)
111   format(1x,2h! ,18x,2(1h(,f8.3,2h, ,f8.3,1h)),16x,1h!)
      write(outunit,105) hxposs(lev),hyposs(lev),possk(lev)
105   format(1x,7h! hrow=,f9.6,7h, hcol=,f9.6,8h, ktime=,f9.6,27x,1h!)
      write(outunit,100)
c
      if (.not. outgrd) go to 99
c output the grid
      mitot   = node(ndihi,mptr) - node(ndilo,mptr) + 1 + 2*nghost
      mjtot   = node(ndjhi,mptr) - node(ndjlo,mptr) + 1 + 2*nghost
      loc     = node(store1,mptr)
      locaux  = node(storeaux,mptr)
      call outval(alloc(loc),nvar,mitot,mjtot,mptr,outgrd,
     1            naux,alloc(locaux))
 99   return
      end
