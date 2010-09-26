c
c --------------------------------------------------
c
      subroutine birect(mptr1)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

c
c :::::::::::::  BIRECT :::::::::::::::::::::::::::::::::::::::
c check each grid, starting with mptr1 (either newstl or lstart)
c to see that it has no more than max1d points in either dimensions.
c needed so that scratch array space in stepgrid not exceeded.
c
c also check for too small grids - but has never happened.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      mptr  = mptr1
      level = node(nestlevel,mptr)
      hx    = hxposs(level)
      hy    = hyposs(level)
c
10    continue
      cxlo    = rnode(cornxlo,mptr)
      cxhi    = rnode(cornxhi,mptr)
      cylo    = rnode(cornylo,mptr)
      cyhi    = rnode(cornyhi,mptr)
      nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
      ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
      minsize = 2*nghost
c
c check number of rows first - if too many, bisect grid with vertical
c line down the middle. make sure new grid corners are anchored
c on coarse grid point. make sure if bisecting coarse grid that
c new grids have even number of points
c
      if (nx + 2*nghost .gt. max1d) then
 
        nxl    = nx/2
        if (level .gt. 1) then 
           lratio = intratx(level-1)
        else
           lratio = 2
        endif
        nxl = (nxl/lratio)*lratio 
        nxr    = nx - nxl 
        cxmid  = cxlo + nxl*hx
 
        mptrnx = nodget(dummy)
        node(levelptr,mptrnx) = node(levelptr,mptr)
        node(levelptr,mptr)   = mptrnx
 
        rnode(cornxhi,mptr) = cxmid
        node(ndihi,mptrnx)  = node(ndihi,mptr)
        node(ndihi,mptr)    = node(ndilo,mptr) + nxl - 1
        node(ndilo,mptrnx)  = node(ndihi,mptr) + 1
        node(ndjhi,mptrnx)  = node(ndjhi,mptr)
        node(ndjlo,mptrnx)  = node(ndjlo,mptr)

        rnode(cornxlo,mptrnx)    = cxmid
        rnode(cornylo,mptrnx)    = cylo
        rnode(cornyhi,mptrnx)    = cyhi
        rnode(cornxhi,mptrnx)    = cxhi
        rnode(timemult,mptrnx)   = rnode(timemult,mptr)
        node(nestlevel,mptrnx)   = node(nestlevel,mptr)

        go to 10
c
c check number of columns next - if too many, bisect grid with horizontal
c line down the middle
c
      else if (ny + 2*nghost .gt. max1d) then
 
        nyl    = ny/2
        if (level .gt. 1) then 
           lratio = intraty(level-1)
        else
           lratio = 2
        endif
        nyl = (nyl/lratio)*lratio
        nyr    = ny - nyl 
        cymid  =  cylo + nyl*hy
 
        mptrnx = nodget(dummy)
        node(levelptr,mptrnx) = node(levelptr,mptr)
        node(levelptr,mptr)   = mptrnx
 
        rnode(cornyhi,mptr)   = cymid

        node(ndjhi,mptrnx) = node(ndjhi,mptr)
        node(ndjhi,mptr)   = node(ndjlo,mptr) + nyl - 1
        node(ndjlo,mptrnx) = node(ndjhi,mptr) + 1
        node(ndihi,mptrnx) = node(ndihi,mptr)
        node(ndilo,mptrnx) = node(ndilo,mptr)

        rnode(cornxlo,mptrnx)   = cxlo
        rnode(cornylo,mptrnx)   = cymid
        rnode(cornyhi,mptrnx)   = cyhi
        rnode(cornxhi,mptrnx)   = cxhi
        node(nestlevel,mptrnx)  = node(nestlevel,mptr)
        rnode(timemult,mptrnx)  = rnode(timemult,mptr)
        go to 10
c
c  grid ok - check the next
c
      else
        mptr = node(levelptr,mptr)
        if (mptr.ne.0) go to 10
 
      endif
 
      return
      end
