c
c --------------------------------------------------------------------
c
       subroutine intcopy(val,mitot,mjtot,nvar,ilo,ihi,jlo,jhi,level,
     &                    iputst,jputst)

       implicit double precision (a-h, o-z)

       dimension val(mitot,mjtot,nvar)

       include "call.i"

       iadd(i,j,ivar) = loc + i - 1 + mi*((ivar-1)*mj+j-1)

c ::::::::::::::::::::::::::: INTCOPY :::::::::::::::::::::::::::::::
c
c    find intersecting grids at the same level. copy data from
c    old grid to val. 
c    old grid has "nghost" ghost cells - passed in nodal common block.
c    new grid has no ghost cells - indices describe entire patch.
c    iputst, jputst: where to copy values into. may not be in
c                    location corresponding to ilo,ihi,etc. if
c                    the patch has been periodically wrapped.
c    
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


       mptr = lstart(level)

 10    if (mptr .eq. 0) go to 99
          iglo = node(ndilo,mptr)
          ighi = node(ndihi,mptr)
          jglo = node(ndjlo,mptr)
          jghi = node(ndjhi,mptr)

c         # does it intersect?
          ixlo = max(iglo,ilo)
          ixhi = min(ighi,ihi)
          jxlo = max(jglo,jlo)
          jxhi = min(jghi,jhi)

          if (ixlo .le. ixhi .and. jxlo .le. jxhi) then
              loc  = node(store1,mptr)
              nx   = ighi - iglo + 1
              ny   = jghi - jglo + 1
              mi   = nx + 2*nghost
              mj   = ny + 2*nghost
              do 20 j    = jxlo, jxhi
              do 20 ivar = 1, nvar
              do 30 i    = ixlo, ixhi
                  val(iputst+i-ilo,jputst+j-jlo,ivar) =
     1                alloc(iadd(i-iglo+nghost+1,j-jglo+nghost+1,ivar))
 30           continue
 20           continue
          endif
          mptr = node(levelptr, mptr)
          go to 10

 99   return
      end
