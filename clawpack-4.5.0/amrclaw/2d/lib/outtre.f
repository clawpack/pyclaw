c
c --------------------------------------------------------------
c
      subroutine outtre(mlev,outgrd,nvar,naux)
c
      implicit double precision (a-h,o-z)
      logical  outgrd

      include  "call.i"
c
c ::::::::::::::::::::::: OUTTRE :::::::::::::::::::::::::::::::::::
c
c outtre - output subtree
c input parameters:
c    mlev   - ptr to subtree to output i.e., start at level(mlev)
c    outgrd - if true, output the values on the grid
c tree is output from 'level' to finest level.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      write (outunit,1)
1     format(1x,14hthe subtree is)
c
      level = node(nestlevel, mlev)
10    if (level .gt. lfine) go to 99
          mptr    = lstart(level)
 20       if (mptr .eq. 0) go to 30
              call outmsh(mptr,outgrd,nvar,naux)
              mptr = node(levelptr, mptr)
          go to 20
 30       continue
          write(outunit,2)level,iregst(level),jregst(level),
     .                    iregend(level),jregend(level)
 2        format(/,"grids at level ",i5," go from ",2i9," to",2i9,/)
          level = level + 1
      go to 10
c
 99   return
      end
