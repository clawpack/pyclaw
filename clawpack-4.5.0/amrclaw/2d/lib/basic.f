c
c ------------------------------------------------------
c
      subroutine basic (time,lst,lend)
c
      implicit double precision (a-h,o-z)

      include  "call.i"
c
c :::::::::::::::::::::: BASIC :::::::::::::::::::::::::
c  basic = outputs basic information needed by the other graphics
c          output routines (valout) at the given time,
c          write the entire levellist, from level 1 to lfine,
c          and the tree structure from level lst to lend.
c :::::::::::::::::::::::::::::;::::::::::::::::::::::::
c
      write(pltunit1,100) time
100   format(8h*TIME = ,f10.5)
      write(pltunit1,101) lfine, (lstart(i),i=1,lfine), nghost
101   format(10i6)
      write(pltunit1,105) xupper,yupper,xlower,ylower
105   format(4e15.8)
      write(pltunit1,102) lst, lend
102   format(2i6)
      write(pltunit1,106)(hxposs(i),i=1,lfine)
      write(pltunit1,106)(hyposs(i),i=1,lfine)
106   format(6e13.8)
c
      level = lst
 10   if (level .gt. lend) go to 99
          mptr = lstart(level)
 20       if (mptr .eq. 0) go to 30
              write(pltunit1,103) mptr, (node(i,mptr),i=1,nsize)
              write(pltunit1,104) (rnode(i,mptr),i=1,rsize)
103           format(10i7)
104           format(5e14.8)
              mptr = node(levelptr,mptr)
          go to 20
 30       level = level + 1
      go to 10
c
 99   return
      end
