c
c ----------------------------------------------------------
c
      subroutine putsp(lbase,level,nvar,naux)
c
      implicit double precision (a-h,o-z)

      include  "call.i"
c
c ::::::::::::::::::::::::::::::: PUTSP :::::::::::::::::::::::::
c
c reclaim list space in nodes cfluxptr and ffluxptr for grids at level
c
c first compute max. space allocated in node cfluxptr.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      if (level .eq. lfine) go to 30
c
      mptr  = lstart(level)
 20      call reclam(node(cfluxptr,mptr), 5*listsp(level))
         node(cfluxptr,mptr) = 0
      mptr  = node(levelptr,mptr)
      if (mptr .ne. 0) go to 20
c
 30    if (level .eq. lbase) go to 99
      mptr = lstart(level)
 40       nx    = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny    = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          ikeep = nx/intratx(level-1)
          jkeep = ny/intraty(level-1)
          lenbc = 2*(ikeep+jkeep)
c         twice perimeter since saving plus or minus fluxes 
c         plus coarse solution storage
          call reclam(node(ffluxptr,mptr), 2*nvar*lenbc+naux*lenbc)
          mptr  = node(levelptr,mptr)
          if (mptr .ne. 0) go to 40
c
 99    return
       end
