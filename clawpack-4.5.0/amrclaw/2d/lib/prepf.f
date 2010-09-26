c
c ----------------------------------------------------------
c
      subroutine prepf(level,nvar,naux)
c
      implicit double precision (a-h,o-z)

      include  "call.i"
c
c ::::::::::::::::::::::::::: PREPF :::::::::::::::::::::::::::::
c
c prepare new fine grids to save fluxes after each integration step
c for future conservative fixing of coarse grids
c save all boundary fluxes of fine grid (even if phys. bndry.) -
c but only save space for every intrat. (remember - 4 fluxes).
c
c :::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::::::
c
      mkid = lstart(level)
 10   if (mkid .eq. 0) go to 99
          nx    = node(ndihi,mkid)-node(ndilo,mkid) + 1
          ny    = node(ndjhi,mkid)-node(ndjlo,mkid) + 1
          ikeep = nx/intratx(level-1)
          jkeep = ny/intraty(level-1)
          lenbc = 2*(ikeep+jkeep)
c
c         get twice the storage, one for plus or minus fluxes, the
c         other for coarse solution for wave fixing. also for aux vars.
c
          node(ffluxptr,mkid) = igetsp(2*nvar*lenbc+naux*lenbc)
          ist   = node(ffluxptr,mkid)

          do 20 i = 1, 2*nvar*lenbc + naux*lenbc
 20          alloc(ist+i-1) = 0.d0
          mkid = node(levelptr,mkid)
          go to 10

 99    return
       end
