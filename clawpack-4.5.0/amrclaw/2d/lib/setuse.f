c
c ----------------------------------------------------------------
c
       subroutine setuse(listbc,maxsp,ispot,mkid,
     1                   ilo, ihi, jlo, jhi,
     2                   iclo,ichi,jclo,jchi,kflag)
c
c :::::::::::::::::::::::: SETUSE ::::::::::::::::::::::::::::::::
c
c set up boundary list for coarse grid, to be used by fluxsv. 
c loop around boundary of fine grids to do this.  each entry has
c     i, j, side #, fine grid #, loc in fine grid list for fluxes.
c  for example, side 1 of fine grid fixes side 3 of coarse grid,
c  so coarse grid list will store the # 3.
c  wrt coarse grid, the sides are:
c              2
c           1     3       that is, right edge of a coarse cell = 3
c              4                    top  edge of a coarse cell = 2
c
c  # lkid is the index into the fine grid's saved fluxes.
c  # the fine grid will save all its fluxes all around its
c  # perimeter. lkid tells where the coarse grid should
c  # taking them from. (no ghost cells in this index, but 
c  # it is 1-based for indexing array, not - based for
c  # integer index of grid location).
c
c  changed 11/11/08: spheredom for periodically mapped spherical
c      grids. could affect top and bottom if fine grid touches
c      edge of domain in y direction. if calling with spheredom
c      (and not yperdom) then grid is NOT periodically mapped.
c  need kflag to indicate spherically mapped now - otherwise
c  cant tell the difference, dont skip appropropriate loops
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      implicit double precision (a-h,o-z)
      dimension listbc(5,maxsp)

      include "call.i"

      ibc = ispot
      ist  = iclo - 1
      iend = ichi + 1
      jst  = jclo - 1
      jend = jchi + 1
c
c  left side (of fine grid, right side of coarse cell)
c
      if (ist .lt. ilo .or. kflag .ne. 1) go to 20
         lkid     = max(jlo,jclo) - jclo + 1
         do 10 j  = max(jlo,jclo), min(jhi,jchi)
            ispot              = ispot + 1
            listbc(1,ispot)    = ist-ilo+nghost+1
            listbc(2,ispot)    = j-jlo+nghost+1
            listbc(3,ispot)    = 3
            listbc(4,ispot)    = mkid
            listbc(5,ispot)    = lkid
            lkid               = lkid + 1
 10      continue
c
c   top side (of fine grid, bottom of coarse cell)
c
 20    if (kflag .eq. 1) then  ! regular interior case
         if (jend .gt. jhi) go to 40
         lkid       = (jchi-jclo+1) + max(ilo,iclo)-iclo + 1
         do 30 i    = max(ilo,iclo), min(ihi,ichi)
            ispot              = ispot + 1
            listbc(1,ispot)    = i-ilo+nghost+1
            listbc(2,ispot)    = jend-jlo+nghost+1
            listbc(3,ispot)    = 4
            listbc(4,ispot)    = mkid
            listbc(5,ispot)    = lkid
c	    write(outunit,595)ispot,(listbc(ipl,ispot),ipl=1,5)
 595        format("   entry ",i5," has ", 5i5)
            lkid               = lkid + 1
 30      continue
       else  if (kflag .eq. 2) then !spherical
c top side of a fine grid is also top side of a coarse cell due to mapping
c        write(outunit,*)":fixing top cells with fine grid ",mkid
c        original code was insanely complicated. look at all indices and decide.
         level  = node(nestlevel,mkid) - 1
         lkid   = (jchi-jclo+1)+ 1  ! starts here wrt fine grid. may not use on coarse grid
         do 31 i    = iclo, ichi  
            iwrap = iregsz(level) - i -1
            if (iwrap .ge. ilo .and. iwrap .le. ihi) then
               ispot              = ispot + 1
               listbc(1,ispot)    = iwrap - ilo + nghost + 1
               listbc(2,ispot)    = jend  - jlo + nghost   ! note adjustment of j (one less)
               listbc(3,ispot)    = 5 ! affects TOP of mapped coarse cell in diff. way
               listbc(4,ispot)    = mkid
               listbc(5,ispot)    = lkid
c	    write(outunit,595)ispot,(listbc(ipl,ispot),ipl=1,5)
          endif
            lkid               = lkid + 1  ! increment fine list loc even if not used
 31      continue

       endif
c
c  right side (of fine grid, left of coarse cell)
c (numbered from bottom to top, so not continuous in lkid numbering)
c
 40    if (iend .gt. ihi .or. kflag .ne. 1) go to 60
       lkid     = (ichi-iclo+1)+(jchi-jclo+1)
     .               + max(jlo,jclo) - jclo + 1
       do 50 j  = max(jlo,jclo), min(jhi,jchi)
          ispot              = ispot + 1
          listbc(1,ispot)    = iend-ilo+nghost+1
          listbc(2,ispot)    = j-jlo+nghost+1
          listbc(3,ispot)    = 1
          listbc(4,ispot)    = mkid
          listbc(5,ispot)    = lkid
          lkid   = lkid + 1
 50    continue
c
c  bottom side (of fine grid, top of coarse cell, unless spheredom)
c (numbered left to right, so not continuous in lkid numbering)
c
 60    if (kflag .eq. 1) then
          if (jst .lt. jlo) go to 80
          lkid =  2*(jchi-jclo+1)+(ichi-iclo+1) + max(ilo,iclo)-iclo + 1
          do 70 i      = max(ilo,iclo), min(ihi,ichi)
             ispot              = ispot + 1
             listbc(1,ispot)    = i-ilo+nghost+1
             listbc(2,ispot)    = jst-jlo+nghost+1
             listbc(3,ispot)    = 2
             listbc(4,ispot)    = mkid
             listbc(5,ispot)    = lkid
             lkid   = lkid + 1
 70      continue
       else  ! spherical
c bottom side of fine grid affects bottom of coarse cell 
c         fine grids saves fluxes in usual way
c         coarse grid only needs to change where to use them
          if (kflag .ne. 3) go to 80
c          write(outunit,*)":fixing bottom cells with fine grid ",mkid
          level  = node(nestlevel,mkid)-1
          lkid   = 2*(jchi-jclo+1) + (ichi-iclo+1) + 1
          do 71 i      = iclo, ichi
             iwrap = iregsz(level) - i - 1
             if (iwrap .ge. ilo .and. iwrap .le. ihi) then
                ispot              = ispot + 1
                listbc(1,ispot)    = iwrap - ilo + nghost + 1
                listbc(2,ispot)    = nghost+1   ! grid bottom is at zero index
                listbc(3,ispot)    = 6   ! affects BOTTOM of mapped coarse cell in diff. way
                listbc(4,ispot)    = mkid
                listbc(5,ispot)    = lkid 
c 	    write(outunit,595)ispot,(listbc(ipl,ispot),ipl=1,5)
             endif
             lkid   = lkid + 1
 71      continue       

       endif
c
 80    continue 
       return
       end
