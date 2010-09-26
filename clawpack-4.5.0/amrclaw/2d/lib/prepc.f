c
c ----------------------------------------------------------
c
      subroutine prepc(level,nvar)
c
      implicit double precision (a-h,o-z)

      include  "call.i"
c
c :::::::::::::::::::: PREPC ::::::::::::::::::::::::::::::::::::::
c
c this routine called because regridding just changed the fine grids.
c modify coarse grid boundary lists to store fluxes in appropriate
c fine grids lists.
c assume new fine grids have node(cfluxptr) initialized to null
c
c  first compute max. possible number of list cells. allocate
c  initially so that one pass through is enough.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      maxsp  = 0
      mkid   = lstart(level+1)
 10   if (mkid .eq. 0) go to 20
         ikeep  = (node(ndihi,mkid)-node(ndilo,mkid)+1)/intratx(level)
         jkeep  = (node(ndjhi,mkid)-node(ndjlo,mkid)+1)/intraty(level)
         maxsp  = maxsp + 2*(ikeep+jkeep)
      mkid = node(levelptr,mkid)
      go to 10
 20   listsp(level) = maxsp
      if (maxsp .eq. 0) go to 99
c
      hxpar   = hxposs(level)
      hypar   = hyposs(level)
      hxkid   = hxposs(level+1)
      hykid   = hyposs(level+1)
      imax    = iregsz(level) - 1
      jmax    = jregsz(level) - 1

      mpar = lstart(level)
 30   if (mpar .eq. 0) go to 99
c
       ispot   = 0
       ilo     = node(ndilo,mpar)
       jlo     = node(ndjlo,mpar)
       ihi     = node(ndihi,mpar)
       jhi     = node(ndjhi,mpar)
       locbc   = igetsp(5*maxsp)
c      #  initialize list to 0 (0 terminator indicates end of bc list)
       do 35 i = 1,5*maxsp
 35      alloc(locbc+i-1) = 0.d0
       node(cfluxptr,mpar) = locbc
c
       mkid = lstart(level+1)
 40    if (mkid .eq. 0) go to 60

          iclo = node(ndilo,mkid)/intratx(level)
          jclo = node(ndjlo,mkid)/intraty(level)
          ichi = node(ndihi,mkid)/intratx(level)
          jchi = node(ndjhi,mkid)/intraty(level)

          iplo = max(ilo,iclo)
          jplo = max(jlo,jclo)
          iphi = min(ihi,ichi)
          jphi = min(jhi,jchi)

c   regular intersections (will check in setuse that no duplication)
c   this first call is only interior interfaces. 

          if (iplo .le. iphi+1 .and. jplo .le. jphi+1) then
               kflag = 1 ! interior stuff, no mappings
                call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,iclo,ichi,jclo,jchi,kflag)
          endif

c   for fine grids touching periodic boundary on right
          if  (xperdom .and. ilo .eq. 0 .and. ichi .eq. imax) then
              kflag = 1 ! periodic in x
              call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,iclo-iregsz(level),ichi-iregsz(level),
     3          jclo,jchi,kflag)
           endif

c   for fine grids touching periodic boundary on left
          if  (xperdom .and. iclo .eq. 0 .and. ihi .eq. imax) then
              kflag = 1
              call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,iclo+iregsz(level),ichi+iregsz(level),
     3          jclo,jchi,kflag)
          endif

c   for fine grids touching periodic boundary on top
          if  (yperdom .and. jlo .eq. 0 .and. jchi .eq. jmax) then
                kflag = 1
                call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,iclo,ichi,
     3          jclo-jregsz(level),jchi-jregsz(level),kflag)
          endif

c   for fine grids touching periodic boundary on bottom
          if  (yperdom .and. jclo .eq. 0 .and. jhi .eq. jmax)  then
              kflag = 1
              call setuse(alloc(locbc),maxsp,ispot,mkid,
     2          ilo,ihi,jlo,jhi,iclo,ichi,
     3          jclo+jregsz(level),jchi+jregsz(level),kflag)
          endif

c   for fine grids touching boundary on top in spherically mapped case
c   and coarse grid touches top too. see if (mapped) x extent overlap.
          if  (spheredom .and. jhi .eq. jmax .and. jchi .eq. jmax) then
               kflag = 2
c              write(dbugunit,*)" for coarse grid ",mpar
               iwrap2 = iregsz(level) - iclo - 1  !higher mapped index
               iwrap1 = iregsz(level) - ichi - 1  !lower mapped index
               if (max(ilo,iwrap1) .le. min(ihi,iwrap2)) then
                  call setuse(alloc(locbc),maxsp,ispot,mkid,
     1                        ilo,ihi,jlo,jhi,iclo,ichi,
     2                        jclo,jchi,kflag)
               endif
          endif

c   fine grids touching boundary on bottom for spherically mapped case
c   coarse grid touches bottom too. see if (mapped) x extents overlap
          if  (spheredom .and. jclo .eq. 0 .and. jlo .eq. 0) then
               kflag = 3
               iwrap2 = iregsz(level) - iclo - 1  !higher mapped index
               iwrap1 = iregsz(level) - ichi - 1  !lower mapped index
               if (max(ilo,iwrap1) .le. min(ihi,iwrap2)) then
                  call setuse(alloc(locbc),maxsp,ispot,mkid,
     1                        ilo,ihi,jlo,jhi,iclo,ichi,
     2                        jclo,jchi,kflag)
               endif
          endif

 50     mkid = node(levelptr,mkid)
        go to 40
c
c  done with subgrid cycle. if no cells would need fixing, all done
c  else cycle through again to set up list with info. for bc processing
c
 60     continue
c
c  for now, leave unused space allocated to the grid. alternative is to
c  return (maxsp-ispot) amt starting at loc node(cfluxptr,mpar)+ispot.
c
       mpar = node(levelptr,mpar)
       go to 30
c
 99    return
       end
