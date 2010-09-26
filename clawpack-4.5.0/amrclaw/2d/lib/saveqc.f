c
c  ================================================================
      subroutine saveqc(level,nvar,naux)
c  ================================================================
c
      implicit double precision (a-h,o-z)

      include  "call.i"
      logical sticksout
c
c ::::::::::::::::::::::::: SAVEQC :::::::::::::::::::::::::::::::::
c prepare new fine grids to save fluxes after each integration step
c for future conservative fix-up on coarse grids.
c save all boundary fluxes of fine grid (even if on a  phys. bndry.) -
c but only save space for every intrat of them. 
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      levc = level - 1
      hxc  = hxposs(levc)
      hyc  = hyposs(levc)

      mkid = lstart(level)
 10   if (mkid .eq. 0) go to 99
          nx    = node(ndihi,mkid)-node(ndilo,mkid) + 1
          ny    = node(ndjhi,mkid)-node(ndjlo,mkid) + 1
          ikeep = nx/intratx(level-1)
          jkeep = ny/intraty(level-1)
          lenbc = 2*(ikeep+jkeep)
          ist   = node(ffluxptr,mkid)
          time = rnode(timemult,mkid)

c         make coarsened enlarged patch for conservative fixup
          ilo = node(ndilo,mkid)
          jlo = node(ndjlo,mkid)
          ihi = node(ndihi,mkid)
          jhi = node(ndjhi,mkid)
          iclo = ilo/intratx(level-1) - 1
          jclo = jlo/intraty(level-1) - 1
          ichi = (ihi+1)/intratx(level-1)
          jchi = (jhi+1)/intraty(level-1)
          nrow = ichi-iclo+1
          ncol = jchi-jclo+1
          xl   = rnode(cornxlo,mkid) - hxc
          yb   = rnode(cornylo,mkid) - hyc
          xr   = rnode(cornxhi,mkid) + hxc
          yt   = rnode(cornyhi,mkid) + hyc
          loctmp = igetsp(nrow*ncol*(nvar+naux))
          loctx  = loctmp + nrow*ncol*nvar
          locaux = node(storeaux,mkid)

          if (iclo .lt. 0 .or. ichi .eq. iregsz(levc) .or.
     &        jclo .lt. 0 .or. jchi .eq. jregsz(levc)) then
            sticksout = .true.
          else
            sticksout = .false.
          endif

          if (sticksout .and. (xperdom.or.yperdom.or.spheredom)) then
             iperim = nrow*ncol
             locflip = igetsp(iperim*(nvar+naux))
             call preicall(alloc(loctmp),alloc(loctx),nrow,ncol,nvar,
     .                     naux,iclo,ichi,jclo,jchi,level-1,locflip)
             call reclam(locflip,iperim*(nvar+naux))
          else 
             call icall(alloc(loctmp),alloc(loctx),nrow,ncol,nvar,naux,
     .                   iclo,ichi,jclo,jchi,level-1,1,1,sticksout)
          endif
          call bc2amr(alloc(loctmp),alloc(loctx),nrow,ncol,nvar,naux,
     .                hxc,hyc,level,time,
     .                xl,xr,yb,yt,
     .                xlower,ylower,xupper,yupper,
     .                xperdom,yperdom,spheredom)
          call cstore(alloc(loctmp),nrow,ncol,nvar,
     .                alloc(ist+nvar*lenbc),lenbc,naux,alloc(loctx),
     .                alloc(ist+2*nvar*lenbc))
          call reclam(loctmp,nrow*ncol*(nvar+naux))

          mkid = node(levelptr,mkid)
          go to 10
 99    return
       end
