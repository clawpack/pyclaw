c
c -----------------------------------------------------------
c
      subroutine update (level, nvar)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      iadd(i,j,ivar)  = loc     + i - 1 + mitot*((ivar-1)*mjtot+j-1)
      iaddf(i,j,ivar) = locf    + i - 1 + mi*((ivar-1)*mj  +j-1)
      iaddfaux(i,j)   = locfaux + i - 1 + mi*((mcapa-1)*mj + (j-1))
      iaddcaux(i,j)   = loccaux + i - 1 + mitot*((mcapa-1)*mjtot+(j-1))
c
c :::::::::::::::::::::::::: UPDATE :::::::::::::::::::::::::::::::::
c update - update all grids at level 'level'.
c          this routine assumes cell centered variables.
c          the update is done from 1 level finer meshes under it.
c input parameter:
c    level  - ptr to the only level to be updated. levels coarser than
c             this will be at a diffeent time.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      lget = level
      if (uprint) write(outunit,100) lget
100   format(19h    updating level ,i5)
c
c  grid loop for each level
c
      dt     = possk(lget)

      mptr = lstart(lget)
 20   if (mptr .eq. 0) go to 85
         loc     = node(store1,mptr)
         loccaux = node(storeaux,mptr)
         nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         mitot   = nx + 2*nghost
         mjtot   = ny + 2*nghost
         ilo     = node(ndilo,mptr)
         jlo     = node(ndjlo,mptr)
         ihi     = node(ndihi,mptr)
         jhi     = node(ndjhi,mptr)
c
         if (node(cfluxptr,mptr) .eq. 0) go to 25
         locuse = igetsp(mitot*mjtot)
         call upbnd(alloc(node(cfluxptr,mptr)),alloc(loc),nvar,
     1              mitot,mjtot,listsp(lget),alloc(locuse),mptr)
         call reclam(locuse,mitot*mjtot)
c
c  loop through all intersecting fine grids as source updaters.
c
 25      mkid = lstart(lget+1)
 30        if (mkid .eq. 0) go to 80
           iclo   = node(ndilo,mkid)/intratx(lget)
           jclo   = node(ndjlo,mkid)/intraty(lget)
           ichi   = node(ndihi,mkid)/intratx(lget)
           jchi   = node(ndjhi,mkid)/intraty(lget)

           mi      = node(ndihi,mkid)-node(ndilo,mkid) + 1 + 2*nghost
           mj      = node(ndjhi,mkid)-node(ndjlo,mkid) + 1 + 2*nghost
           locf    = node(store1,mkid)
           locfaux = node(storeaux,mkid)
c
c  calculate starting and ending indices for coarse grid update, if overlap
c
         iplo = max(ilo,iclo)
         jplo = max(jlo,jclo)
         iphi = min(ihi,ichi)
         jphi = min(jhi,jchi)

         if (iplo .gt. iphi .or. jplo .gt. jphi) go to 75
c
c  calculate starting index for fine grid source pts.
c
         iff    = iplo*intratx(lget) - node(ndilo,mkid) + nghost + 1
         jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
         totrat = intratx(lget) * intraty(lget)
 
         do 71 i = iplo-ilo+nghost+1, iphi-ilo+nghost+1
         do 70 j = jplo-jlo+nghost+1, jphi-jlo+nghost+1
           if (uprint) then
              write(outunit,101) i,j,mptr,iff,jff,mkid
 101          format(' updating pt. ',2i4,' of grid ',i3,' using ',2i4,
     1               ' of grid ',i4)
              write(outunit,102)(alloc(iadd(i,j,ivar)),ivar=1,nvar)
 102          format(' old vals: ',4e12.4)
           endif
c
c
c  update using intrat fine points in each direction
c
           do 35 ivar = 1, nvar
 35           alloc(iadd(i,j,ivar)) = 0.d0
c
           if (mcapa .eq. 0) then
               do 50 jco  = 1, intraty(lget)
               do 50 ico  = 1, intratx(lget)
               do 40 ivar = 1, nvar
                 alloc(iadd(i,j,ivar))= alloc(iadd(i,j,ivar)) + 
     1                        alloc(iaddf(iff+ico-1,jff+jco-1,ivar))
 40              continue
 50            continue
            do 60 ivar = 1, nvar
 60          alloc(iadd(i,j,ivar)) = alloc(iadd(i,j,ivar))/totrat
               
           else

               do 51 jco  = 1, intraty(lget)
               do 51 ico  = 1, intratx(lget)
               capa = alloc(iaddfaux(iff+ico-1,jff+jco-1))
               do 41 ivar = 1, nvar
                 alloc(iadd(i,j,ivar))= alloc(iadd(i,j,ivar)) + 
     1                       alloc(iaddf(iff+ico-1,jff+jco-1,ivar))*capa
 41              continue
 51            continue
            do 61 ivar = 1, nvar
 61          alloc(iadd(i,j,ivar)) = alloc(iadd(i,j,ivar))/
     1                               (totrat*alloc(iaddcaux(i,j)))
           endif
c
            if (uprint) write(outunit,103)(alloc(iadd(i,j,ivar)),
     .                                     ivar=1,nvar)
 103        format(' new vals: ',4e12.4)
c
           jff = jff + intraty(lget)
 70        continue
           iff = iff + intratx(lget)
           jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
 71        continue
c
 75         mkid = node(levelptr,mkid)
            go to 30
c
 80         mptr = node(levelptr, mptr)
            go to 20
c
 85       continue
c
 99   return
      end
