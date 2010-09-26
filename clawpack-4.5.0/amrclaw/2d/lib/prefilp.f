c
c --------------------------------------------------------------
c
      recursive subroutine prefilrecur(level,nvar,
     1                       valbig,aux,naux,time,mitot,mjtot,
     1                       nrowst,ncolst,
     3                       ilo,ihi,jlo,jhi)

c
      implicit double precision (a-h,o-z)

      include  "call.i"

      dimension valbig(mitot,mjtot,nvar), aux(mitot,mjtot,naux)
      dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)

      dimension scratch(max(mitot,mjtot)*nghost*nvar)
      dimension scratchaux(max(mitot,mjtot)*nghost*naux)
      iadd(i,j,ivar)  = locflip + i - 1 + nr*((ivar-1)*nc+j-1)
      iaddscratch(i,j,ivar)  = i  + nr*((ivar-1)*nc+j-1)


c
c  :::::::::::::: PREFILRECUR :::::::::::::::::::::::::::::::::::::::::::
c     For periodic boundary conditions more work needed to fill the
c     piece of the boundary. This routine was
c     called because the patch sticks out of the domain,
c     and has periodic bc.s preprocess the patch before calling
c     filpatch to shift the patch periodically back into the domain.
c
c     Inputs to this routine:
c     xl, xr, yb, yt = the location in physical space of
c     corners of a patch.
c     ilo,ihi,jlo,jhi = the location in index space of this patch.
c
c     Outputs from this routine:
c     The values around the border of the grid are inserted
c     directly into the enlarged valbig array for this piece.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c
c     # will divide patch into 9 possibilities (some empty): 
c       x sticks out left, x interior, x sticks out right
c       same for y. for example, the max. would be
c       i from (ilo,-1), (0,iregsz(level)-1), (iregsz(level),ihi)
        
!--        write(dbugunit,*)" entering prefilrecur with level ",level
!--        write(dbugunit,*)"     and patch indices ilo,ihi,jlo,jhi ",
!--     &             ilo,ihi,jlo,jhi

        ist(1)  = ilo
        ist(2)  = 0
        ist(3)  = iregsz(level)
        iend(1) = -1
        iend(2) = iregsz(level)-1
        iend(3) = ihi

        jst(1)  = jlo
        jst(2)  = 0
        jst(3)  = jregsz(level)
        jend(1) = -1
        jend(2) = jregsz(level)-1
        jend(3) = jhi

        ishift(1) = iregsz(level)
        ishift(2) = 0
        ishift(3) = -iregsz(level)
        jshift(1) = jregsz(level)
        jshift(2) = 0
        jshift(3) = -jregsz(level)

       do 20 i = 1, 3
          i1 = max(ilo,  ist(i))
          i2 = min(ihi, iend(i))
       do 10 j = 1, 3
          j1 = max(jlo,  jst(j))
          j2 = min(jhi, jend(j))


          if ((i1 .le. i2) .and. (j1 .le. j2)) then ! part of patch in this region
c
c           there is something to fill. j=2 case is interior, no special
c           mapping needed even if spherical bc
            if (.not. spheredom .or. j .eq.2 ) then
              iputst = (i1 - ilo) + nrowst
              jputst = (j1 - jlo) + ncolst
              call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,
     1                      iputst,jputst,
     2                      i1+ishift(i),i2+ishift(i),
     3                      j1+jshift(j),j2+jshift(j))
            else
              nr = i2 - i1 + 1
              nc = j2 - j1 + 1
              ng = 0    ! no ghost cells in this little patch. fill everything.

              jbump = 0
              if (j1 < 0)   jbump = abs(j1)  ! bump up so new bottom is 0
              if (j2 >= jregsz(level)) jbump = -(j2+1-jregsz(level)) ! bump down

c             next 2 lines would take care of periodicity in x
              iwrap1 = i1 + ishift(i)
              iwrap2 = i2 + ishift(i)
c             next 2 lines take care of mapped sphere bcs
              iwrap1 = iregsz(level) - iwrap1 -1
              iwrap2 = iregsz(level) - iwrap2 -1
c             swap so that smaller one is left index, etc since mapping reflects
              tmp = iwrap1
              iwrap1 = iwrap2
              iwrap2 = tmp

              jwrap1 = j1 + jbump
              xlwrap = xlower + iwrap1*hxposs(level)
              ybwrap = ylower + jwrap1*hyposs(level)
              

c              locflip = igetsp(nr*nc*(nvar+naux))

c              write(dbugunit,11) nr,nc,mitot,mjtot
 11           format("scratch space nr*nc = ",i5,i5," grid is ",i5,i5)
c              locflipaux = locflip + nr*nc*nvar
!--              if (naux>0) call setaux(nr,nc,ng,nr,nc,xlwrap,ybwrap,
!--     1                                hxposs(level),hyposs(level),naux,
!--     2                                alloc(locflipaux))
              if (naux>0) call setaux(nr,nc,ng,nr,nc,xlwrap,ybwrap,
     1                                hxposs(level),hyposs(level),naux,
     2                                scratchaux)

c              write(dbugunit,101) i1,i2,j1,j2
c              write(dbugunit,102) iwrap1,iwrap2,j1+jbump,j2+jbump
 101          format(" actual patch from i:",2i5," j :",2i5)
 102          format(" filrecur called w i:",2i5," j :",2i5)
!--              call filrecur(level,nvar,alloc(locflip),alloc(locflipaux),
!--     1                      naux,time,nr,nc,1,1,iwrap1,iwrap2,
!--     2                      j1+jbump,j2+jbump)
              call filrecur(level,nvar,scratch,scratchaux,
     1                      naux,time,nr,nc,1,1,iwrap1,iwrap2,
     2                      j1+jbump,j2+jbump)

c             copy back using weird mapping for spherical folding
              do 15 ii = i1, i2
              do 15 jj = j1, j2
c              write(dbugunit,100) nrowst+ii-ilo,ncolst+jj-jlo,nr-(ii-i1),
c     1                            nc-jj+j1
 100          format(" filling loc ",2i5," with ",2i5)

              do 14 ivar = 1, nvar
                 valbig(nrowst+(ii-ilo),ncolst+(jj-jlo),ivar) = 
     1                 scratch(iaddscratch(nr-(ii-i1),nc-(jj-j1),ivar))
 14           continue
c              write(dbugunit,103)(valbig(nrowst+(ii-ilo),ncolst+(jj-jlo),
c     &                                  ivar),ivar=1,nvar)
 103          format(" new val is ",4e15.7)
 15           continue
             
c             call reclam(locflip, (nr*nc*(nvar+naux)))
            endif

          endif

 10    continue
 20    continue
      


      return
      end
