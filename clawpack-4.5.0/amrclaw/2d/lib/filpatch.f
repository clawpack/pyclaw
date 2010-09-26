c
c ---------------------------------------------------------------
c
        recursive subroutine filrecur(level,nvar,valbig,aux,naux,
     1                      time,mitot,mjtot,
     2                      nrowst,ncolst,ilo,ihi,jlo,jhi)

c :::::::::::::::::::::::::::: FILPATCH :::::::::::::::::::::::::;
c
c  fill the portion of valbig from rows  nrowst
c                             and  cols  ncolst
c  the patch can also be described by the corners (xlp,ybp) by (xrp,ytp).
c  vals are needed at time time , and level level,
c
c  first fill with  values obtainable from the level level
c  grids. if any left unfilled, then enlarge remaining rectangle of
c  unfilled values by 1 (for later linear interp), and recusively
c  obtain the remaining values from  coarser levels.
c
c :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;

      implicit double precision (a-h,o-z)

      include  "call.i"

      logical   set, sticksout
      dimension valbig(mitot,mjtot,nvar), aux(mitot,mjtot,naux)

c  use stack-based scratch arrays instead of alloc, since dont really
c  need to save beyond these routines, and to allow dynamic memory resizing
c
c     use 1d scratch arrays that are potentially the same size as 
c     current grid, since may not coarsen.
c     need to make it 1d instead of 2 and do own indexing, since
c     when pass it in to subroutines they treat it as having different
c     dimensions than the max size need to allocate here
c
!--      dimension valcrse((ihi-ilo+2)*(jhi-jlo+2)*nvar)  ! NB this is a 1D array 
!--      dimension auxcrse((ihi-ilo+2)*(jhi-jlo+2)*naux)  ! the +2 is to expand on coarse grid to enclose fine
c ### turns out you need 3 rows, forget offset of 1 plus one on each side
      dimension valcrse((ihi-ilo+3)*(jhi-jlo+3)*nvar)  ! NB this is a 1D array 
      dimension auxcrse((ihi-ilo+3)*(jhi-jlo+3)*naux)  ! the +3 is to expand on coarse grid to enclose fine
c
      dimension flaguse(ihi-ilo+1,jhi-jlo+1)

c      iadflag(i,j) =  locuse + i-1+(j-1)*nrowp  ! no longer used
c      ivalc(i,j,ivar) = loccrse + (i - 1) + nrowc*(j - 1)
      ivalc(i,j,ivar) = i + nrowc*(j - 1)
     &                    + nrowc*ncolc*(ivar-1)
      sticksout(iplo,iphi,jplo,jphi)  =
     &            (iplo .lt. 0 .or. jplo .lt. 0 .or.
     &             iphi .ge. iregsz(levc) .or. jphi .ge. jregsz(levc))


!--      write(*,*)" entering filrecur with level ",level
!--      write(*,*)"     and patch indices ilo,ihi,jlo,jhi ",
!--     &             ilo,ihi,jlo,jhi
c
c We begin by filling values for grids at level level. If all values can be
c filled in this way, we return;

        nrowp   = ihi - ilo + 1
        ncolp   = jhi - jlo + 1
c        locuse  = igetsp(nrowp*ncolp)
        hxf     = hxposs(level)
        hyf     = hyposs(level)
        xlp     = xlower + ilo*hxf
        xrp     = xlower + (ihi+1)*hxf
        ybp     = ylower + jlo*hyf
        ytp     = ylower + (jhi+1)*hyf

        call intfil
     &  (valbig,mitot,mjtot,time,flaguse,nrowst,ncolst,
     &   ilo,ihi,jlo,jhi,level,nvar,naux)
c     &  (valbig,mitot,mjtot,time,locuse,nrowst,ncolst,

c
c Trimbd returns set = true if all of the entries are filled (=1.).
c set = false, otherwise. If set = true, then no other levels are
c are required to interpolate, and we return.
c
c Note that the used array is filled entirely in intfil, i.e. the
c marking done there also takes  into account the points filled by
c the boundary conditions. bc2amr will be called later, after all 4
c boundary pieces filled.

c        call trimbd(alloc(locuse),nrowp,ncolp,set,il,ir,jb,jt)
        call trimbd(flaguse,nrowp,ncolp,set,il,ir,jb,jt)

        if (set) go to 90 ! all done except for bcs
c
c otherwise make recursive calls to coarser levels to fill remaining unset points
c
        if (level .eq. 1) then
           write(outunit,*)" error in filrecur - level 1 not set"
           write(outunit,900) nrowst,ncolst
           write(*,*)" error in filrecur - level 1 not set"
           write(*,*)" should not need more recursion "
           write(*,*)" to set patch boundaries"
           write(*,900) nrowst,ncolst
900        format("start at row: ",i4," col ",i4)
           stop
        endif

c set = false. we will have to interpolate some values from coarser
c levels. We begin by initializing the level level arrays, so that we can use
c purely recursive formulation for interpolating.


        levc = level - 1
        hxc  = hxposs(levc)
        hyc  = hyposs(levc)

        isl  = il + ilo - 1
        isr  = ir + ilo - 1
        jsb  = jb + jlo - 1
        jst  = jt + jlo - 1
c
c       coarsen
        lratiox = intratx(levc)
        lratioy = intraty(levc)
        iplo   = (isl-lratiox+nghost*lratiox)/lratiox - nghost
        jplo   = (jsb-lratioy+nghost*lratioy)/lratioy - nghost
        iphi   = (isr+lratiox)/lratiox
        jphi   = (jst+lratioy)/lratioy

        xlc  =  xlower + iplo*hxc
        ybc  =  ylower + jplo*hyc
        xrc  =  xlower + (iphi+1)*hxc
        ytc  =  ylower + (jphi+1)*hyc

        nrowc   =  iphi - iplo + 1
        ncolc   =  jphi - jplo + 1
        ntot    = nrowc*ncolc*(nvar+naux)
c        write(*,876) nrowc,ncolc, ihi-ilo+2,jhi-jlo+2
c        write(*,876) nrowc,ncolc, ihi-ilo+3,jhi-jlo+3
 876    format(" needed coarse grid size ",2i5," allocated ",2i5)
        if (nrowc .gt. ihi-ilo+3 .or. ncolc .gt. jhi-jlo+3) then
            write(*,*)" did not make big enough work space in filrecur "
            write(*,*)" need coarse space with nrowc,ncolc ",nrowc,ncolc
            write(6,*)" made space for ilo,ihi,jlo,jhi ",ilo,ihi,jlo,jhi
            stop
        endif
c        loccrse = igetsp(ntot)
c        locauxc = loccrse + nrowc*ncolc*nvar
        if (naux.gt.0) then
              maxmx = nrowc - 2*nghost
              mx = maxmx
              maxmy = ncolc - 2*nghost
              my = maxmy
              xl = xlc + nghost*hxc
              yb = ybc + nghost*hyc
              call setaux(maxmx,maxmy,nghost,mx,my,xl,yb,hxc,hyc,
     &                    naux,auxcrse)
c     &                    naux,alloc(locauxc))
        endif

        if ((xperdom .or. (yperdom .or. spheredom)) .and.
     &       sticksout(iplo,iphi,jplo,jphi)) then
             call prefilrecur(levc,nvar,valcrse,auxcrse,
     1                        naux,time,nrowc,ncolc,1,1,
     2                        iplo,iphi,jplo,jphi)
        else
c           call filrecur(levc,nvar,alloc(loccrse),alloc(locauxc),naux,
           call filrecur(levc,nvar,valcrse,auxcrse,naux,
     1                   time,nrowc,ncolc,1,1,
     2                   iplo,iphi,jplo,jphi)
        endif

c       interpolate back up

20      continue

        do 100 iff = 1,nrowp
          ic = 2 + (iff-(isl-ilo)-1)/lratiox
          eta1 = (-0.5d0+dble(mod(iff-1,lratiox)))/dble(lratiox)

        do 100 jf  = 1,ncolp
          jc = 2 + (jf -(jsb-jlo)-1)/lratioy
          eta2 = (-0.5d0+dble(mod(jf -1,lratioy)))/dble(lratioy)

c           flag = alloc(iadflag(iff,jf))
           flag = flaguse(iff,jf)
           if (flag .eq. 0.0) then

c               xif = xlp + (.5 + float(iff-1))*hxf
c               yjf = ybp + (.5 + float(jf -1))*hyf

c               ic=idint((xif-xlc+.5*hxc)/hxc)
c               jc=idint((yjf-ybc+.5*hyc)/hyc)

c               xc = xlc + (float(ic) - .5)*hxc
c               yc = ybc + (float(jc) - .5)*hyc

c               eta1 = (xif - xc)/hxc
c               eta2 = (yjf - yc)/hyc

                do 101 ivar = 1,nvar
!--
!--                   valp10 = alloc(ivalc(ic+1,jc,ivar))
!--                   valm10 = alloc(ivalc(ic-1,jc,ivar))
!--                   valc   = alloc(ivalc(ic  ,jc,ivar))
!--                   valp01 = alloc(ivalc(ic  ,jc+1,ivar))
!--                   valm01 = alloc(ivalc(ic  ,jc-1,ivar))

                   valp10 = valcrse(ivalc(ic+1,jc,ivar))
                   valm10 = valcrse(ivalc(ic-1,jc,ivar))
                   valc   = valcrse(ivalc(ic  ,jc,ivar))
                   valp01 = valcrse(ivalc(ic  ,jc+1,ivar))
                   valm01 = valcrse(ivalc(ic  ,jc-1,ivar))

                   dupc = valp10 - valc
                   dumc = valc   - valm10
                   ducc = valp10 - valm10
                   du   = dmin1(dabs(dupc),dabs(dumc))
                   du   = dmin1(2.d0*du,.5d0*dabs(ducc))
                   fu = dmax1(0.d0,dsign(1.d0,dupc*dumc))

                   dvpc = valp01 - valc
                   dvmc = valc   - valm01
                   dvcc = valp01 - valm01
                   dv   = dmin1(dabs(dvpc),dabs(dvmc))
                   dv   = dmin1(2.d0*dv,.5d0*dabs(dvcc))
                   fv = dmax1(0.d0,dsign(1.d0,dvpc*dvmc))

                   valint = valc + eta1*du*dsign(1.d0,ducc)*fu
     .                           + eta2*dv*dsign(1.d0,dvcc)*fv

c                  valc00 = alloc(ivalc(ic,jc,ivar))
c                  valc10 = alloc(ivalc(ic+1,jc,ivar))
c                  valc01 = alloc(ivalc(ic,jc+1,ivar))
c                  valc11 = alloc(ivalc(ic+1,jc+1,ivar))
c                  valint = (1. - eta2)*
c    &               ((1. - eta1)*valc00 + eta1*valc10)
c    &               + eta2*((1. - eta1)*valc01 + eta1*valc11)

                   valbig(iff+nrowst-1,jf+ncolst-1,ivar) = valint

101             continue

           endif

100     continue

c        call reclam(loccrse,ntot)

 90     continue
c
c  set bcs, whether or not recursive calls needed. set any part of patch that stuck out
c
        call bc2amr(valbig,aux,mitot,mjtot,nvar,naux,
     1              hxf,hyf,level,time,
     2              xlp,xrp,ybp,ytp,
     3              xlower,ylower,xupper,yupper,
     4              xperdom,yperdom,spheredom)


c        call reclam(locuse,nrowp*ncolp)

        return
        end
