c
c -----------------------------------------------------------
c
      subroutine flglvl(nvar,naux,lcheck,nxypts,index,lbase,ldom2,
     .                  npts,t0)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

c
c :::::::::::::::::::: FLGLVL :::::::::::::::::::::::::::::::::
c
c flglvl = controls the error estimation/flagging bad pts. for
c          an entire level of grids.  returns pointer into alloc
c          where the (x,y) coordinations of the flagged pts. are.
c input parameters:
c           lcheck = level to be flagged
c output parameters:
c           nxypts = no. of flagged pts. total
c           index  = starting index in alloc of the flagged pts.
c                    (which occupy 2*nxypts locations).
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
      nxypts = 0
c
c   reserve space for entire domain worth of flagged points at
c   level lcheck. bits would be better, but integer will do
c   dom2 - holds domain flags
c   dom  - holds flagged pts.
c   dom3 - scratch
c
      isize = iregsz(lcheck)
      jsize = jregsz(lcheck)
      ibytesPerDP = 8
      ldom  = igetsp((isize+2)*(jsize+2)/ibytesPerDP+1)
c    
c   prepare domain in ldom2 (so can use ldom as scratch array before 
c   putting in the flags)
c
      idim = iregsz(lbase)
      jdim = jregsz(lbase)
      call domprep(alloc(ldom2),lbase,idim,jdim)

      call domshrink(alloc(ldom2),alloc(ldom),idim,jdim)

      do 6 lev = lbase+1, lcheck
         call domup(alloc(ldom2),alloc(ldom),idim,jdim,
     1              intratx(lev-1)*idim,intraty(lev-1)*jdim,lev-1)
         idim = intratx(lev-1)*idim
         jdim = intraty(lev-1)*jdim
         call domshrink(alloc(ldom2),alloc(ldom),idim,jdim)
 6    continue
c     # finish by transferring from iflags to iflags2
      call domcopy(alloc(ldom2),alloc(ldom),isize,jsize)
c
      numbad = 0
c     always call spest to set up stuff (initialize iflags, fill locbig)
c      call spest(nvar,naux,lcheck,alloc(ldom),isize,jsize,t0)
c     ### modified to pass in ldom instead of alloc(ldom) - called iflags in spest -
c     ###  since spest calls igetsp, if alloc is resized and moved, need relative
c     ### indexing, or iflags would have invalid address on the inside
      call spest(nvar,naux,lcheck,ldom,isize,jsize,t0)
      if (tol .gt. 0.) call errest(nvar,naux,lcheck)

      if (ibuff .gt. 0) then ! get scratch storage for bufnst
         ibytesPerDP = 8
         ldom3 = igetsp((isize+2)*(jsize+2)/ibytesPerDP+1)  ! incase need to resize
      endif
      call bufnst(nvar,naux,numbad,lcheck,alloc(ldom),isize,jsize,ldom3)
      if (ibuff .gt. 0) then ! return scratch storage for bufnst
         call reclam(ldom3,(isize+2)*(jsize+2)/ibytesPerDP+1)
      endif

      nxypts = nxypts + numbad
c
c  colate flagged pts into flagged points array
c
      if (nxypts .gt. 0) then
          index = igetsp(2*nxypts)
          call colate(alloc(index),nxypts,lcheck,
     1                alloc(ldom),alloc(ldom2),isize,jsize,npts)
      else 
         npts = nxypts
      endif

      call reclam(ldom,  (isize+2)*(jsize+2)/ibytesPerDP+1) 

      return
      end
