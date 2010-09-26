c
c ---------------------------------------------------------
c
      subroutine restrt(nsteps,time,nvar)
c
      implicit double precision (a-h,o-z)

      logical   ee

      include  "call.i"

      dimension intrtx(maxlv),intrty(maxlv),intrtt(maxlv)
      character*16  chkname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id
c
c     # HDF: Arrays which will hold output arrays.
c
      dimension  iqout(14), qout(4)
c
c     # HDF: Declare external HDF functions
c
      integer  sfstart, sfend
      external sfstart, sfend
c
c     # HDF: Set up HDF constants
c
      integer 	DFACC_READ, DFACC_WRITE, DFACC_CREATE 
      parameter(DFACC_READ = 1, DFACC_WRITE = 2, DFACC_CREATE = 4)

      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
c :::::::::::::::::::::::::::: RESTRT ::::::::::::::::::::::::::::::::
c read back in the check point files written by subr. check.
c
c some input variables might have changed, and also the
c alloc array could have been written with a smaller size at checkpoint
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     ###  make the file name showing the time step
c
      chkname = 'restart.data.hdf'
c     
c     # HDF: open hdf restart file.
c     
      sd_id = sfstart(chkname,DFACC_READ)
      if (sd_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to open HDF file',
     &        ' (call to sfstart in restrt_hdf.f)'
         STOP
      end if

      index = 0
      call read_integer_vector(sd_id,14,index,'iqout ',iqout)
      lenmax  = iqout(1)
      lendim  = iqout(2)
      isize   = iqout(3)
      lenf    = iqout(4)
      ibuf    = iqout(5)
      mstart  = iqout(6)
      ndfree  = iqout(7)
      lfine   = iqout(8)
      iorder  = iqout(9)
      mxnold  = iqout(10)
      kcheck1 = iqout(11)
      nsteps  = iqout(12)
      matlabu = iqout(13)
      lentot  = iqout(14)

      index = index + 1
      call read_integer_vector(sd_id,maxlv,index,'icheck',icheck)
      index = index + 1
      call read_integer_vector(sd_id,maxlv,index,'lstart',lstart)
      index = index + 1
      call read_integer_vector(sd_id,maxlv,index,'newstl',newstl)
      index = index + 1
      call read_integer_vector(sd_id,maxlv,index,'listsp',listsp)
      index = index + 1
      call read_integer_vector(sd_id,maxlv,index,'intratx',intrtx)
      index = index + 1
      call read_integer_vector(sd_id,maxlv,index,'intraty',intrty)
      index = index + 1
      call read_integer_vector(sd_id,maxlv,index,'kratio',intrtt)
      index = index + 1
      call read_integer_vector(sd_id,maxlv,index,'iregsz',iregsz)
      index = index + 1
      call read_integer_vector(sd_id,maxlv,index,'jregsz',jregsz)
      index = index + 1
      call read_integer_array(sd_id,lfdim,2,index,'lfree ',lfree)
      index = index + 1
      call read_integer_array(sd_id,rsize,maxgr,index,'node  ',node)

c     # HDF: DFNT_FLOAT64

      index = index + 1
      call read_double_vector(sd_id,4,index,'qout  ',qout)
      tl  = qout(1)
      time = qout(2)
      evol = qout(3)
      rvol = qout(4)

      index = index + 1
      call read_double_vector(sd_id,maxlv,index,'hxposs',hxposs)
      index = index + 1
      call read_double_vector(sd_id,maxlv,index,'hyposs',hyposs)
      index = index + 1
      call read_double_vector(sd_id,maxlv,index,'possk ',possk)
      index = index + 1
      call read_double_vector(sd_id,10,index,'rvoll ',rvoll)
      index = index + 1
      call read_double_vector(sd_id,lendim,index,'alloc ',alloc)
      
      index = index + 1
      call read_double_array(sd_id,rsize,maxgr,index,'rnode ',rnode)
c
c     # HDF: Close HDF file.
c
      istat = sfend(sd_id)
      if (istat.eq.FAIL) then
         WRITE(*,*) 'Failed to close SDS',
     &        ' (call to sfend in restrt_hdf.f)'
         STOP
      end if

      write(outunit,100) nsteps,time
      write(6,100) nsteps,time
  100 format(/,' RESTARTING the calculation after ',i5,' steps',
     1      /,'  (time = ',e15.7,')')

      do i = 1, mxnold-1
        if ( (intratx(i) .ne. intrtx(i)) .or.
     .       (intraty(i) .ne. intrty(i)) .or.
     .       (kratio(i) .ne.  intrtt(i)) ) then
        write(outunit,*) 
     .  " not allowed to change existing refinement ratios on Restart"
        write(*,*)
     .  " not allowed to change existing refinement ratios on Restart"
        write(outunit,*)" Old ratios:"
        write(*,*)      " Old ratios:"
        write(outunit,903)(intrtx(j),j=1,mxnold-1)
        write(*,903)      (intrtx(j),j=1,mxnold-1)
        write(outunit,903)(intrty(j),j=1,mxnold-1)
        write(*,903)      (intrty(j),j=1,mxnold-1)
        write(outunit,903)(intrtt(j),j=1,mxnold-1)
        write(*,903)      (intrtt(j),j=1,mxnold-1)
 903    format(6i3)
        stop
       endif
      end do
c
c     adjust free list of storage in case size has changed.
c
      idif = memsize - isize
      if (idif .gt. 0) then
         lfree(lenf,1) = isize + 2
         call reclam(isize+1,idif)
      else if (idif .lt. 0) then
         write(outunit,900) isize, memsize
         write(*,900)       isize, memsize
  900    format(' size of alloc not allowed to shrink with restart ',/,
     .         ' old size ',i7,' current size  ',i7)
         stop
      endif
c
c     adjust storage in case mxnest has changed - only allow it to increase,
c     and only at non-multiples of error estimation on old mxnest.
c
      if (mxnest .eq. mxnold) go to 99

      if (mxnest .lt. mxnold) then
         if (lfine .lt. mxnest) then
            go to 99
         else
            write(outunit,901) mxnold, mxnest
            write(*,      901) mxnold, mxnest
  901       format(' only allow mxnest to increase: ',/,
     &            '  old mxnest ',i4, ' new mxnest ',i4)
            stop
         endif
      endif

c     see if simple enough situation to allow changing mxnest
      ee = .false.
      do 10 level = 1, mxnold
         if (icheck(level) .ge. kcheck) then
            ee = .true.
            kmust = icheck(level)
         endif
   10 continue
      if (ee) then
         write(outunit,902) mxnold, mxnest, kmust
         write(*      ,902) mxnold, mxnest, kmust
  902    format(/,' only allow changes in mxnest (from ',
     &         i4,' to ',i4,')',/,
     &         ' when not time to error estimate: ',/,
     &         ' please run a few more steps before changing ',/,
     &         ' so that # of steps not greater then kcheck',/,
     &         ' or make kcheck > ',i4 )
         stop
      else
c        #  add second storage location to previous mxnest level
         mptr = lstart(mxnold)
   15    if (mptr .eq. 0) go to 25
         mitot = node(ndihi,mptr)-node(ndilo,mptr)+1+2*nghost
         mjtot = node(ndjhi,mptr)-node(ndjlo,mptr)+1+2*nghost
         node(store2,mptr) = igetsp(mitot*mjtot*nvar)
         mptr = node(levelptr,mptr)
         go to 15
   25    continue
      endif
c
c     # add new info. to spatial and counting arrays
   99 level = lfine + 1
   35 if (level .gt. mxnest) go to 45
      hxposs(level) = hxposs(level-1) / dble(intratx(level-1))
      hyposs(level) = hyposs(level-1) / dble(intraty(level-1))
      possk (level) = possk (level-1) / dble(kratio(level-1))
      iregsz(level) = iregsz(level-1) * dble(intratx(level-1))
      jregsz(level) = jregsz(level-1) * dble(intraty(level-1))
      level         = level + 1
      go to 35
   45 continue
c
c
      return
      end
c=================================================================
c     # HDF: NOTE THAT ALL INTEGER HDF CHECKPOINTING ROUTINES FOR
c     #      BOTH WRITING AND READING ARE INCLUDED HERE, ALTHOUGH
c     #      THE INTEGER WRITING ROUTINES ARE CALLED IN check_hdf.f.
c     #      THIS IS DONE TO PREVENT COMPILER ERRORS SINCE THE SAME
c     #      HDF ROUTINES ARE USED TO READ BOTH INTEGER AND DOUBLE
c     #      PRECISION ARRAYS.
c=================================================================
      subroutine read_integer_vector(sd_id,idims,index,qname,iout)
      implicit double precision (a-h,o-z)
      character*6  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id
      dimension  iout(idims), istart(1), istride(1), iedges(1)
c
c     # HDF: Declare external HDF functions
c
      integer  sfcreate, sfrdata, sfselect, sfendacc
      external sfcreate, sfrdata, sfselect, sfendacc
c
c     # HDF: Set up HDF constants
c
      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
c     # HDF: Select dataset in HDF file.
c
      sds_id = sfselect(sd_id,index)
      if (sds_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to select data set for  variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfselect in restrt_hdf.f)'
         STOP
      end if
c
c     # HDF: Set up dimensions for integer vector.
c
      istart(1) = 0
      iedges(1) = idims
      istride(1) = 1
c
c     # HDF: read integer vector from hdf file.
c
      istat = sfrdata(sds_id,istart,istride,iedges,iout)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to read variable ', qname,
     &        ' from restart HDF file'
         WRITE(*,*) '(call to sfrdata in restrt_hdf.f)'
         STOP
      end if
c
c     # HDF: End access to integer vector in HDF file.
c
      istat = sfendacc(sds_id)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to end access to variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfendacc in restrt_hdf.f)'
         STOP
      end if

      return
      end

c==========================================================
      subroutine read_integer_array(sd_id,idim1,idim2,index,qname,iout)
      implicit double precision (a-h,o-z)
      character*6  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id
      dimension  istart(2), istride(2), iedges(2)
      dimension  iout(idim1,idim2)
c
c     # HDF: Declare external HDF functions
c
      integer  sfcreate, sfrdata, sfselect, sfendacc
      external sfcreate, sfrdata, sfselect, sfendacc
c
c     # HDF: Set up HDF constants
c
      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
c     # HDF: Select dataset in HDF file.
c
      sds_id = sfselect(sd_id,index)
      if (sds_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to select data set for  variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfselect in restrt_hdf.f)'
         STOP
      end if
c
c     # HDF: Set up dimensions for integer array.
c
      istart(1) = 0
      istart(2) = 0
      iedges(1) = idim1
      iedges(2) = idim2
      istride(1) = 1
      istride(2) = 1
c
c     # HDF: read integer array from hdf file.
c
      istat = sfrdata(sds_id,istart,istride,iedges,iout)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to read variable ', qname,
     &        ' from restart HDF file'
         WRITE(*,*) '(call to sfrdata in restrt_hdf.f)'
         STOP
      end if
c
c     # HDF: End access to integer array in HDF file.
c
      istat = sfendacc(sds_id)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to end access to variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfendacc in restrt_hdf.f)'
         STOP
      end if

      return
      end

c==========================================================
      subroutine dump_integer_vector(sd_id,idims,qname,iout)
      implicit double precision (a-h,o-z)
      character*6  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id
      dimension  iout(idims), istart(1), istride(1), iedges(1), idim(1)
c
c     # HDF: Declare external HDF functions
c
      integer  sfcreate, sfwdata, sfscompress, sfendacc
      external sfcreate, sfwdata, sfscompress, sfendacc
c
c     # HDF: Set up HDF constants
c
      integer   DFNT_FLOAT64, DFNT_INT32
      parameter(DFNT_FLOAT64 = 6, DFNT_INT32 = 24)

      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
c     # HDF: Set up compression constants for HDF file.
c
      integer   COMP_CODE_DEFLATE, DEFLATE_LEVEL
      parameter (COMP_CODE_DEFLATE = 4, DEFLATE_LEVEL = 6)
c
c     # HDF: Create integer vector in HDF file.
c     #      NOTE THAT THE RANK MUST BE ONE.
c
      irank = 1 
      idim(1) = idims
      sds_id = sfcreate(sd_id,qname,DFNT_INT32,irank,idim)
      if (sds_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to create variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfcreate in check_hdf.f)'
         STOP
      end if
c
c     # HDF: Set up dimensions for integer vector.
c
      istart(1) = 0
      iedges(1) = idims
      istride(1) = 1
c
c     # HDF: set compression mode and write data to hdf file.
c
      istat=sfscompress(sds_id,COMP_CODE_DEFLATE,DEFLATE_LEVEL)
      istat = sfwdata(sds_id,istart,istride,iedges,iout)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to write variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfwdata in check_hdf.f)'
         STOP
      end if
c
c     # HDF: End access to integer vector in HDF file.
c
      istat = sfendacc(sds_id)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to end access to variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfendacc in check_hdf.f)'
         STOP
      end if

      return
      end

c==========================================================
      subroutine dump_integer_array(sd_id,idim1,idim2,qname,iout)
      implicit double precision (a-h,o-z)
      character*6  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id
      dimension  idims(2), istart(2), istride(2), iedges(2)
      dimension  iout(idim1,idim2)
c
c     # HDF: Declare external HDF functions
c
      integer  sfcreate, sfwdata, sfscompress, sfendacc
      external sfcreate, sfwdata, sfscompress, sfendacc
c
c     # HDF: Set up HDF constants
c
      integer   DFNT_FLOAT64, DFNT_INT32
      parameter(DFNT_FLOAT64 = 6, DFNT_INT32 = 24)

      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
c     # HDF: Set up compression constants for HDF file.
c
      integer   COMP_CODE_DEFLATE, DEFLATE_LEVEL
      parameter (COMP_CODE_DEFLATE = 4, DEFLATE_LEVEL = 6)
c
c     # HDF: Create integer array in HDF file.
c     #      NOTE THAT THE RANK MUST BE TWO.
c
      irank = 2
      idims(1) = idim1
      idims(2) = idim2
      sds_id = sfcreate(sd_id,qname,DFNT_INT32,irank,idims)
      if (sds_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to create variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfcreate in check_hdf.f)'
         STOP
      end if
c
c     # HDF: Set up dimensions for integer array.
c
      istart(1) = 0
      istart(2) = 0
      iedges(1) = idims(1)
      iedges(2) = idims(2)
      istride(1) = 1
      istride(2) = 1
c
c     # HDF: set compression mode and write data to hdf file.
c
      istat=sfscompress(sds_id,COMP_CODE_DEFLATE,DEFLATE_LEVEL)
      istat = sfwdata(sds_id,istart,istride,iedges,iout)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to write variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfwdata in check_hdf.f)'
         STOP
      end if
c
c     # HDF: End access to integer array in HDF file.
c
      istat = sfendacc(sds_id)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to end access to variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfendacc in check_hdf.f)'
         STOP
      end if

      return
      end


