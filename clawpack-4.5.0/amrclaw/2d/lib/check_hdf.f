c
c ---------------------------------------------------------
c
      subroutine check(nsteps,time,nvar,naux)
c
c :::::::::::::::::::::: CHECK ::::::::::::::::::::::::::::::::;
c   check point routine - can only call at end of coarse grid cycle
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      implicit double precision (a-h,o-z)
      include  "call.i"

      character*16  chkname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id
c
c     # HDF: Arrays which will hold output arrays.
      dimension  iqout(14), qout(4)
c
c     # HDF: Declare external HDF functions
c
      integer  sfstart, sfend
      external sfstart, sfend
c
c     # HDF: Set up HDF constants
c
      integer 	DFACC_CREATE
      parameter(DFACC_CREATE = 4)

      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
c     ###  make the file name showing the time step
c
      chkname = 'fort.chk'
     &     // char(ichar('0') + mod(nsteps/1000,10)) 
     &     // char(ichar('0') + mod(nsteps/100,10)) 
     &     // char(ichar('0') + mod(nsteps/10,10)) 
     &     // char(ichar('0') + mod(nsteps,10))
     &     // '.hdf'
c     
c     # HDF: create hdf file.
c     
      sd_id = sfstart(chkname, DFACC_CREATE)
      if (sd_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to create HDF file (call to sfstart)'
         STOP
      end if
c     
c     # HDF: create array of integers for storage.
c     
c     # HDF: DFNT_INT32
      iqout(1) = lenmax
      iqout(2) = lendim
      iqout(3) = memsize
      iqout(4) = lenf
      iqout(5) = ibuff
      iqout(6) = mstart
      iqout(7) = ndfree
      iqout(8) = lfine
      iqout(9) = iorder
      iqout(10) = mxnest
      iqout(11) = kcheck
      iqout(12) = nsteps
      iqout(13) = matlabu
      iqout(14) = lentot
      call dump_integer_vector(sd_id,14,'iqout ',iqout)
      call dump_integer_vector(sd_id,maxlv,'icheck',icheck)
      call dump_integer_vector(sd_id,maxlv,'lstart',lstart)
      call dump_integer_vector(sd_id,maxlv,'newstl',newstl)
      call dump_integer_vector(sd_id,maxlv,'listsp',listsp)
      call dump_integer_vector(sd_id,maxlv,'intratx',intratx)
      call dump_integer_vector(sd_id,maxlv,'intraty',intraty)
      call dump_integer_vector(sd_id,maxlv,'kratio',kratio)
      call dump_integer_vector(sd_id,maxlv,'iregsz',iregsz)
      call dump_integer_vector(sd_id,maxlv,'jregsz',jregsz)

      call dump_integer_array(sd_id,lfdim,2,'lfree ',lfree)
      call dump_integer_array(sd_id,rsize,maxgr,'node  ',node)

c     # HDF: DFNT_FLOAT64

      qout(1) = tol
      qout(2) = time
      qout(3) = evol
      qout(4) = rvol
      call dump_double_vector(sd_id,4,'qout  ',qout)
      call dump_double_vector(sd_id,maxlv,'hxposs',hxposs)
      call dump_double_vector(sd_id,maxlv,'hyposs',hyposs)
      call dump_double_vector(sd_id,maxlv,'possk ',possk)
      call dump_double_vector(sd_id,10,'rvoll ',rvoll)
      call dump_double_vector(sd_id,lendim,'alloc ',alloc)
      
      call dump_double_array(sd_id,rsize,maxgr,'rnode ',rnode)
c
c     # HDF: Close HDF file.
c
      istat = sfend(sd_id)
      if (istat.eq.FAIL) then
         WRITE(*,*) 'Failed to close SDS (call to sfend)'
         STOP
      end if

      return
      end

c=================================================================
c     # HDF: NOTE THAT ALL DOUBLE PRECISION HDF CHECKPOINTING ROUTINES
c     #      FOR BOTH WRITING AND READING ARE INCLUDED HERE, ALTHOUGH
c     #      THE DOUBLE PRECISION READING ROUTINES ARE CALLED IN 
c     #      restrt_hdf.f.  THIS IS DONE TO PREVENT COMPILER ERRORS 
c     #      SINCE THE SAME HDF ROUTINES ARE USED TO READ BOTH INTEGER
c     #      AND DOUBLE PRECISION ARRAYS.
c=================================================================
      subroutine dump_double_vector(sd_id,idims,qname,out)
      implicit double precision (a-h,o-z)
      character*6  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id
      dimension  out(idims), istart(1), istride(1), iedges(1), idim(1)
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
c     # HDF: Create double vector in HDF file.
c     #      NOTE THAT THE RANK MUST BE ONE.
c
      irank = 1
      idim(1) = idims
      sds_id = sfcreate(sd_id,qname,DFNT_FLOAT64,irank,idim)
      if (sds_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to create variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfcreate in check_hdf.f)'
         STOP
      end if
c
c     # HDF: Set up dimensions for double vector.
c
      istart(1) = 0
      iedges(1) = idims
      istride(1) = 1
c
c     # HDF: set compression mode and write data to hdf file.
c
      istat=sfscompress(sds_id,COMP_CODE_DEFLATE,DEFLATE_LEVEL)
      istat = sfwdata(sds_id,istart,istride,iedges,out)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to write variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfwdata in check_hdf.f)'
         STOP
      end if
c
c     # HDF: End access to double vector in HDF file.
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
      subroutine dump_double_array(sd_id,idim1,idim2,qname,out)
      implicit double precision (a-h,o-z)
      character*6  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id
      dimension  idims(2), istart(2), istride(2), iedges(2)
      dimension  out(idim1,idim2)
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
c     # HDF: Create double array in HDF file.
c     #      NOTE THAT THE RANK MUST BE TWO.
c
      irank = 2
      idims(1) = idim1
      idims(2) = idim2
      sds_id = sfcreate(sd_id,qname,DFNT_FLOAT64,irank,idims)
      if (sds_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to create variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfcreate in check_hdf.f)'
         STOP
      end if
c
c     # HDF: Set up dimensions for double array.
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
      istat = sfwdata(sds_id,istart,istride,iedges,out)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to write variable ', qname,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfwdata in check_hdf.f)'
         STOP
      end if
c
c     # HDF: End access to double array in HDF file.
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
      subroutine read_double_vector(sd_id,idims,index,qname,out)
      implicit double precision (a-h,o-z)
      character*6  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id
      dimension  out(idims), istart(1), istride(1), iedges(1)
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
c     # HDF: Set up dimensions for double vector.
c
      istart(1) = 0
      iedges(1) = idims
      istride(1) = 1
c
c     # HDF: read double vector from hdf file.
c
      istat = sfrdata(sds_id,istart,istride,iedges,out)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to read variable ', qname,
     &        ' from restart HDF file'
         WRITE(*,*) '(call to sfrdata in restrt_hdf.f)'
         STOP
      end if
c
c     # HDF: End access to double vector in HDF file.
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
      subroutine read_double_array(sd_id,idim1,idim2,index,qname,out)
      implicit double precision (a-h,o-z)
      character*6  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id
      dimension  istart(2), istride(2), iedges(2)
      dimension  out(idim1,idim2)
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
c     # HDF: Set up dimensions for double array.
c
      istart(1) = 0
      istart(2) = 0
      iedges(1) = idim1
      iedges(2) = idim2
      istride(1) = 1
      istride(2) = 1
c
c     # HDF: read double array from hdf file.
c
      istat = sfrdata(sds_id,istart,istride,iedges,out)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to read variable ', qname,
     &        ' from restart HDF file'
         WRITE(*,*) '(call to sfrdata in restrt_hdf.f)'
         STOP
      end if
c
c     # HDF: End access to double array in HDF file.
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
