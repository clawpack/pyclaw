c     
c     =====================================================
      subroutine restart(maxmx,maxmy,meqn,mbc,mx,my,
     &     xlower,ylower,dx,dy,q)
c     =====================================================
c     
c     # Initialize q using values from an old HDF output file.
c     # Copy the HDF output file to restart.data.hdf
c     # and call this routine from qinit.
c     
c     # See http://hdf.ncsa.uiuc.edu/ for more info on HDF.
c     
c     # Written 2003 by Peter Blossey 
c     
      implicit double precision (a-h,o-z)
c     
      parameter   (nDim = 2)
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      character*14 fname
c     
c     # HDF: Declare variables that describe datasets and HDF files.
c     
      integer    sd_id, sds_id, sds_start, sds_edges, sds_stride
      dimension  sds_start(nDim), sds_edges(nDim), sds_stride(nDim) 

c     # HDF: x- and y- dimensions are reversed when output to HDF file.
      dimension  qbuf(21), qout(my,mx)
c     
c     # HDF: Declare external HDF functions
c     
      integer  sfstart, sfselect, sfrdata, sfendacc, sfend
      external sfstart, sfselect, sfrdata, sfendacc, sfend
c     
c     # HDF: Set up HDF constants
c     
      integer 	DFACC_READ
      parameter(DFACC_READ = 1)

      integer   SUCCEED, FAIL
      parameter(SUCCEED = 0, FAIL = -1)
c
      common /restrt_block/ tinitial, iframe
c     
c     # first create the file name and open file
c     
      fname = 'fort.q'
     &     // char(ichar('0') + mod(iframe/1000,10)) 
     &     // char(ichar('0') + mod(iframe/100,10)) 
     &     // char(ichar('0') + mod(iframe/10,10)) 
     &     // char(ichar('0') + mod(iframe,10))
     &     // '.hdf'
      write(*,*) 'Restarting from ', fname
c     
c     # HDF: open hdf restart file.
c     
      sd_id = sfstart(fname,DFACC_READ)
      if (sd_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to open HDF file',
     &        ' (call to sfstart in restart_hdf.f)'
         STOP
      end if
c     
c     # HDF: Select grid parameter dataset in HDF file.
c     
      index = 0
      sds_id = sfselect(sd_id,index)
      if (sds_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to select data set for  variable ', index,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfselect in restrt_hdf.f)'
         STOP
      end if
c     
c     # HDF: Set up dimensions for double vector.
c     
      sds_start(1) = 0
      sds_edges(1) = 21
      sds_stride(1) = 1
c     
c     # HDF: read double vector from hdf file.
c     
      istat = sfrdata(sds_id,sds_start,sds_stride,sds_edges,qbuf)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to read variable ', index,
     &        ' from restart HDF file'
         WRITE(*,*) '(call to sfrdata in restrt_hdf.f)'
         STOP
      end if
c     
c     # HDF: End access to double vector in HDF file.
c     
      istat = sfendacc(sds_id)  
      if (istat.eq.FAIL) THEN
         WRITE(*,*) 'Failed to end access to variable ', index,
     &        ' in restart HDF file'
         WRITE(*,*) '(call to sfendacc in restrt_hdf.f)'
         STOP
      end if
c     
c     # Check parameters against those read in from claw3ez.data
c     
      mx_in = qbuf(6)
      my_in = qbuf(7)
      if (mx_in .ne. mx .or. my_in .ne. my) then
         stop 'rstart.f : grid dimensions not compatible'
      endif
c     
      meqn_in = qbuf(4)
      if (meqn_in .ne. meqn) then
         stop 'rstart.f : meqn not compatible'
      endif
c
c     # Read in new starting time.
c
      tinitial = qbuf(3)
c
c
      write(*,*) 'Restarting from file ', fname, ' at time ', tinitial
      write(*,*) 'Initial condition will not be written to a MATLAB ',
     &     'output file'
      write(*,*)
c     
c     # Loop over fields in q
c     
      do m = 1,meqn
c     
c     # HDF: Select dataset in HDF file.
c     
         index = index + 1
         sds_id = sfselect(sd_id,index)
         if (sds_id.eq.FAIL) THEN
            WRITE(*,*) 'Failed to select data set number ', index,
     &           ' in restart HDF file'
            WRITE(*,*) '(call to sfselect in restrt_hdf.f)'
            STOP
         end if
c     
c     # HDF: Set up dimensions for double array.
c     
         sds_start(1) = 0
         sds_start(2) = 0

         sds_edges(1) = my
         sds_edges(2) = mx

         sds_stride(1) = 1
         sds_stride(2) = 1
c     
c     # HDF: read double array from hdf file.
c     
         istat = sfrdata(sds_id,sds_start,sds_stride,sds_edges,qout)  
         if (istat.eq.FAIL) THEN
            WRITE(*,*) 'Failed to read variable ', index,
     &           ' from restart HDF file'
            WRITE(*,*) '(call to sfrdata in restrt_hdf.f)'
            STOP
         end if
c     
c     # HDF: End access to double array in HDF file.
c     
         istat = sfendacc(sds_id)  
         if (istat.eq.FAIL) THEN
            WRITE(*,*) 'Failed to end access to variable ', index,
     &           ' in restart HDF file'
            WRITE(*,*) '(call to sfendacc in restrt_hdf.f)'
            STOP
         end if
c     
c     # Put the data into q(:,:,m).
c
         do j = 1,my
            do i = 1,mx
               q(i,j,m) = qout(j,i)
            end do
         end do
c     
      end do
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
