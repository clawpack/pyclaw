c
c
c =========================================================
      subroutine out1(maxmx,meqn,mbc,mx,xlower,dx,q,t,iframe,aux,maux)
c =========================================================
c
c
c     # Output the results for a general system of conservation laws
c     # in 1 dimension to a hdf file as a scientific data set.
c     # See http://hdf.ncsa.uiuc.edu/ for more info on HDF.
c     # The results in this output file can be plotted in MATLAB
c     # using the "plotclaw1" script.
c
c     # Revised 2003 by Peter Blossey
c     # Adapted code written by Sorin Mitran to standard F77 CLAWPACK style
c
      implicit double precision (a-h,o-z)
c
      parameter   (nDim = 1)
      dimension    q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, maux)
      character*14 fname
      character*13 qname2
      character*9  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id, sds_rank
      integer    sds_dims, sds_start, sds_edges, sds_stride
      dimension  sds_dims(nDim), sds_start(nDim)
      dimension  sds_edges(nDim), sds_stride(nDim) 
      dimension  qbuf(21), qout(mx)
c
c     # HDF: Declare external HDF functions
c
      integer  sfstart, sfcreate, sfwdata, sfscompress, sfendacc, sfend
      external sfstart, sfcreate, sfwdata, sfscompress, sfendacc, sfend
c
c     # HDF: Set up HDF constants
c
      integer 	DFACC_READ, DFACC_WRITE, DFACC_CREATE 
      parameter(DFACC_READ = 1, DFACC_WRITE = 2, DFACC_CREATE = 4)

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
c     # Write the results to the file fort.q<iframe>
c     # Use format required by matlab script  plotclaw2.m
c     # The same format is used by the amrclaw package.
c     # Here it's adapted to output just the single grid.
c
c     # first create the file name and open file
c
      fname = 'fort.q'
     &     // char(ichar('0') + mod(iframe/1000,10)) 
     &     // char(ichar('0') + mod(iframe/100,10)) 
     &     // char(ichar('0') + mod(iframe/10,10)) 
     &     // char(ichar('0') + mod(iframe,10))
     &     // '.hdf'
c
c     # Specify grid number and create a string which will describe this
c     # grid in HDF file.  This could be an input for simulations with
c     # multiple grids, as in AMRCLAW.  
c
      ngrids_out=1
      qname = 'grid_'
     &     // char(ichar('0') + mod(ngrids_out/1000,10)) 
     &     // char(ichar('0') + mod(ngrids_out/100,10)) 
     &     // char(ichar('0') + mod(ngrids_out/10,10)) 
     &     // char(ichar('0') + mod(ngrids_out,10))
c
c     # HDF: create hdf file.
c
      sd_id = sfstart(fname, DFACC_CREATE)
      if (sd_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to create HDF file (call to sfstart)'
         STOP
      end if
c
c     # HDF: create a data set for parameters describing q in HDF file.
c
      sds_rank = 1
      sds_dims(1) = 21
      
      sds_id = sfcreate(sd_id,qname,DFNT_FLOAT64,sds_rank,sds_dims)
      if (sds_id.eq.FAIL) THEN
         WRITE(*,*) 'Failed to create scientific data set in HDF file'
         STOP
      end if
c
c     # HDF: set up parameters describing data set.
c
      sds_start(1)  = 0
      sds_edges(1)  = sds_dims(1)
      sds_stride(1) = 1        
      qbuf(1) = ngrids_out
      qbuf(2) = nDim
      qbuf(3) = t
      qbuf(4) = meqn
      qbuf(5) = 1.
      qbuf(6) = mx
      qbuf(7) = 0.
      qbuf(8) = 0.
      qbuf(9) = 0.
      qbuf(10) = xlower
      qbuf(11) = 0.
      qbuf(12) = 0.
      qbuf(13) = 0.
      qbuf(14) = xlower+mx*dx
      qbuf(15) = 0.
      qbuf(16) = 0.
      qbuf(17) = 0.
      qbuf(18) = dx
      qbuf(19) = 0.
      qbuf(20) = 0.
      qbuf(21) = 0.
      istat = sfwdata(sds_id,sds_start,sds_stride,sds_edges,qbuf)  
      istat = sfendacc(sds_id)  
c
c     # Loop over fields in q
c
      do m = 1,meqn
c
c     # HDF: create a data set for parameters describing q in HDF file.
c
         qname2 = qname // '_'
     &        // char(ichar('0') + mod(m/100,10)) 
     &        // char(ichar('0') + mod(m/10,10)) 
     &        // char(ichar('0') + mod(m,10))
c
         sds_rank = nDim
         sds_dims(1) = mx
c      
         sds_id = sfcreate(sd_id,qname,DFNT_FLOAT64,sds_rank,sds_dims)
         if (sds_id.eq.FAIL) THEN
            WRITE(*,*) 'Failed to create data set in HDF file'
            STOP
         end if
c
c     # HDF: set up parameters describing data set.
c
         sds_start(1)  = 0
         sds_edges(1)  = sds_dims(1)
         sds_stride(1) = 1        

c     # Copy current field of q into qout.
c
         do i = 1,mx
            qout(i) = q(i,m)
         end do
c
c     # HDF: set compression mode and write data to hdf file.
c
         istat=sfscompress(sds_id,COMP_CODE_DEFLATE,DEFLATE_LEVEL)
         istat = sfwdata(sds_id,sds_start,sds_stride,sds_edges,qout)
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to write SDS (call to sfwdata)'
            STOP
         end if
c
c     # HDF: Close the data set
c
         istat = sfendacc(sds_id)  
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to close SDS (call to sfendacc)'
            STOP
         end if
      end do
c
c     # HDF: Close HDF file.
c
      istat = sfend(sd_id)
      if (istat.eq.FAIL) then
         WRITE(*,*) 'Failed to close ', fname, ' (call to sfend)'
         STOP
      end if

      return
      end
