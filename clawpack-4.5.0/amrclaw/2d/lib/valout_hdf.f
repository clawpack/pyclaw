c
c -----------------------------------------------------
c
      subroutine valout (lst, lend, time, nvar, naux)
c
      implicit double precision (a-h,o-z)
      include  "call.i"

      parameter   (nDim = 2)
      character*14  fname
      character*13 qname2
      character*9  qname
c
c     # HDF: Declare variables that describe datasets and HDF files.
c
      integer    sd_id, sds_id, sds_rank
      integer    sds_dims, sds_start, sds_edges, sds_stride
      dimension  sds_dims(nDim), sds_start(nDim)
      dimension  sds_edges(nDim), sds_stride(nDim) 

c     # HDF: Arrays which will hold output arrays.
      integer    ntotal_output
      parameter (ntotal_output = 1000000)
      dimension  qbuf(21), qout(ntotal_output)
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
      iadd(i,j,ivar) = loc + i - 1 + mitot*((ivar-1)*mjtot+j-1)
c
c ::::::::::::::::::::::::::: VALOUT ::::::::::::::::::::::::::::::::::;
c valout = graphics output of soln values for contour or surface plots.
c          can output for matlab or ncar graphics post-processing
c :::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::::::::::;
c

c     ### NCAR graphics output

      if (ncarout) then

        call basic (time, lst, lend )
c
        write(pltunit1,100)  nvar
100     format(10h*VALS     ,i10)
c
        level = lst
10      if (level .gt. lend) go to 60
            mptr = lstart(level)
20          if (mptr .eq. 0) go to 50
                nx = node(ndihi,mptr)-node(ndilo,mptr) + 1
                ny = node(ndjhi,mptr)-node(ndjlo,mptr) + 1
                mitot = nx + 2*nghost
                mjtot = ny + 2*nghost
                loc = node(store1,mptr)
                call outvar(alloc(loc),mitot,mjtot,nvar,mptr,nghost)
                mptr = node(levelptr,mptr)
            go to 20
50          continue
            level = level + 1
        go to 10
c
60    continue
      endif


c     ### MATLAB graphics output
c
c     # HDF: Modified 8/2003 by Peter Blossey to generate HDF
c     #      (version 4) output files.  These portable, binary
c     #      format output files allow for self-describing data
c     #      sets as well as data compression.
c
c     #      For more information about HDF, see
c     #      http://hdf.ncsa.uiuc.edu/
c
      if (matlabout) then
c        ###  make the file names and open output files
         iframe = matlabu
         fname = 'fort.q'
     &        // char(ichar('0') + mod(iframe/1000,10)) 
     &        // char(ichar('0') + mod(iframe/100,10)) 
     &        // char(ichar('0') + mod(iframe/10,10)) 
     &        // char(ichar('0') + mod(iframe,10))
     &        // '.hdf'

         level = lst
         ngrids = 0
c
c     # HDF: create hdf file.
c
         sd_id = sfstart(fname, DFACC_CREATE)
         if (sd_id.eq.FAIL) THEN
            WRITE(*,*) 'Failed to create HDF file (call to sfstart)'
            STOP
         end if

 65      if (level .gt. lfine) go to 90
         mptr = lstart(level)
 70      if (mptr .eq. 0) go to 80
         ngrids  = ngrids + 1
         nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         loc     = node(store1, mptr)
         locaux  = node(storeaux,mptr)
         mitot   = nx + 2*nghost
         mjtot   = ny + 2*nghost
         xlow = rnode(cornxlo,mptr)
         ylow = rnode(cornylo,mptr)

c     
c     # HDF: Check whether grid is larger than statically allocated
c     #      output array.
c     
         if ((nx*ny).gt.ntotal_output) then
            WRITE(*,*) 'Grid number ', mptr, ' on level number ',
     &           level, ' exceeds ntotal_output'
            WRITE(*,*) 'Increase ntotal_output in valout_hdf.f'
            STOP 
         end if
c     
c     # HDF: create a data set for parameters describing q in HDF file.
c     
         qname = 'grid_'
     &        // char(ichar('0') + mod(mptr/1000,10)) 
     &        // char(ichar('0') + mod(mptr/100,10)) 
     &        // char(ichar('0') + mod(mptr/10,10)) 
     &        // char(ichar('0') + mod(mptr,10))
         
         sds_rank = 1
         sds_dims(1) = 21
         
         sds_id = sfcreate(sd_id,qname,DFNT_FLOAT64,sds_rank,sds_dims)
         if (sds_id.eq.FAIL) THEN
            WRITE(*,*) 
     &           'Failed to create scientific data set in HDF file'
            STOP
         end if
c     
c     # HDF: set up parameters describing data set.
c     
         sds_start(1)  = 0
         sds_edges(1)  = sds_dims(1)
         sds_stride(1) = 1        
         qbuf(1) = ngrids
         qbuf(2) = nDim
         qbuf(3) = time
         qbuf(4) = nvar
         qbuf(5) = level
         qbuf(6) = nx
         qbuf(7) = ny
         qbuf(8) = 0.
         qbuf(9) = 0.
         qbuf(10) = xlow
         qbuf(11) = ylow
         qbuf(12) = 0.
         qbuf(13) = 0.
         qbuf(14) = xlow+nx*hxposs(level)
         qbuf(15) = ylow+ny*hyposs(level)
         qbuf(16) = 0.
         qbuf(17) = 0.
         qbuf(18) = hxposs(level)
         qbuf(19) = hyposs(level)
         qbuf(20) = 0.
         qbuf(21) = 0.
         istat = sfwdata(sds_id,sds_start,sds_stride,sds_edges,qbuf)  
         istat = sfendacc(sds_id)  
c     
c     # Loop over fields in q
c     
         do ivar = 1,nvar
c     
c     # HDF: create a data set for parameters describing q in HDF file.
c     
            qname2 = qname // '_'
     &           // char(ichar('0') + mod(ivar/100,10)) 
     &           // char(ichar('0') + mod(ivar/10,10)) 
     &           // char(ichar('0') + mod(ivar,10))
c     
c     # HDF: Reverse dimensions when storing arrays because of different
c     #      conventions between c and FORTRAN as to which dimension should
c     #      be stored first.  Reversing the dimensions here will make
c     #      the x-direction first when reading into MATLAB.
c     
            sds_rank = nDim
            sds_dims(1) = ny
            sds_dims(2) = nx
c     
            sds_id = sfcreate(sd_id,qname2,DFNT_FLOAT64,sds_rank,
     &           sds_dims)
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

            sds_start(2)  = 0
            sds_edges(2)  = sds_dims(2)
            sds_stride(2) = 1        
c     
c     # Copy current field of q into qout.
c     
            ioffset = 0
            do j = nghost+1, mitot-nghost
               do i = nghost+1, mjtot-nghost
                  qout(ioffset+i-nghost) = alloc(iadd(j,i,ivar))
               end do
               ioffset = ioffset + ny
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

         mptr = node(levelptr, mptr)
         go to 70
 80      level = level + 1
         go to 65

 90      continue

c     
c     # HDF: Close HDF file.
c     
         istat = sfend(sd_id)
         if (istat.eq.FAIL) then
            WRITE(*,*) 'Failed to close ', fname, ' (call to sfend)'
            STOP
         end if
c     
c     # Output message saying that write is complete.
c     
         write(6,601) matlabu,time
 601     format('AMRCLAW: Frame ',i4,
     &        ' output files done at time t = ', d12.5,/)

         matlabu = matlabu + 1
      endif  

      return
      end
