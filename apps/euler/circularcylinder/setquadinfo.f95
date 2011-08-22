      subroutine setquadinfo(maxmx, maxmy, mbc,mx,my,xlower,ylower,
     &      dx,dy,maxlevel,level,refratios,mcapa,mbc_mapped, maux,aux)

      implicit none

      integer maxmx, maxmy, mbc, mx, my, maux ,mcapa
      double precision xlower, ylower, dx, dy
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)
      double precision normals(2,2), quad(0:1,0:1,2), slength(2)

      double precision sum_length(2)
      double precision dxf, dyf, xef, yef, xe, ye, xp, yp
      double precision sum_area, xcorner, ycorner, area
      integer rfactor, i, j, ir, jr, icell, jcell
      integer maxlevel, level, mbc_mapped(4)

c     # 'Refratios' should have at least 'maxlevel' entries.  For
c     # CLAWPACK, this array is allocated in driver.f.  For
c     # ChomboClaw, this is allocated internally.  This
c     # array has not been allocatd for AMRClaw.
      integer refratios(maxlevel), iref

      if (maux .lt. mcapa) then
         write(6,*) 'setquadinfo : maux must be >= 7 when using '
         write(6,*) 'mapped grids with routine setquadinfo.f'
         stop
      endif
      if (mcapa .ne. 7) then
         write(6,'(A,A,I2,A)') 'setquadinfo : mcapa should ',
     &         'be set to 7 when using setquadinfo.'
         stop
      endif


c     # Compute the refinement factor for this level. Note that
c     # the values assigned to 'level' and 'maxlevel', depend on
c     # how the value for the coarsest level is
c     # defined.  For AMRClaw, the coarsest level is level '1'.
c     # For ChomboClaw, this level is level '0'.  For this reason,
c     # we don't want to use the variable 'level' as an index
c     # into the array 'refratios'.  We use 'iref' instead.
      rfactor = 1
      iref = 1
      do ir = level,maxlevel-1
         rfactor = rfactor*refratios(iref)
         iref = iref + 1
      enddo
      dxf = dx/rfactor
      dyf = dy/rfactor

c     # Note that we allow the user to indicate how far into the ghost
c     # cell region to call the mapping.  For some mappings (the polar
c     # grid for example), it is difficult to define the mapping
c     # in a meaningful way for the ghost cell values
c     # corresponding to r < 0.
      do j = 1-mbc_mapped(3),my+mbc_mapped(4)
         ye = ylower + (j-1)*dy
         do i = 1-mbc_mapped(1),mx+mbc_mapped(2)
            xe = xlower + (i-1)*dx

            sum_area = 0.d0
            sum_length(1) = 0.d0
            sum_length(2) = 0.d0
            do ir = 1,rfactor
               xef = xe + (ir-1)*dxf
               do jr = 1,rfactor
                  yef = ye + (jr - 1)*dyf
                  do icell = 0,1
                     xcorner = xef + icell*dxf
                     do jcell = 0,1
                        ycorner = yef + jcell*dyf
                        call mapc2p(xcorner,ycorner,xp,yp)
                        quad(icell,jcell,1) = xp
                        quad(icell,jcell,2) = yp
                     enddo
                  enddo

c                 # Compute area, lengths, and normals for the
c                 # cell whose corners are stored in 'quad'.
                  call compute_info2(quad,normals,slength,area)

c                 # Accumulate areas and lengths;  normal information
c                 # isn't used here, but recomputed for coarse cell
c                 # below.
                  sum_area = sum_area + area
                  if (ir .eq. 1) then
                     sum_length(1) = sum_length(1) + slength(1)
                  endif
                  if (jr .eq. 1) then
                     sum_length(2) = sum_length(2) + slength(2)
                  endif
               enddo
            enddo

c           # Now get normal vector at interfaces
c           # Use coarse grid straight edge approximation
            do icell = 0,1
               xcorner = xe + icell*dx
               do jcell = 0,1
                  ycorner = ye + jcell*dy
                  call mapc2p(xcorner,ycorner,xp,yp)
                  quad(icell,jcell,1) = xp
                  quad(icell,jcell,2) = yp
               enddo
            enddo

c          # Call same routine as above, this time only using
c           # normals information
           call compute_info2(quad,normals,slength,area)

c          # Store information in aux array.  Of course, this can be
c          # stored in any arrangement desired.
c          # Normal at x face - first row in normals matrix
           aux(i,j,1) = normals(1,1)
           aux(i,j,2) = normals(1,2)
           aux(i,j,3) = sum_length(1)/dy

c          # Normal at y face - second row in normals matrix.
           aux(i,j,4) = normals(2,1)
           aux(i,j,5) = normals(2,2)
           aux(i,j,6) = sum_length(2)/dx

c          Capacity
           if (sum_area .le. 0.d0) then
              write(6,*) 'setquadinfo : capacity <= 0'
              write(6,*) 'area = ', sum_area
              write(6,*) 'i,j = ', i,j
              write(6,'(A,A,A,A)') 'This usually means that there ',
     &              'is a problem with your mapping function. ',
     &              'Check to be sure your mapping is valid in ',
     &              'ghost cell regions. '
              write(6,*) 'Area is being to set abs(area)'
              sum_area = abs(sum_area)
c              stop
           endif
           aux(i,j,mcapa) = sum_area/(dx*dy)
        enddo
      enddo

      end

c     # This is where the area, length and normals are computed
      subroutine compute_info2(quad,normals,slength,area)
      implicit none

      double precision normals(2,2),quad(0:1,0:1,2), slength(2)
      double precision v(2), vnorm
      double precision w(2), wnorm
      double precision area, xpcorn(5), ypcorn(5)
      integer ic

c     # compute normals to left edge
      v(1) = quad(0,1,1) - quad(0,0,1)
      v(2) = quad(0,1,2) - quad(0,0,2)
      vnorm = sqrt(v(1)*v(1) + v(2)*v(2))
      slength(1) = vnorm
      if (vnorm .eq. 0.d0) then
         normals(1,1) = 1.d0
         normals(1,2) = 0.d0
      else
         normals(1,1) = v(2)/vnorm  !! Normal should point in positive x dir.
         normals(1,2) = -v(1)/vnorm
      endif
c
      w(1) = quad(1,0,1) - quad(0,0,1)
      w(2) = quad(1,0,2) - quad(0,0,2)
      wnorm = sqrt(w(1)*w(1) + w(2)*w(2))
      slength(2) = wnorm
      if (wnorm .eq. 0.d0) then
         normals(2,1) = 1.d0
         normals(2,2) = 0.d0
      else
         normals(2,1) = -w(2)/wnorm  !! normal should point in positive y dir.
         normals(2,2) = w(1)/wnorm
      endif

      xpcorn(1) = quad(0,0,1)
      xpcorn(2) = quad(0,1,1)
      xpcorn(3) = quad(1,1,1)
      xpcorn(4) = quad(1,0,1)
      xpcorn(5) = quad(0,0,1)

      ypcorn(1) = quad(0,0,2)
      ypcorn(2) = quad(0,1,2)
      ypcorn(3) = quad(1,1,2)
      ypcorn(4) = quad(1,0,2)
      ypcorn(5) = quad(0,0,2)

      area = 0.d0
      do ic=1,4
         area = area + 0.5d0 * (ypcorn(ic)+ypcorn(ic+1)) *
     &         (xpcorn(ic+1)-xpcorn(ic))
      enddo

      end
