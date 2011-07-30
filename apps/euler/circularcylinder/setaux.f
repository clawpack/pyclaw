c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
      implicit none

      integer maxmx, maxmy, mbc, mx,my, maux
      double precision xlower, ylower, dx,dy
      integer maxlevel, level, refratios(20), mcapa
      integer i,j, ij

      double precision kappa_max, kappa_min, kappa_avg
      integer mapped_mbc(4)

      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,maux)

      common /comlevel/ maxlevel, level, refratios

      mcapa = 7

c     # we can call the mapping for the ghost cells.
      do i = 1,4
         mapped_mbc(i) = 2
      enddo

      call setquadinfo(maxmx, maxmy, mbc,mx,my,xlower,ylower,
     &      dx,dy,maxlevel,level,refratios,mcapa,
     &      mapped_mbc, maux, aux)

      kappa_max = 0
      kappa_min = 100
      kappa_avg = 0;
      do j = 1,my
         do i = 1,mx
            kappa_max = max(aux(i,j,mcapa),kappa_max)
            kappa_min = min(aux(i,j,mcapa),kappa_min)
            kappa_avg = kappa_avg + aux(i,j,mcapa)
         enddo
      enddo

      write(6,'(A,F16.8)') 'Max kappa : ', kappa_max
      write(6,'(A,F16.8)') 'Min kappa : ', kappa_min
      write(6,'(A,F16.8)') 'Avg kappa : ', kappa_avg/(mx*my)
      write(6,'(A,F16.8)') 'Ratio     : ', kappa_max/kappa_min
      write(6,*) ' '


      end
