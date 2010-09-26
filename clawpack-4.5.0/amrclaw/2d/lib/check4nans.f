c
c
c =========================================================
      subroutine check4nans(maxmx,maxmy,meqn,mbc,mx,my,q,t,
     &                      ichecknan)
c =========================================================
c
c     # check for NANs in solution q
c
      implicit double precision (a-h,o-z)
      dimension   q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)

c     write(*,*) 'Checking for NANs at ichecknan = ',ichecknan
c     write(*,*) '  maxmx = ',maxmx,'  maxmy = ',maxmy,'  meqn = ',meqn
      
      do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
          do m=1,meqn
            if (.not.(q(i,j,m).eq.q(i,j,m))) then
c             # true if q(i,j,m) = NAN
	      write(*,*) 'SOLUTION ERROR --- ABORTING CALCULATION'
	      write(*,*) 'At ichecknan = ',ichecknan
              write(*,*) '   mx,my,t:',mx,my,t
              write(*,*) '   i,j,m:',i,j,m
              write(*,*) '   q(i,j,m) = ',q(i,j,m)
	      stop
              endif
            enddo
          enddo
        enddo
c
c     # uncomment the next line if desired when debugging:
c     write(*,*) 'No NANs at ichecknan = ',ichecknan,' at t = ',t

      return
      end
