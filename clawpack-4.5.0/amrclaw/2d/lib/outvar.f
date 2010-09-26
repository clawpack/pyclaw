c
c -------------------------------------------------------------
c
      subroutine outvar(rect,mitot,mjtot,nvar,mptr,ng)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      dimension rect(mitot,mjtot,nvar)

c ::::::::::::::: OUTVAR ::::::::::::::::::::::::::::::::::
c
c  dump soln for graphics 
c
c  only output max - 1 rows and cols, since with cell centered
c  variables there is one extra cell outside the grid.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      write(pltunit1,100) mptr
 100  format('*SOLN     ',i10,' is the grid - all variables')
c
      do 20 ivar = 1, nvar
         write(pltunit1,101) ((rect(i,j,ivar),i=ng+1,mitot-ng),
     .                                 j=ng+1,mjtot-ng)
 101     format(5e12.6)
 20   continue
c
      return
      end
