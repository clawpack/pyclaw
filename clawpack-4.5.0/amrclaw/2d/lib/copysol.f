c
c ----------------------------------------------------------
c
      subroutine copysol(valbig,val,nvar,mitot,mjtot,nghost,
     1                   midub,mjdub,ngbig)
c
      implicit double precision (a-h,o-z)

      dimension  valbig(midub,mjdub,nvar), val(mitot,mjtot,nvar)
c
c copy solution into grid with different number ghsot cells
c
       do 10 ivar = 1, nvar
       do 10 j = nghost+1, mjtot-nghost
       do 10 i = nghost+1, mitot-nghost
          valbig(i-nghost+ngbig,j-nghost+ngbig,ivar) = val(i,j,ivar)
 10    continue
c
       return
       end
