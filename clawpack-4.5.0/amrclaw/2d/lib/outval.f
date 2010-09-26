c
c =======================================================================
      subroutine outval(val,nvar,mitot,mjtot,mptr,outgrd,naux,aux)
c =======================================================================
c
      implicit double precision (a-h,o-z)

      dimension  val(mitot,mjtot,nvar)
      dimension  aux(mitot,mjtot,naux)
      logical    outgrd

      include  "call.i"

c ::::::::::::::::::::::OUTVAL :::::::::::::::::::::::::::::::
c print solution and aux. variables to output. 
c if cell outside domain, don't print soln. value - nothing
c currently in ghost cells.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      if (.not. outgrd) go to 99
      level = node(nestlevel,mptr)
      hx    = hxposs(level)
      hy    = hyposs(level)
      cornx = rnode(cornxlo,mptr) -  nghost*hx
      corny = rnode(cornylo,mptr) -  nghost*hy
c
      do 25 j=nghost+1,mjtot-nghost
      do 20 i=nghost+1,mitot-nghost

          x  = cornx + hx*(dble(i)-.5d0)
          y  = corny + hy*(dble(j)-.5d0)
          write(outunit,107) x,y,i,j,(val(i,j,ivar),ivar=1,nvar)
 107      format(2hx=,f6.3,3h y=,f5.3,3h,i=,i3,3h,j=,i3,' a=',
     *           5(e9.3,1x))
          if (naux.gt.0) write(outunit,108) (aux(i,j,iaux),iaux=1,naux)
 108      format(1x,'aux = ',7(e9.3,1x))

 20   continue
 25   continue

 99   return
      end
