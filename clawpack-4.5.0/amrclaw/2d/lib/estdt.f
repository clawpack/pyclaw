c
c-----------------------------------------------------------------------
c
      subroutine estdt(val,mitot,mjtot,nvar,dx,dy,dt,nghost,aux,
     &                 naux,cfl)
c
c :::::::::::::::::::::::: ESTDT :::::::::::::::::::::::::::;
c  estimate the initial time step for the given values
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

       implicit double precision (a-h, o-z)
       dimension val(mitot,mjtot,nvar)
       dimension aux(mitot,mjtot,naux)
c
c
       return
       end
