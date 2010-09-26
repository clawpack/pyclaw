


c
c ---------------------------------------------------
c
       subroutine coarsen(valdub,midub,mjdub,valbgc,mi2tot,mj2tot,nvar)
      
       implicit double precision (a-h, o-z)

       dimension  valdub(midub, mjdub, nvar)
       dimension  valbgc(mi2tot,mj2tot,nvar)

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
c coarsen = coarsen the fine grid data (with double the usual
c           number of ghost cells to prepare coarsened
c           grid for error estimation.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       do 10 ivar = 1, nvar
       do 10 j = 1, mj2tot

          jfine = 2*(j-1) + 1

          do 10 i = 1, mi2tot

             ifine = 2*(i-1) + 1
             valbgc(i,j,ivar) = (valdub(ifine,jfine,ivar) +
     &                           valdub(ifine+1,jfine,ivar)+
     &                           valdub(ifine,jfine+1,ivar) +
     &                           valdub(ifine+1,jfine+1,ivar))/4.d0
10     continue

       return
       end
