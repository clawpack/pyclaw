c
c ----------------------------------------------------------------
c
       subroutine auxcoarsen(auxdub,midub,mjdub,auxbgc,
     1                       mi2tot,mj2tot,naux,auxtype)
      
       implicit double precision (a-h, o-z)

       dimension     auxdub(midub, mjdub, naux)
       dimension     auxbgc(mi2tot,mj2tot,naux)
       character*10  auxtype(naux)

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
c coarsen = coarsen the fine grid auxiliary data (with double the usual
c           number of ghost cells to prepare coarsened data
c           for error estimation.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       do 50 iaux = 1, naux

       if (auxtype(iaux) .eq. "center" .or. 
     .     auxtype(iaux) .eq. "capacity") then
         do 20 j = 1, mj2tot
            jfine = 2*(j-1) + 1
            do 20 i = 1, mi2tot
               ifine = 2*(i-1) + 1
               auxbgc(i,j,iaux) = (auxdub(ifine,jfine,iaux) +
     &                             auxdub(ifine+1,jfine,iaux)+
     &                             auxdub(ifine,jfine+1,iaux) +
     &                             auxdub(ifine+1,jfine+1,iaux))/4.d0
20       continue

       elseif (auxtype(iaux) .eq. "xleft") then 
         do 10 j = 1, mj2tot
            jfine = 2*(j-1) + 1
            do 10 i = 1, mi2tot
               ifine = 2*(i-1) + 1
               auxbgc(i,j,iaux) = (auxdub(ifine,jfine,iaux) +
     &                          auxdub(ifine,jfine+1,iaux)) /2.d0 
10       continue

       elseif (auxtype(iaux) .eq. "yleft") then 
         do 15 j = 1, mj2tot
            jfine = 2*(j-1) + 1
            do 15 i = 1, mi2tot
               ifine = 2*(i-1) + 1
               auxbgc(i,j,iaux) = (auxdub(ifine,jfine,iaux) +
     &                             auxdub(ifine+1,jfine,iaux))/2.d0
15       continue

       endif

50     continue

       return
       end
