c
        subroutine fixcapaq(val,aux,mitot,mjtot,valc,auxc,mic,mjc,
     &                      nvar,naux,levc)

      implicit double precision (a-h,o-z)

      include "call.i"
c
c :::::::::::::::::::::::  FIXCAPAQ ::::::::::::::::::::::::::::::
c  new fine grid solution q was linearly interpolated. but want
c  to conserve kappa*q, not q. calculate the discrepancy
c  in kappa*q using this q, and modify q to account for it.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      dimension   val(mitot,mjtot,nvar), valc(mic,mjc,nvar)
      dimension   aux(mitot,mjtot,naux), auxc(mic,mjc,naux)

      dcapamax = 0.d0
      lratiox  = intratx(levc)
      lratioy  = intraty(levc)

      do 10 ic = 2, mic-1
      do 10 jc = 2, mjc-1


       do 15 ivar = 1, nvar

       capaqfine = 0.d0

       do 20 ico = 1, lratiox
       ifine = (ic-2)*lratiox + nghost + ico
       do 20 jco = 1, lratioy
         jfine = (jc-2)*lratioy + nghost + jco
         capaqfine = capaqfine + aux(ifine,jfine,mcapa)*
     &                           val(ifine,jfine,ivar)
20     continue

       dcapaq = auxc(ic,jc,mcapa)*valc(ic,jc,ivar)-
     &          capaqfine/(lratiox*lratioy)
       dcapamax = dmax1(dcapamax,dabs(dcapaq))
      
       do 30 ico = 1, lratiox
       ifine = (ic-2)*lratiox + nghost + ico
       do 30 jco = 1, lratioy
         jfine = (jc-2)*lratioy + nghost + jco
         val(ifine,jfine,ivar) = val(ifine,jfine,ivar) +
     &                           dcapaq/aux(ifine,jfine,mcapa)
30     continue

15     continue

10     continue

c      write(6,*)" max discrepancy ", dcapamax

       return
       end
