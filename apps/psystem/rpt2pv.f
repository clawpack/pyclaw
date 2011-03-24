c
c
c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  imp,asdq,bmasdq,bpasdq)
c     =====================================================
c
c     # Riemann solver in the transverse direction.
c     # This is a dummy routine that returns zeros and is only intended
c     # to illustrate the format of this routine.  See various example
c     # directories for better examples.

c     # This dummy routine can be used if transverse solves are not being
c     # used, i.e. if method(3)=0.
c
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c
      implicit double precision (a-h,o-z)
      dimension     ql(meqn,1-mbc:maxm+mbc)
      dimension     qr(meqn,1-mbc:maxm+mbc)
      dimension   asdq(meqn,1-mbc:maxm+mbc)
      dimension bmasdq(meqn,1-mbc:maxm+mbc)
      dimension bpasdq(meqn,1-mbc:maxm+mbc)
      dimension   aux1(3,1-mbc:maxm+mbc)
      dimension   aux2(3,1-mbc:maxm+mbc)
      dimension   aux3(3,1-mbc:maxm+mbc)
      dimension s (2,1-mbc:maxm+mbc)

      do 10 i = 2-mbc, mx+mbc
         if (imp.eq.1) then !left going fluctuation
	     ix = i-1
	   else !right going fluctuation
	     ix = i
	   endif

c material properties 
         pjm=aux1(1,ix) !density at (ix, j-1) where ix=i or i-1
         pj=aux2(1,ix)  
         pjp=aux3(1,ix)
         Ejm=aux1(2,ix)
         Ej=aux2(2,ix)
         Eip=aux3(2,ix)

c linearity of the material
         linearity_matjm=aux1(3,ix)
         linearity_matj=aux2(3,ix)
         linearity_matjp=aux3(3,ix)

c eps (strain) at different rows
         epsjm=aux1(4,ix)
         epsj=aux2(4,ix)
         epsjp=aux3(4,ix)

c sigmap
         sigmapjm=sigmap(epsjm,Ejm,linearity_matjm)
         sigmapj=sigmap(epsj,Ej,linearity_matj)
         sigmapjp=sigmap(epsjp,Ejp,linearity_matjp)

c computation of components of eigenvectors 
         r11=1/dsqrt(sigmapjm*pjm) 
         r13=-1/dsqrt(sigmapjp*pjp) 

c shock speeds (eigenvalues of A and B)
         s(1,i)=-dsqrt(sigmapjm/pjm)  !lambda_1, lambda_2=0
         s(2,i)=dsqrt(sigmapjp/pjp)     !lambda_3

         if(ixy.eq.1) then !x direction, use eigenvectors of matrix B
            !for the decomposition of the fluctuations
c computation of coefficients of the decomposition of the fluctuations. 
c I called them gamma since beta is used in the decomposition
c of the flux difference (notation used in paper for the f-wave)
            gamma1=(asdq(1,i)-r13*asdq(3,i))/(r11-r13)
            gamma3=(-asdq(1,i)+r11*asdq(3,i))/(r11-r13)
c computation of the fluctuations 
            !down going part of the normal fluctuation
            bmasdq(1,i)=s(1,i)*gamma1*r11
            bmasdq(2,i)=s(1,i)*gamma1*0
            bmasdq(3,i)=s(1,i)*gamma1*1
           !up going part of the normal fluctuation
            bpasdq(1,i)=s(2,i)*gamma3*r13
            bpasdq(2,i)=s(2,i)*gamma3*0
            bpasdq(3,i)=s(2,i)*gamma3*1

         else !y direction, use eigenvectors of matrix A
            !for the decomposition of the fluctuations
            gamma1=(asdq(1,i)-r13*asdq(2,i))/(r11-r13)
            gamma3=(-asdq(1,i)+r11*asdq(2,i))/(r11-r13)
c computation of fluctuations
            !"down" going part
            bmasdq(1,i)=s(1,i)*gamma1*r11
            bmasdq(2,i)=s(1,i)*gamma1*1
            bmasdq(3,i)=s(1,i)*gamma1*0
            !"up" going part
            bpasdq(1,i)=s(1,i)*gamma3*r13
            bpasdq(2,i)=s(1,i)*gamma3*1
            bpasdq(3,i)=s(1,i)*gamma3*0
         endif


 10      continue

      return
      end
