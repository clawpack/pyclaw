      subroutine rp1(maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  fwave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for the nonlinear elastic equations in 1d,
c     #  variable coefficients
c     #   eps_t - (m/rho(x))_x = 0
c     #   m_t - sigma(eps,x)_x =0
c     # where eps=strain, m=rho*u=momentum
c
c     # aux(i,1) = rho(i)
c     # function sigma(eps,i) gives stress-strain relation in ith cell
c     # function sigmap(eps,i) gives d/(d eps) of sigma
c     #    For linear:   sigma(eps,i) = K_i * eps, and sigmap(eps,i) = K_i
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, fwave contains the waves as jumps in f,
c     #            s the speeds,
c     #
c     #            amdq = A^- Delta q, 
c     #            apdq = A^+ Delta q,
c     #                   the decomposition of the flux difference
c     #                       f(qr(i-1)) - f(ql(i))
c     #                   into leftgoing and rightgoing parts respectively.
c     #
c
c     # Note that the ith Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      dimension auxl(1-mbc:maxm+mbc, 3)
      dimension auxr(1-mbc:maxm+mbc, 3)
      dimension fwave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension apdq(1-mbc:maxm+mbc, meqn)
      dimension amdq(1-mbc:maxm+mbc, meqn)
c
c     local arrays
c     ------------
      dimension delta(2)
c
c
c     # split the jump in q at each interface into waves
c
      do 20 i = 2-mbc, mx+mbc
         rhoi = auxl(i,1)
         rhoim = auxr(i-1,1)
         epsi = ql(i,1)
         epsim = qr(i-1,1)
         urhoi = ql(i,2)
         urhoim = qr(i-1,2)

c        #linearize on each side:

         bulki = sigmap(epsi,i,auxl(i,2),auxl(i,3))
         bulkim = sigmap(epsim,i-1,auxr(i-1,2),auxr(i-1,3))
         ci = dsqrt(bulki / rhoi)
         cim = dsqrt(bulkim / rhoim)
         zi = ci*rhoi
         zim = cim*rhoim

         du = urhoi/rhoi - urhoim/rhoim
         dsig = sigma(epsi,i,auxl(i,2),auxl(i,3)) 
     &          - sigma(epsim,i-1,auxr(i-1,2),auxr(i-1,3))
         b1 = -(zi*du + dsig) / (zim + zi)
         b2 = -(zim*du - dsig) / (zim + zi)
c         a1 = b1 / (-cim)
c         a2 = b2 / ci

c         estarim = epsim + a1
c         estari = epsi - a2

c         ui = urhoi / rhoi
c         uim = urhoim / rhoim

c         ustar = urhoim/rhoim + (estarim - epsim) * cim
c               = urhoi/rhoi - (estari - epsi) * ci
c               = uim - b1
        
c
c        # Compute the waves.
c
         fwave(i,1,1) = b1
         fwave(i,2,1) = b1 * zim
         s(i,1) = -cim
c
         fwave(i,1,2) = b2
         fwave(i,2,2) = b2*(-zi)
         s(i,2) = ci

   20    continue
c
c
c     # compute the leftgoing and rightgoing fluctuations:
c     # Note s(i,1) < 0   and   s(i,2) > 0.
c
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            amdq(i,m) = fwave(i,m,1)
            apdq(i,m) = fwave(i,m,2)
  220       continue
c
      return
      end


c     --------------------------------------------
      double precision function sigma(eps,i,coefA,coefB)
c     --------------------------------------------
      implicit double precision (a-h,o-z)

c     # stress-strain relation in ith cell


c     # nonlinear in both layers:
c     sigma = (coefB + coefA)*(dexp(eps) - 1.d0)

c     # nonlinear layer 0:
c     sigma = coefB*eps + coefA*(dexp(eps) - 1.d0)

c     # nonlinear layer 1:
c     sigma = coefB*(dexp(eps) - 1.d0) + coefA*eps

      if (coefA.gt.0.d0) then
c         # layer A:
c         sigma = coefA*eps 
          sigma = dexp(coefA*eps) - 1.d0
c          sigma = coefA*eps + 0.5d0*eps**2
        else
c         # layer B:
          sigma = coefB*eps 
c         sigma = dexp(coefB*eps) - 1.d0
        endif

      return
      end


c     --------------------------------------------
      double precision function sigmap(eps,i,coefA,coefB)
c     --------------------------------------------
      implicit double precision (a-h,o-z)

c     # derivative of stress-strain relation in ith cell

c     # nonlinear in both layers:
c     sigmap = (coefB + coefA)*dexp(eps)

c     # nonlinear layer 0:
c     sigmap = coefB + coefA*dexp(eps)

c     # nonlinear layer 1:
c     sigmap = coefB*dexp(eps) + coefA

      if (coefA.gt.0.d0) then
c         # layer A:
c         sigmap = coefA
          sigmap = coefA*dexp(coefA*eps)
c         sigmap = coefA + eps
        else
c         # layer B:
          sigmap = coefB
c         sigmap = coefB*dexp(coefB*eps)
        endif


      return
      end
