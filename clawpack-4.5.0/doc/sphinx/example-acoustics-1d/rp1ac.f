
c
c
c     =====================================================
      subroutine rp1(maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  wave,s,amdq,apdq)
c     =====================================================
c
c     # begin_html  
c       [www.clawpack.org CLAWPACK] Riemann solver
c       for constant coefficient acoustics in 1 space dimension.
c          \[   q_t + A q_x = 0.  \]
c       where
c          \[   q(x,t) = \vector{ p(x,t)\\ u(x,t)}           \]
c       and the coefficient matrix is          
c          \[   A = \begin{matrix}
c                     0      & K\\
c                     1/\rho & 0 
c                     \end{matrix}.
c          \]
c
c       The parameters $\rho = $ density and $K =$ bulk modulus
c       are set in [code: setprob.data] and the sound speed $c$ 
c       and impedance $Z$ are determined from these in [code: setprob.f].
c
c       <b>On input:</b>
c         <ul>
c         <li> ql contains the state vector at the left edge of each cell,
c         <li> qr contains the state vector at the right edge of each cell,
c         <li> auxl, auxr are not used in this Riemann solver.
c       </ul>
c
c
c       Note that the i'th Riemann problem has left state qr(i-1,:)
c                                          and right state ql(i,:).
c       From the basic clawpack routine step1, rp1 is called with ql = qr = q.

c       <b>On output:</b>
c         <ul>
c         <li> wave contains the waves, 
c         <li> s the speeds, 
c         <li> amdq the  left-going flux difference  $\A^-\Delta Q$, 
c         <li> apdq the right-going flux difference  $\A^+\Delta Q$ 
c       </ul>
c
c
c       For additional documentation on Riemann solvers rp1, see
c       [claw/doc/rp1]
c
c       For details on solution of the Riemann problem for acoustics,
c       see Chapter 3 of [www.clawpack.org/book.html FVMHP].
c
c
c     # end_html
c
c
c
      implicit double precision (a-h,o-z)
c
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
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
c     # density, bulk modulus, and sound speed, and impedence of medium:
c     # (should be set in setprob.f)
      common /cparam/ rho,bulk,cc,zz   
c
c
c     # begin_html
c       Split the jump in $Q$ at each interface into waves.
c       First find $\alpha^1$ and $\alpha^2$, the coefficients of the 
c       2 eigenvectors:
c         \[ \delta = \alpha^1 \vector{ -Z \\ 1} + 
c                     \alpha^2 \vector{  Z \\ 1}      \]
c       
c       Note that the eigendecomposition of $A$ is 
c       $A = R \Lambda R^{-1}$, with
c          \[   R = \begin{matrix}
c                        -Z & Z \\
c                         1 & 1
c                        \end{matrix}, \quad
c               \Lambda = \begin{matrix}
c                        -c & 0 \\
c                         0 & c
c                        \end{matrix}, \quad
c               R^{-1} = \frac{1}{2Z}
c                        \begin{matrix}
c                        -1 & Z \\
c                         1 & Z 
c                        \end{matrix}, \quad
c          \]
c     # end_html
c
      do 20 i = 2-mbc, mx+mbc
         delta(1) = ql(i,1) - qr(i-1,1)
         delta(2) = ql(i,2) - qr(i-1,2)
         a1 = (-delta(1) + zz*delta(2)) / (2.d0*zz)
         a2 =  (delta(1) + zz*delta(2)) / (2.d0*zz)
c
c        # Compute the waves.
c
         wave(i,1,1) = -a1*zz
         wave(i,2,1) = a1
         s(i,1) = -cc
c
         wave(i,1,2) = a2*zz
         wave(i,2,2) = a2
         s(i,2) = cc
c
   20    continue
c
c
c     # begin_html
c       Compute the leftgoing and rightgoing fluctuations:
c
c       For this problem we always have $s^1 =-c \lt 0$ and $s^2 = c\gt 0$, so
c       \[
c          \A^-\Delta Q = s^1\W^1, \quad 
c          \A^+\Delta Q = s^2\W^2.
c       \]
c     # end_html
c
      do 220 m=1,meqn
         do 220 i = 2-mbc, mx+mbc
            amdq(i,m) = s(i,1)*wave(i,m,1)
            apdq(i,m) = s(i,2)*wave(i,m,2)
  220       continue
c
      return
      end
