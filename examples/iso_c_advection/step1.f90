module step1m
implicit none

interface
   subroutine rp(mx, mbc, ql, qr, s, wave, amdq, apdq, context)
     use iso_c_binding
     integer(c_int), intent(in), value :: mx, mbc
     real(c_double), intent(in) :: &
          ql(1, 1-mbc:mx+mbc), &
          qr(1, 1-mbc:mx+mbc)
     type(c_ptr), intent(in), value :: context

     real(c_double), intent(out) :: &
          s(1, 1-mbc:mx+mbc), &
          wave(1, 1, 1-mbc:mx+mbc), &
          apdq(1, 1-mbc:mx+mbc), &
          amdq(1, 1-mbc:mx+mbc)
   end subroutine rp
end interface

contains

  ! ===================================================================
  subroutine step1(num_eqn,num_waves,num_ghost,num_aux,mx,q,aux,dx, &
       dt,method,mthlim,cfl,f,wave,s,amdq,apdq,dtdx,use_fwave,c_rp1,context) &
       bind(C, name='step1_c')
    ! ===================================================================

    !     # Take one time step, updating q.

    !     method(1) = 1   ==>  Godunov method
    !     method(1) = 2   ==>  Slope limiter method
    !     mthlim(p)  controls what limiter is used in the pth family


    !     amdq, apdq, wave, s, and f are used locally:

    !     amdq(1-num_ghost:mx+num_ghost, num_eqn) = left-going flux-differences
    !     apdq(1-num_ghost:mx+num_ghost, num_eqn) = right-going flux-differences
    !        e.g. amdq(m,i) = m'th component of A^- \Delta q from i'th Riemann
    !                         problem (between cells i-1 and i).

    !     wave(1-num_ghost:mx+num_ghost, num_eqn, num_waves) = waves from solution of
    !                                           Riemann problems,
    !            wave(m,mw,i) = mth component of jump in q across
    !                           wave in family mw in Riemann problem between
    !                           states i-1 and i.

    !     s(1-num_ghost:mx+num_ghost, num_waves) = wave speeds,
    !            s(m,iw) = speed of wave in family mw in Riemann problem between
    !                      states i-1 and i.

    !     f(1-num_ghost:mx+num_ghost, num_eqn) = correction fluxes for second order method
    !            f(m,i) = mth component of flux at left edge of ith cell
    !     --------------------------------------------------------------------

  use iso_c_binding

    logical(c_bool), intent(in) :: use_fwave
    integer(c_int), intent(in) :: &
         num_eqn, num_waves, num_ghost, num_aux, mx, &
         method(7), mthlim(num_waves)

    real(c_double), intent(in) :: &
         aux(num_aux,1-num_ghost:mx+num_ghost)

    real(c_double), intent(inout) :: &
         q(num_eqn,1-num_ghost:mx+num_ghost)

    type(c_funptr), intent(in)     :: c_rp1
    type(c_ptr), intent(in)        :: context

    real(c_double), intent(out) :: &
         cfl, &
         s(num_waves,1-num_ghost:mx+num_ghost), &
         wave(num_eqn, num_waves,1-num_ghost:mx+num_ghost), &
         amdq(num_eqn,1-num_ghost:mx+num_ghost), &
         apdq(num_eqn,1-num_ghost:mx+num_ghost), &
         dtdx(1-num_ghost:mx+num_ghost), &
         f(num_eqn,1-num_ghost:mx+num_ghost), &
         dx, &
         dt

    logical :: limit
    procedure(rp), pointer :: rp1
    real(c_double) :: dtdxave
    integer(c_int) :: i, m, mw, index_capa

    call c_f_procpointer(c_rp1, rp1)

!     # check if any limiters are used:
    limit = .false.
    do mw=1,num_waves
       if (mthlim(mw) > 0) limit = .TRUE.
    end do

    index_capa = method(6)
    do i=1-num_ghost,mx+num_ghost
        if (index_capa > 0) then
            if (aux(index_capa,i) <= 0.d0) then
                write(6,*) 'Error -- capa must be positive'
                stop
            endif
            dtdx(i) = dt / (dx*aux(index_capa,i))
        else
            dtdx(i) = dt/dx
        endif
    end do

!     # solve Riemann problem at each interface
!     -----------------------------------------

    call rp1(mx,num_ghost,q,q,s,wave,amdq,apdq,context)

!     # Modify q for Godunov update:
!     # Note this may not correspond to a conservative flux-differencing
!     # for equations not in conservation form.  It is conservative if
!     # amdq + apdq = f(q(i)) - f(q(i-1)).

    forall(i=1:mx+1, m=1:num_eqn)
       q(m,i) = q(m,i) - dtdx(i)*apdq(m,i)
       q(m,i-1) = q(m,i-1) - dtdx(i-1)*amdq(m,i)
    end forall


!     # compute maximum wave speed:
    cfl = 0.d0
    do mw=1,num_waves
        do i=1,mx+1
        !          # if s>0 use dtdx(i) to compute CFL,
        !          # if s<0 use dtdx(i-1) to compute CFL:
            cfl = dmax1(cfl, dtdx(i)*s(mw,i), -dtdx(i-1)*s(mw,i))
        end do
    end do

    if (method(2) == 1) then
       return
    endif
!     # compute correction fluxes for second order q_{xx} terms:
!     ----------------------------------------------------------

    forall(i=1-num_ghost:mx+num_ghost, m=1:num_eqn)
       f(m,i) = 0.d0
    end forall

!      # apply limiter to waves:
    if (limit) call limiter(mx,num_eqn,num_waves,num_ghost,mx, &
         wave,s,mthlim)

    if (use_fwave.eqv. .FALSE. ) then
       do i=1,mx+1
          do m=1,num_eqn
             do mw=1,num_waves
                dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
                f(m,i) = f(m,i) + 0.5d0 * dabs(s(mw,i)) &
                     * (1.d0 - dabs(s(mw,i))*dtdxave) * wave(m,mw,i)
             end do
          end do
       end do
    else
        do i=1,mx+1
           do m=1,num_eqn
              do mw=1,num_waves
                 dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
                 f(m,i) = f(m,i) + 0.5d0 * dsign(1.d0,s(mw,i)) &
                      * (1.d0 - dabs(s(mw,i))*dtdxave) * wave(m,mw,i)
              end do
           end do
        end do
     endif



     ! update q by differencing correction fluxes
     ! ============================================

     ! (Note:  Godunov update has already been performed above)

    forall(i=1:mx+1, m=1:num_eqn)
       q(m,i) = q(m,i) - dtdx(i) * (f(m,i+1) - f(m,i))
    end forall

  end subroutine step1
end module step1m
