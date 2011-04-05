module reconstruct
    double precision, allocatable, private :: dq1m(:)
    double precision, allocatable, private :: uu(:,:),dq(:,:)
    double precision, allocatable, private :: uh(:,:,:),gg(:,:),hh(:,:),u(:,:,:)
    double precision, allocatable, private :: evl(:,:,:),evr(:,:,:)
    double precision, private  :: epweno = 1.e-36
    logical :: recon_alloc = .False.

contains

    subroutine alloc_recon_workspace(maxnx,mbc,meqn,mwaves,lim_type,char_decomp)
        integer,intent(in) :: maxnx,mbc,meqn,mwaves,char_decomp,lim_type

        select case(lim_type)
            case(1)
            select case(char_decomp)
                case(1) ! Storage for tvd2_wave()
                    allocate(uu(1-mbc:maxnx+mbc,mwaves))
                case(2) ! Storage for tvd2_char()
                    allocate(dq(1-mbc:maxnx+mbc,meqn))
                    allocate( u(1-mbc:maxnx+mbc,meqn,2))
                    allocate(hh(1-mbc:maxnx+mbc,-1:1))
            end select
            case(2)
            select case(char_decomp)
                case(0)
                    allocate(uu(maxnx+2*mbc,2))
                    allocate( dq1m(maxnx+2*mbc))
                case(2) ! Storage for weno5_char
                    allocate(dq(maxnx+2*mbc,meqn))
                    allocate(uu(maxnx+2*mbc,2))
                    allocate(hh(maxnx+2*mbc,-2:2))
                case(3) ! Storage for weno5_trans
                    allocate(dq(maxnx+2*mbc,meqn))
                    allocate(gg(maxnx+2*mbc,meqn))
                    allocate( u(maxnx+2*mbc,meqn,2))
                    allocate(hh(maxnx+2*mbc,-2:2))
                    allocate(uh(maxnx+2*mbc,meqn,2))
            end select
            recon_alloc = .True.
        end select

    end subroutine alloc_recon_workspace


    subroutine dealloc_recon_workspace(lim_type,char_decomp)
        integer,intent(in) :: lim_type,char_decomp

        select case(lim_type)
            case(1)
            select case(char_decomp)
                case(1) ! Storage for tvd2_wave()
                    deallocate(uu)
                case(2) ! Storage for tvd2_char()
                    deallocate(dq)
                    deallocate( u)
                    deallocate(hh)
            end select
            case(2)
             select case(char_decomp)
                case(0)
                    deallocate(uu)
                    deallocate(dq1m)
                case(2) ! Storage for weno5_char
                    deallocate(dq)
                    deallocate(uu)
                    deallocate(hh)
                case(3) ! Storage for weno5_trans
                    deallocate(dq)
                    deallocate(gg)
                    deallocate( u)
                    deallocate(hh)
                    deallocate(uh)
            end select
            recon_alloc = 0
        end select
    end subroutine dealloc_recon_workspace

    ! ===================================================================
    subroutine weno5(q,ql,qr,meqn,maxnx,mbc)
    ! ===================================================================

        implicit double precision (a-h,o-z)

        double precision, intent(in) :: q(meqn,maxnx+2*mbc)
        double precision, intent(out) :: ql(meqn,maxnx+2*mbc),qr(meqn,maxnx+2*mbc)

        integer, parameter :: mbc=3
        integer :: meqn, mx2

        mx2  = size(q,1); meqn = size(q,2)

        !loop over all equations (all components).  
        !the reconstruction is performed component-wise;
        !no characteristic decomposition is used here

        do m=1,meqn

            forall (i=2:mx2)
                ! compute and store the differences of the cell averages
                dq1m(i)=q(i,m)-q(i-1,m)
            end forall

            ! the reconstruction

            do m1=1,2

                ! m1=1: construct ql
                ! m1=2: construct qr

                im=(-1)**(m1+1)
                ione=im; inone=-im; intwo=-2*im
  
                do i=mbc,mx2-mbc+1
  
                    t1=im*(dq1m(i+intwo)-dq1m(i+inone))
                    t2=im*(dq1m(i+inone)-dq1m(i      ))
                    t3=im*(dq1m(i      )-dq1m(i+ione ))
  
                    tt1=13.*t1**2+3.*(   dq1m(i+intwo)-3.*dq1m(i+inone))**2
                    tt2=13.*t2**2+3.*(   dq1m(i+inone)+   dq1m(i      ))**2
                    tt3=13.*t3**2+3.*(3.*dq1m(i      )-   dq1m(i+ione ))**2
       
                    tt1=(epweno+tt1)**2
                    tt2=(epweno+tt2)**2
                    tt3=(epweno+tt3)**2
                    s1 =tt2*tt3
                    s2 =6.*tt1*tt3
                    s3 =3.*tt1*tt2
                    t0 =1./(s1+s2+s3)
                    s1 =s1*t0
                    s3 =s3*t0
  
                    uu(i,m1) = (s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3. &
                             +(-q(i-2,m)+7.*(q(i-1,m)+q(i,m))-q(i+1,m))/12.

                end do
            end do

           qr(mbc-1:mx2-mbc,  m)=uu(mbc:mx2-mbc+1,1)
           ql(mbc  :mx2-mbc+1,m)=uu(mbc:mx2-mbc+1,2)

        end do

      return
      end subroutine weno5

    ! ===================================================================
    subroutine weno5_char(q,ql,qr,evl,evr)
    ! ===================================================================

        ! This one uses characteristic decomposition
        !  evl, evr are left and right eigenvectors at each interface

        implicit double precision (a-h,o-z)

        double precision, intent(in) :: q(:,:)
        double precision, intent(out) :: ql(:,:),qr(:,:)
        double precision, intent(in) :: evl(:,:,:),evr(:,:,:)

        integer, parameter :: mbc=3
        integer :: meqn, mx2

        mx2  = size(q,1); meqn = size(q,2)

        ! loop over all equations (all components).  
        ! the reconstruction is performed using characteristic decomposition

        forall(m=1:meqn,i=2:mx2)
            ! compute and store the differences of the cell averages
            dq(i,m)=q(i,m)-q(i-1,m)
        end forall

        forall(m=1:meqn,i=3:mx2-1)
            ! Compute the part of the reconstruction that is
            ! stencil-independent
            qr(i-1,m) = (-q(i-2,m)+7.*(q(i-1,m)+q(i,m))-q(i+1,m))/12.
            ql(i,m)   = qr(i-1,m)
        end forall

        do ip=1,meqn

            ! Project the difference of the cell averages to the
            ! 'm'th characteristic field

        
            do m2 = -2,2
               do  i = mbc+1,mx2-2
                  hh(i,m2) = 0.d0
                  do m=1,meqn 
                    hh(i,m2) = hh(i,m2)+ evl(i,ip,m)*dq(i+m2,m)
                  enddo
               enddo
            enddo


            ! the reconstruction

            do m1=1,2

                ! m1=1: construct ql
                ! m1=2: construct qr

                im=(-1)**(m1+1)
                ione=im
                inone=-im
                intwo=-2*im
  
                do i=mbc,mx2-mbc+1
      
                    t1=im*(hh(i,intwo)-hh(i,inone))
                    t2=im*(hh(i,inone)-hh(i,0    ))
                    t3=im*(hh(i,0    )-hh(i,ione ))
      
                    tt1=13.*t1**2+3.*(   hh(i,intwo)-3.*hh(i,inone))**2
                    tt2=13.*t2**2+3.*(   hh(i,inone)+   hh(i,0    ))**2
                    tt3=13.*t3**2+3.*(3.*hh(i,0    )-   hh(i,ione ))**2

                    tt1=(epweno+tt1)**2
                    tt2=(epweno+tt2)**2
                    tt3=(epweno+tt3)**2
                    s1 =tt2*tt3
                    s2 =6.*tt1*tt3
                    s3 =3.*tt1*tt2
                    t0 =1./(s1+s2+s3)
                    s1 =s1*t0
                    s3 =s3*t0
      
                    uu(i,m1) = ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.

                end do !end loop over interfaces
            end do !end loop over which side of interface

                ! Project to the physical space:
            do m = 1,meqn
                do i=mbc,mx2-mbc+1
                    qr(i-1,m) = qr(i-1,m) + evr(i,m,ip)*uu(i,1)
                    ql(i  ,m) = ql(i  ,m) + evr(i,m,ip)*uu(i,2)
                enddo
            enddo
        enddo !end loop over waves

      return
      end subroutine weno5_char

    ! ===================================================================
    subroutine weno5_trans(q,ql,qr,evl,evr)
    ! ===================================================================

        ! Transmission-based WENO reconstruction

        implicit double precision (a-h,o-z)

        double precision, intent(in) :: q(:,:)
        double precision, intent(out) :: ql(:,:),qr(:,:)
        double precision, intent(in) :: evl(:,:,:),evr(:,:,:)

        integer, parameter :: mbc=3
        integer :: meqn, mx2

        mx2  = size(q,1); meqn = size(q,2)


        ! the reconstruction is performed using characteristic decomposition

        do m=1,meqn
            ! compute and store the differences of the cell averages
            forall (i=2:mx2)
                dq(i,m)=q(i,m)-q(i-1,m)
            end forall
        enddo

        ! Find wave strengths at each interface
        ! 'm'th characteristic field
        do mw=1,meqn
            do i = 2,mx2
                gg(i,mw) = 0.d0
                do m=1,meqn
                    gg(i,mw) = gg(i,mw)+ evl(i,mw,m)*dq(i,m)
                enddo
            enddo
        enddo

        do mw=1,meqn
            ! Project the waves to the
            ! 'm'th characteristic field

            do m1 = -2,2
                do  i = mbc+1,mx2-2
                    hh(i,m1) = 0.d0
                    do m=1,meqn 
                        hh(i,m1) = hh(i,m1)+evl(i,mw,m)* &
                                    gg(i+m1,mw)*evr(i+m1,m,mw)
                    enddo
                enddo
            enddo

            ! the reconstruction

            do m1=1,2
                ! m1=1: construct ql
                ! m1=2: construct qr
                im=(-1)**(m1+1)
                ione=im; inone=-im; intwo=-2*im
  
                do i=mbc,mx2-mbc+1
  
                    t1=im*(hh(i,intwo)-hh(i,inone))
                    t2=im*(hh(i,inone)-hh(i,0    ))
                    t3=im*(hh(i,0    )-hh(i,ione ))
  
                    tt1=13.*t1**2+3.*(   hh(i,intwo)-3.*hh(i,inone))**2
                    tt2=13.*t2**2+3.*(   hh(i,inone)+   hh(i,0    ))**2
                    tt3=13.*t3**2+3.*(3.*hh(i,0    )-   hh(i,ione ))**2
       
                    tt1=(epweno+tt1)**2
                    tt2=(epweno+tt2)**2
                    tt3=(epweno+tt3)**2
                    s1 =tt2*tt3
                    s2 =6.*tt1*tt3
                    s3 =3.*tt1*tt2
                    t0 =1./(s1+s2+s3)
                    s1 =s1*t0
                    s3 =s3*t0
  
                    u(i,mw,m1) = ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.

                enddo
            enddo
        enddo

        ! Project to the physical space:

        do m1 =  1,2
            do m =  1, meqn
                do i = mbc,mx2-mbc+1
                    uh(i,m,m1) =( -q(i-2,m) + 7*( q(i-1,m)+q(i,m) ) &
                                         - q(i+1,m) )/12.
                    do mw=1,meqn 
                        uh(i,m,m1) = uh(i,m,m1) + evr(i,m,mw)*u(i,mw,m1)
                    enddo
                enddo
            enddo
        enddo

        qr(mbc-1:mx2-mbc,1:meqn) = uh(mbc:mx2-mbc+1,1:meqn,1)
        ql(mbc:mx2-mbc+1,1:meqn) = uh(mbc:mx2-mbc+1,1:meqn,2)

        return
    end subroutine weno5_trans

    ! ===================================================================
    subroutine weno5_wave(q,ql,qr,wave)
    ! ===================================================================

        !  Fifth order WENO reconstruction, based on waves
        !  which are later interpreted as slopes.

        implicit double precision (a-h,o-z)

        double precision, intent(in) :: q(:,:)
        double precision, intent(out) :: ql(:,:),qr(:,:)
        double precision, intent(in) :: wave(:,:,:)
        double precision u(2)

        integer, parameter :: mbc=3
        integer :: meqn, mx2

        mx2  = size(q,1); meqn = size(q,2); mwaves=size(wave,3)

        ! loop over interfaces (i-1/2)
        do i=2,mx2
            ! Compute the part of the reconstruction that is stencil-independent
            do m=1,meqn
                qr(i-1,m) = (-q(i-2,m)+7.*(q(i-1,m)+q(i,m))-q(i+1,m))/12.
                ql(i,m)   = qr(i-1,m)
            enddo
            ! the reconstruction is performed in terms of waves
            do mw=1,mwaves
                ! loop over which side of x_i-1/2 we're on
                do m1=1,2
                    ! m1=1: construct q^-_{i-1/2}
                    ! m1=2: construct q^+_{i-1/2}
                    im=(-1)**(m1+1)
                    ione=im; inone=-im; intwo=-2*im
  
                    wnorm2 = wave(i      ,1,mw)*wave(i,1,mw)
                    theta1 = wave(i+intwo,1,mw)*wave(i,1,mw)
                    theta2 = wave(i+inone,1,mw)*wave(i,1,mw)
                    theta3 = wave(i+ione ,1,mw)*wave(i,1,mw)
                    do m=2,meqn
                        wnorm2 = wnorm2 + wave(i,      m,mw)*wave(i,m,mw)
                        theta1 = theta1 + wave(i+intwo,m,mw)*wave(i,m,mw)
                        theta2 = theta2 + wave(i+inone,m,mw)*wave(i,m,mw)
                        theta3 = theta3 + wave(i+ione ,m,mw)*wave(i,m,mw)
                    enddo

                    t1=im*(theta1-theta2)
                    t2=im*(theta2-wnorm2)
                    t3=im*(wnorm2-theta3)
  
                    tt1=13.*t1**2+3.*(theta1   -3.*theta2)**2
                    tt2=13.*t2**2+3.*(theta2   +   wnorm2)**2
                    tt3=13.*t3**2+3.*(3.*wnorm2-   theta3)**2
       
                    tt1=(epweno+tt1)**2
                    tt2=(epweno+tt2)**2
                    tt3=(epweno+tt3)**2
                    s1 =tt2*tt3
                    s2 =6.*tt1*tt3
                    s3 =3.*tt1*tt2
                    t0 =1./(s1+s2+s3)
                    s1 =s1*t0
                    s3 =s3*t0
  
                    if(wnorm2.gt.1.e-14) then
                        u(m1) = ( s1*(t2-t1) + (0.5*s3-0.25)*(t3-t2) ) /3.
                        wnorm2=1.d0/wnorm2
                    else
                        u(m1) = 0.d0
                        wnorm2=0.d0
                    endif
                enddo !end loop over which side of interface
                do m=1,meqn
                    qr(i-1,m) = qr(i-1,m) +  u(1)*wave(i,m,mw)*wnorm2
                    ql(i  ,m) = ql(i  ,m) +  u(2)*wave(i,m,mw)*wnorm2
                enddo
            enddo !loop over waves
        enddo !loop over interfaces

    end subroutine weno5_wave

    ! ===================================================================
    subroutine weno5_fwave(q,ql,qr,fwave,s)
    ! ===================================================================
    !
    !  Fifth order WENO reconstruction, based on f-waves
    !  that are interpreted as slopes.
!

      implicit double precision (a-h,o-z)

      double precision, intent(in) :: q(:,:)
      double precision, intent(inout) :: fwave(:,:,:), s(:,:)
      double precision, intent(out) :: ql(:,:), qr(:,:)
      double precision  u(20,2)

      integer, parameter :: mbc=3
      integer :: meqn, mx2

      mx2= size(q,1); meqn = size(q,2); mwaves=size(fwave,3)

      ! convert fwaves to waves by dividing by the sound speed
      ! We do this in place to save memory
      ! and get away with it because the waves are never used again
      forall(i=1:mx2,mw=1:mwaves,m=1:meqn)
          fwave(i,m,mw)=fwave(i,m,mw)/s(i,mw)
      end forall

      ! loop over interfaces (i-1/2)
      do i=2,mx2
        ! Compute the part of the reconstruction that is
        !  stencil-independent
        do m=1,meqn
          qr(i-1,m) = q(i-1,m)
          ql(i  ,m) = q(i,m)
        enddo
        ! the reconstruction is performed in terms of waves
        do mw=1,mwaves
         ! loop over which side of x_i-1/2 we're on
          do m1=1,2
            ! m1=1: construct q^-_{i-1/2}
            ! m1=2: construct q^+_{i-1/2}

            im=(-1)**(m1+1)
            ione=im; inone=-im; intwo=-2*im
  
            ! compute projections of waves in each family
            ! onto the corresponding wave at the current interface
            wnorm2 = fwave(i      ,1,mw)*fwave(i,1,mw)
            theta1 = fwave(i+intwo,1,mw)*fwave(i,1,mw)
            theta2 = fwave(i+inone,1,mw)*fwave(i,1,mw)
            theta3 = fwave(i+ione ,1,mw)*fwave(i,1,mw)
            do m=2,meqn
              wnorm2 = wnorm2 + fwave(i,      m,mw)*fwave(i,m,mw)
              theta1 = theta1 + fwave(i+intwo,m,mw)*fwave(i,m,mw)
              theta2 = theta2 + fwave(i+inone,m,mw)*fwave(i,m,mw)
              theta3 = theta3 + fwave(i+ione ,m,mw)*fwave(i,m,mw)
            enddo
!
             t1=im*(theta1-theta2)
             t2=im*(theta2-wnorm2)
             t3=im*(wnorm2-theta3)
  
             tt1=13.*t1**2+3.*(theta1   -3.*theta2)**2
             tt2=13.*t2**2+3.*(theta2   +   wnorm2)**2
             tt3=13.*t3**2+3.*(3.*wnorm2-   theta3)**2
       
             tt1=(epweno+tt1)**2
             tt2=(epweno+tt2)**2
             tt3=(epweno+tt3)**2
             s1 =tt2*tt3
             s2 =6.*tt1*tt3
             s3 =3.*tt1*tt2
             t0 =1./(s1+s2+s3)
             s1 =s1*t0
             s3 =s3*t0
  
           if(wnorm2.gt.1.e-14) then
             u(mw,m1) = ( (s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3. &
                       + im*(theta2+6.d0*wnorm2-theta3)/12.d0)
             wnorm2=1.d0/wnorm2
           else
             u(mw,m1) = 0.d0
             wnorm2=0.d0
           endif
          enddo !end loop over which side of interface
          do m=1,meqn
            qr(i-1,m) = qr(i-1,m) +  u(mw,1)*fwave(i,m,mw)*wnorm2
            ql(i  ,m) = ql(i  ,m) +  u(mw,2)*fwave(i,m,mw)*wnorm2
          enddo
        enddo !loop over fwaves
      enddo !loop over interfaces

    end subroutine weno5_fwave

    ! ===================================================================
    subroutine tvd2(q,ql,qr,mthlim)
    ! ===================================================================
    ! Second order TVD reconstruction

        implicit double precision (a-h,o-z)

        double precision, intent(in) :: q(:,:)
        integer, intent(in) :: mthlim(:)
        double precision, intent(out) :: ql(:,:),qr(:,:)
        integer, parameter :: mbc=2
        integer :: meqn, mx2

        mx2  = size(q,1); meqn = size(q,2); 

        ! loop over all equations (all components).  
        ! the reconstruction is performed component-wise

        do m=1,meqn

            ! compute and store the differences of the cell averages

            do i=mbc+1,mx2-mbc
                dqm=dqp
                dqp=q(i+1,m)-q(i,m)
                r=dqp/dqm

                select case(mthlim(m))
                case(1)
                    ! minmod
                    qlimitr = dmax1(0.d0, dmin1(1.d0, r))
                case(2)
                    ! superbee
                    qlimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
                case(3)
                    ! van Leer
                    qlimitr = (r + dabs(r)) / (1.d0 + dabs(r))
                case(4)
                    ! monotonized centered
                    c = (1.d0 + r)/2.d0
                    qlimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
                case(5)
                    ! Cada & Torrilhon simple
                    beta=2.d0
                    xgamma=2.d0
                    alpha=1.d0/3.d0
                    pp=(2.d0+r)/3.d0
                    amax = dmax1(-alpha*r,0.d0,dmin1(beta*r,pp,xgamma))
                    qlimitr = dmax1(0.d0, dmin1(pp,amax))
                end select

           qr(i,m) = q(i,m) + 0.5d0*qlimitr*dqm
           ql(i,m) = q(i,m) - 0.5d0*qlimitr*dqm

         enddo
      enddo

      return
      end subroutine tvd2


    ! ===================================================================
    subroutine tvd2_char(q,ql,qr,mthlim,evl,evr)
    ! ===================================================================

        ! Second order TVD reconstruction for WENOCLAW
        ! This one uses characteristic decomposition

        !  evl, evr are left and right eigenvectors at each interface
        implicit double precision (a-h,o-z)

        double precision, intent(in) :: q(:,:)
        integer, intent(in) :: mthlim(:)
        double precision, intent(out) :: ql(:,:),qr(:,:)
        double precision, intent(in) :: evl(:,:,:),evr(:,:,:)
        integer, parameter :: mbc=2
        integer :: meqn, mx2

        mx2  = size(q,1); meqn = size(q,2); 

        ! loop over all equations (all components).  
        ! the reconstruction is performed using characteristic decomposition

        ! compute and store the differences of the cell averages
        forall(m=1:meqn,i=2:mx2)
            dq(i,m)=q(i,m)-q(i-1,m)
        end forall

        do m=1,meqn

            ! Project the difference of the cell averages to the
            ! 'm'th characteristic field
            do m1 = -1,1
                do  i = mbc+1,mx2-1
                    hh(i,m1) = 0.d0
                    do mm=1,meqn
                        hh(i,m1) = hh(i,m1)+ evl(i,m,mm)*dq(i+m1,mm)
                    enddo
                enddo
            enddo


            ! the reconstruction
            do m1=1,2
                im=(-1)**(m1+1)
                ! m1=1: construct qr_i-1
                ! m1=2: construct ql_i

                do i=mbc+1,mx2-1
                    ! dqp=hh(i,m1-1)
                    ! dqm=hh(i,m1-2)
                    if (dabs(hh(i,m1-2)).gt.1.e-14) then
                        r=hh(i,m1-1)/hh(i,m1-2)
                    else
                        r=0.d0
                    endif
               
                    select case(mthlim(m))
                    case(1)
                        ! minmod
                        slimitr = dmax1(0.d0, dmin1(1.d0, r))
                    case(2)
                        ! superbee
                        slimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
                    case(3)
                        ! van Leer
                        slimitr = (r + dabs(r)) / (1.d0 + dabs(r))
                    case(4)
                        ! monotonized centered
                        c = (1.d0 + r)/2.d0
                        slimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
                    case(5)
                        ! Cada & Torrilhon simple
                        beta=2.d0
                        xgamma=2.d0
                        alpha=1.d0/3.d0
                        pp=(2.d0+r)/3.d0
                        amax = dmax1(-alpha*r,0.d0,dmin1(beta*r,pp,xgamma))
                        slimitr = dmax1(0.d0, dmin1(pp,amax))
                    end select
    
                     u(i,m,m1) = im*0.5d0*slimitr*hh(i,m1-2)

                enddo
            enddo
        enddo

        ! Project to the physical space:
        do m =  1, meqn
            do i = mbc+1,mx2-1
                qr(i-1,m)=q(i-1,m)
                ql(i  ,m)=q(i  ,m)
                do mm=1,meqn 
                    qr(i-1,m) = qr(i-1,m) + evr(i,m,mm)*u(i,mm,1)
                    ql(i  ,m) = ql(i  ,m) + evr(i,m,mm)*u(i,mm,2)
                enddo
            enddo
        enddo
    end subroutine tvd2_char

    ! ===================================================================
    subroutine tvd2_wave(q,ql,qr,wave,s,mthlim)
    ! ===================================================================
        ! Second order TVD reconstruction for WENOCLAW
        ! This one uses projected waves

        implicit double precision (a-h,o-z)

        double precision, intent(in) :: q(:,:)
        integer, intent(in) :: mthlim(:)
        double precision, intent(out) :: ql(:,:),qr(:,:)
        double precision, intent(in) :: wave(:,:,:), s(:,:)
        integer, parameter :: mbc=2
        integer :: meqn, mx2

        mx2  = size(q,1); meqn = size(q,2); mwaves=size(wave,3)

        forall(i=2:mx2,m=1:meqn)
            qr(i-1,m) = q(i-1,m)
            ql(i  ,m) = q(i  ,m)
        end forall

        ! loop over all equations (all components).  
        ! the reconstruction is performed using characteristic decomposition

        do mw=1,mwaves
            dotr = 0.d0
            do i=mbc,mx2-mbc
                wnorm2=0.d0
                dotl=dotr
                dotr=0.d0
                do m=1,meqn
                    wnorm2 = wnorm2 + wave(i,m,mw)**2
                    dotr = dotr + wave(i,m,mw)*wave(i+1,m,mw)
                enddo
                if (i.eq.0) cycle
                if (wnorm2.eq.0.d0) cycle
                
                if (s(i,mw).gt.0.d0) then
                    r = dotl / wnorm2
                else
                    r = dotr / wnorm2
                endif

                select case(mthlim(mw))
                    case(1)
                        ! minmod
                        wlimitr = dmax1(0.d0, dmin1(1.d0, r))
                    case(2)
                        ! superbee
                        wlimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))
                    case(3)
                        ! van Leer
                        wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))
                    case(4)
                        ! monotonized centered
                        c = (1.d0 + r)/2.d0
                        wlimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
                    case(5)
                        ! Cada & Torrilhon simple
                        beta=2.d0
                        xgamma=2.d0
                        alpha=1.d0/3.d0
                        pp=(2.d0+r)/3.d0
                        amax = dmax1(-alpha*r,0.d0,dmin1(beta*r,pp,xgamma))
                        wlimitr = dmax1(0.d0, dmin1(pp,amax))
                end select

                uu(i,mw) = 0.5d0*wlimitr

                do m =  1, meqn
                    qr(i-1,m) = qr(i-1,m) + wave(i,m,mw)*uu(i,mw)
                    ql(i  ,m) = ql(i  ,m) - wave(i,m,mw)*uu(i,mw)
                enddo ! end loop over equations

            enddo
        enddo !end loop over waves

      return
      end subroutine tvd2_wave

end module reconstruct

