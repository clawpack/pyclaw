module reconstruct
    double precision, allocatable  :: dq1m(:)
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
                    allocate(uu(mwaves,1-mbc:maxnx+mbc))
                case(2) ! Storage for tvd2_char()
                    ! Do the array bounds here cause a bug?
                    allocate(dq(meqn,1-mbc:maxnx+mbc))
                    allocate( u(meqn,2,1-mbc:maxnx+mbc))
                    allocate(hh(-1:1,1-mbc:maxnx+mbc))
            end select
            case(2)
            select case(char_decomp)
                case(0)
                    allocate(uu(2,maxnx+2*mbc))
                    allocate( dq1m(maxnx+2*mbc))
                case(2) ! Storage for weno5_char
                    allocate(dq(meqn,maxnx+2*mbc))
                    allocate(uu(2,maxnx+2*mbc))
                    allocate(hh(-2:2,maxnx+2*mbc))
                case(3) ! Storage for weno5_trans
                    allocate(dq(meqn,maxnx+2*mbc))
                    allocate(gg(meqn,maxnx+2*mbc))
                    allocate( u(meqn,2,maxnx+2*mbc))
                    allocate(hh(-2:2,maxnx+2*mbc))
                    allocate(uh(meqn,2,maxnx+2*mbc))
            end select
        end select
        recon_alloc = .True.

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
            recon_alloc = .False.
        end select
    end subroutine dealloc_recon_workspace

    ! ===================================================================
    subroutine weno5(q,ql,qr,meqn,maxnx,mbc)
    ! ===================================================================

        implicit double precision (a-h,o-z)

        double precision, intent(in) :: q(meqn,maxnx+2*mbc)
        double precision, intent(out) :: ql(meqn,maxnx+2*mbc),qr(meqn,maxnx+2*mbc)

        integer :: meqn, mx2

        mx2  = size(q,2); meqn = size(q,1)

        !loop over all equations (all components).  
        !the reconstruction is performed component-wise;
        !no characteristic decomposition is used here

        do m=1,meqn

            forall (i=2:mx2)
                ! compute and store the differences of the cell averages
                dq1m(i)=q(m,i)-q(m,i-1)
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
  
                    uu(m1,i) = (s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3. &
                             +(-q(m,i-2)+7.*(q(m,i-1)+q(m,i))-q(m,i+1))/12.

                end do
            end do

           qr(m,mbc-1:mx2-mbc  )=uu(1,mbc:mx2-mbc+1)
           ql(m,mbc  :mx2-mbc+1)=uu(2,mbc:mx2-mbc+1)

        end do

      return
      end subroutine weno5


end module reconstruct

