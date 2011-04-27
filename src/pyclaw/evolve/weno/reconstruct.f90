module reconstruct
contains


  subroutine smoothness_k3(f, n, sigma)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: sigma(n,0:2)
    integer :: i

    do i=3, n-3
       sigma(i,0) = +3.33333333333333 * f(i+0) * f(i+0) &
            -10.3333333333333 * f(i+0) * f(i+1) &
            +3.66666666666667 * f(i+0) * f(i+2) &
            +8.33333333333333 * f(i+1) * f(i+1) &
            -6.33333333333333 * f(i+1) * f(i+2) &
            +1.33333333333333 * f(i+2) * f(i+2)
       sigma(i,1) = +1.33333333333333 * f(i-1) * f(i-1) &
            -4.33333333333333 * f(i-1) * f(i+0) &
            +1.66666666666667 * f(i-1) * f(i+1) &
            +4.33333333333333 * f(i+0) * f(i+0) &
            -4.33333333333333 * f(i+0) * f(i+1) &
            +1.33333333333333 * f(i+1) * f(i+1)
       sigma(i,2) = +1.33333333333333 * f(i-2) * f(i-2) &
            -6.33333333333333 * f(i-2) * f(i-1) &
            +3.66666666666667 * f(i-2) * f(i+0) &
            +8.33333333333333 * f(i-1) * f(i-1) &
            -10.3333333333333 * f(i-1) * f(i+0) &
            +3.33333333333333 * f(i+0) * f(i+0)
    end do
  end subroutine smoothness_k3

  subroutine weights_left_k3(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:2)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:2,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2
    do i=3, n-3
       accumulator = 0.0
       omega0 = +0.1 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.6 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.3 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
    end do
  end subroutine weights_left_k3
  subroutine reconstruct_left_k3(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:2,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2
    real(kind=8) :: fr0, fr1, fr2
    do i=3, n-3
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       fr0 = +1.83333333333333 * f(i+0) &
            -1.16666666666667 * f(i+1) &
            +0.333333333333333 * f(i+2)
       fr1 = +0.333333333333333 * f(i-1) &
            +0.833333333333333 * f(i+0) &
            -0.166666666666667 * f(i+1)
       fr2 = -0.166666666666667 * f(i-2) &
            +0.833333333333333 * f(i-1) &
            +0.333333333333333 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2
    end do
  end subroutine reconstruct_left_k3
  subroutine weights_right_k3(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:2)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:2,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2
    do i=3, n-3
       accumulator = 0.0
       omega0 = +0.3 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.6 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.1 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
    end do
  end subroutine weights_right_k3
  subroutine reconstruct_right_k3(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:2,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2
    real(kind=8) :: fr0, fr1, fr2
    do i=3, n-3
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       fr0 = +0.333333333333333 * f(i+0) &
            +0.833333333333333 * f(i+1) &
            -0.166666666666667 * f(i+2)
       fr1 = -0.166666666666667 * f(i-1) &
            +0.833333333333333 * f(i+0) &
            +0.333333333333333 * f(i+1)
       fr2 = +0.333333333333333 * f(i-2) &
            -1.16666666666667 * f(i-1) &
            +1.83333333333333 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2
    end do
  end subroutine reconstruct_right_k3
  subroutine smoothness_k4(f, n, sigma)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: sigma(n,0:3)
    integer :: i

    do i=4, n-4
       sigma(i,0) = +8.77916666666667 * f(i+0) * f(i+0) &
            -39.175 * f(i+0) * f(i+1) &
            +29.3416666666667 * f(i+0) * f(i+2) &
            -7.725 * f(i+0) * f(i+3) &
            +45.8458333333333 * f(i+1) * f(i+1) &
            -71.8583333333333 * f(i+1) * f(i+2) &
            +19.3416666666667 * f(i+1) * f(i+3) &
            +29.3458333333333 * f(i+2) * f(i+2) &
            -16.175 * f(i+2) * f(i+3) &
            +2.27916666666667 * f(i+3) * f(i+3)
       sigma(i,1) = +2.27916666666667 * f(i-1) * f(i-1) &
            -10.5083333333333 * f(i-1) * f(i+0) &
            +8.00833333333333 * f(i-1) * f(i+1) &
            -2.05833333333333 * f(i-1) * f(i+2) &
            +14.3458333333333 * f(i+0) * f(i+0) &
            -24.8583333333333 * f(i+0) * f(i+1) &
            +6.675 * f(i+0) * f(i+2) &
            +11.8458333333333 * f(i+1) * f(i+1) &
            -6.84166666666667 * f(i+1) * f(i+2) &
            +1.1125 * f(i+2) * f(i+2)
       sigma(i,2) = +1.1125 * f(i-2) * f(i-2) &
            -6.84166666666667 * f(i-2) * f(i-1) &
            +6.675 * f(i-2) * f(i+0) &
            -2.05833333333333 * f(i-2) * f(i+1) &
            +11.8458333333333 * f(i-1) * f(i-1) &
            -24.8583333333333 * f(i-1) * f(i+0) &
            +8.00833333333333 * f(i-1) * f(i+1) &
            +14.3458333333333 * f(i+0) * f(i+0) &
            -10.5083333333333 * f(i+0) * f(i+1) &
            +2.27916666666667 * f(i+1) * f(i+1)
       sigma(i,3) = +2.27916666666667 * f(i-3) * f(i-3) &
            -16.175 * f(i-3) * f(i-2) &
            +19.3416666666667 * f(i-3) * f(i-1) &
            -7.725 * f(i-3) * f(i+0) &
            +29.3458333333333 * f(i-2) * f(i-2) &
            -71.8583333333333 * f(i-2) * f(i-1) &
            +29.3416666666667 * f(i-2) * f(i+0) &
            +45.8458333333333 * f(i-1) * f(i-1) &
            -39.175 * f(i-1) * f(i+0) &
            +8.77916666666667 * f(i+0) * f(i+0)
    end do
  end subroutine smoothness_k4

  subroutine weights_left_k4(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:3)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:3,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3
    do i=4, n-4
       accumulator = 0.0
       omega0 = +0.0285714285714286 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.342857142857143 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.514285714285714 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.114285714285714 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
    end do
  end subroutine weights_left_k4
  subroutine reconstruct_left_k4(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:3,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3
    real(kind=8) :: fr0, fr1, fr2, fr3
    do i=4, n-4
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       fr0 = +2.08333333333333 * f(i+0) &
            -1.91666666666667 * f(i+1) &
            +1.08333333333333 * f(i+2) &
            -0.25 * f(i+3)
       fr1 = +0.25 * f(i-1) &
            +1.08333333333333 * f(i+0) &
            -0.416666666666667 * f(i+1) &
            +0.0833333333333333 * f(i+2)
       fr2 = -0.0833333333333333 * f(i-2) &
            +0.583333333333333 * f(i-1) &
            +0.583333333333333 * f(i+0) &
            -0.0833333333333333 * f(i+1)
       fr3 = +0.0833333333333333 * f(i-3) &
            -0.416666666666667 * f(i-2) &
            +1.08333333333333 * f(i-1) &
            +0.25 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3
    end do
  end subroutine reconstruct_left_k4
  subroutine weights_right_k4(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:3)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:3,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3
    do i=4, n-4
       accumulator = 0.0
       omega0 = +0.114285714285714 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.514285714285714 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.342857142857143 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.0285714285714286 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
    end do
  end subroutine weights_right_k4
  subroutine reconstruct_right_k4(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:3,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3
    real(kind=8) :: fr0, fr1, fr2, fr3
    do i=4, n-4
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       fr0 = +0.25 * f(i+0) &
            +1.08333333333333 * f(i+1) &
            -0.416666666666667 * f(i+2) &
            +0.0833333333333333 * f(i+3)
       fr1 = -0.0833333333333333 * f(i-1) &
            +0.583333333333333 * f(i+0) &
            +0.583333333333333 * f(i+1) &
            -0.0833333333333333 * f(i+2)
       fr2 = +0.0833333333333333 * f(i-2) &
            -0.416666666666667 * f(i-1) &
            +1.08333333333333 * f(i+0) &
            +0.25 * f(i+1)
       fr3 = -0.25 * f(i-3) &
            +1.08333333333333 * f(i-2) &
            -1.91666666666667 * f(i-1) &
            +2.08333333333333 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3
    end do
  end subroutine reconstruct_right_k4
  subroutine smoothness_k5(f, n, sigma)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: sigma(n,0:4)
    integer :: i

    do i=5, n-5
       sigma(i,0) = +21.4123015873016 * f(i+0) * f(i+0) &
            -128.869246031746 * f(i+0) * f(i+1) &
            +150.560119047619 * f(i+0) * f(i+2) &
            -81.644246031746 * f(i+0) * f(i+3) &
            +17.1287698412698 * f(i+0) * f(i+4) &
            +202.492658730159 * f(i+1) * f(i+1) &
            -488.507142857143 * f(i+1) * f(i+2) &
            +269.535317460317 * f(i+1) * f(i+3) &
            -57.144246031746 * f(i+1) * f(i+4) &
            +301.86369047619 * f(i+2) * f(i+2) &
            -338.173809523809 * f(i+2) * f(i+3) &
            +72.3934523809524 * f(i+2) * f(i+4) &
            +95.8259920634921 * f(i+3) * f(i+3) &
            -41.369246031746 * f(i+3) * f(i+4) &
            +4.49563492063492 * f(i+4) * f(i+4)
       sigma(i,1) = +4.49563492063492 * f(i-1) * f(i-1) &
            -27.8275793650794 * f(i-1) * f(i+0) &
            +32.7684523809524 * f(i-1) * f(i+1) &
            -17.519246031746 * f(i-1) * f(i+2) &
            +3.58710317460317 * f(i-1) * f(i+3) &
            +48.1593253968254 * f(i+0) * f(i+0) &
            -121.42380952381 * f(i+0) * f(i+1) &
            +66.8686507936508 * f(i+0) * f(i+2) &
            -13.9359126984127 * f(i+0) * f(i+3) &
            +80.6136904761905 * f(i+1) * f(i+1) &
            -92.2571428571429 * f(i+1) * f(i+2) &
            +19.685119047619 * f(i+1) * f(i+3) &
            +27.4926587301587 * f(i+2) * f(i+2) &
            -12.0775793650794 * f(i+2) * f(i+3) &
            +1.37063492063492 * f(i+3) * f(i+3)
       sigma(i,2) = +1.37063492063492 * f(i-2) * f(i-2) &
            -10.119246031746 * f(i-2) * f(i-1) &
            +13.4767857142857 * f(i-2) * f(i+0) &
            -7.72757936507937 * f(i-2) * f(i+1) &
            +1.62876984126984 * f(i-2) * f(i+2) &
            +20.8259920634921 * f(i-1) * f(i-1) &
            -59.3404761904762 * f(i-1) * f(i+0) &
            +35.5353174603175 * f(i-1) * f(i+1) &
            -7.72757936507937 * f(i-1) * f(i+2) &
            +45.8636904761905 * f(i+0) * f(i+0) &
            -59.3404761904762 * f(i+0) * f(i+1) &
            +13.4767857142857 * f(i+0) * f(i+2) &
            +20.8259920634921 * f(i+1) * f(i+1) &
            -10.119246031746 * f(i+1) * f(i+2) &
            +1.37063492063492 * f(i+2) * f(i+2)
       sigma(i,3) = +1.37063492063492 * f(i-3) * f(i-3) &
            -12.0775793650794 * f(i-3) * f(i-2) &
            +19.685119047619 * f(i-3) * f(i-1) &
            -13.9359126984127 * f(i-3) * f(i+0) &
            +3.58710317460317 * f(i-3) * f(i+1) &
            +27.4926587301587 * f(i-2) * f(i-2) &
            -92.2571428571429 * f(i-2) * f(i-1) &
            +66.8686507936508 * f(i-2) * f(i+0) &
            -17.519246031746 * f(i-2) * f(i+1) &
            +80.6136904761905 * f(i-1) * f(i-1) &
            -121.42380952381 * f(i-1) * f(i+0) &
            +32.7684523809524 * f(i-1) * f(i+1) &
            +48.1593253968254 * f(i+0) * f(i+0) &
            -27.8275793650794 * f(i+0) * f(i+1) &
            +4.49563492063492 * f(i+1) * f(i+1)
       sigma(i,4) = +4.49563492063492 * f(i-4) * f(i-4) &
            -41.369246031746 * f(i-4) * f(i-3) &
            +72.3934523809524 * f(i-4) * f(i-2) &
            -57.144246031746 * f(i-4) * f(i-1) &
            +17.1287698412698 * f(i-4) * f(i+0) &
            +95.8259920634921 * f(i-3) * f(i-3) &
            -338.173809523809 * f(i-3) * f(i-2) &
            +269.535317460317 * f(i-3) * f(i-1) &
            -81.644246031746 * f(i-3) * f(i+0) &
            +301.86369047619 * f(i-2) * f(i-2) &
            -488.507142857143 * f(i-2) * f(i-1) &
            +150.560119047619 * f(i-2) * f(i+0) &
            +202.492658730159 * f(i-1) * f(i-1) &
            -128.869246031746 * f(i-1) * f(i+0) &
            +21.4123015873016 * f(i+0) * f(i+0)
    end do
  end subroutine smoothness_k5

  subroutine weights_left_k5(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:4)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:4,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4
    do i=5, n-5
       accumulator = 0.0
       omega0 = +0.00793650793650794 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.158730158730159 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.476190476190476 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.317460317460317 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega4 = +0.0396825396825397 / (1.e-36 + sigma(i,4)) / (1.e-36 + sigma(i,4))
       accumulator = accumulator + omega4
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega4 = omega4 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
       omega(i,1,4,1) = omega4
    end do
  end subroutine weights_left_k5
  subroutine reconstruct_left_k5(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:4,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4
    real(kind=8) :: fr0, fr1, fr2, fr3, fr4
    do i=5, n-5
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       omega4 = omega(i,1,4,1)
       fr0 = +2.28333333333333 * f(i+0) &
            -2.71666666666667 * f(i+1) &
            +2.28333333333333 * f(i+2) &
            -1.05 * f(i+3) &
            +0.2 * f(i+4)
       fr1 = +0.2 * f(i-1) &
            +1.28333333333333 * f(i+0) &
            -0.716666666666667 * f(i+1) &
            +0.283333333333333 * f(i+2) &
            -0.05 * f(i+3)
       fr2 = -0.05 * f(i-2) &
            +0.45 * f(i-1) &
            +0.783333333333333 * f(i+0) &
            -0.216666666666667 * f(i+1) &
            +0.0333333333333333 * f(i+2)
       fr3 = +0.0333333333333333 * f(i-3) &
            -0.216666666666667 * f(i-2) &
            +0.783333333333333 * f(i-1) &
            +0.45 * f(i+0) &
            -0.05 * f(i+1)
       fr4 = -0.05 * f(i-4) &
            +0.283333333333333 * f(i-3) &
            -0.716666666666667 * f(i-2) &
            +1.28333333333333 * f(i-1) &
            +0.2 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3 + &
            fr4 * omega4
    end do
  end subroutine reconstruct_left_k5
  subroutine weights_right_k5(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:4)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:4,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4
    do i=5, n-5
       accumulator = 0.0
       omega0 = +0.0396825396825397 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.317460317460317 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.476190476190476 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.158730158730159 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega4 = +0.00793650793650794 / (1.e-36 + sigma(i,4)) / (1.e-36 + sigma(i,4))
       accumulator = accumulator + omega4
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega4 = omega4 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
       omega(i,1,4,1) = omega4
    end do
  end subroutine weights_right_k5
  subroutine reconstruct_right_k5(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:4,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4
    real(kind=8) :: fr0, fr1, fr2, fr3, fr4
    do i=5, n-5
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       omega4 = omega(i,1,4,1)
       fr0 = +0.2 * f(i+0) &
            +1.28333333333333 * f(i+1) &
            -0.716666666666667 * f(i+2) &
            +0.283333333333333 * f(i+3) &
            -0.05 * f(i+4)
       fr1 = -0.05 * f(i-1) &
            +0.45 * f(i+0) &
            +0.783333333333333 * f(i+1) &
            -0.216666666666667 * f(i+2) &
            +0.0333333333333333 * f(i+3)
       fr2 = +0.0333333333333333 * f(i-2) &
            -0.216666666666667 * f(i-1) &
            +0.783333333333333 * f(i+0) &
            +0.45 * f(i+1) &
            -0.05 * f(i+2)
       fr3 = -0.05 * f(i-3) &
            +0.283333333333333 * f(i-2) &
            -0.716666666666667 * f(i-1) &
            +1.28333333333333 * f(i+0) &
            +0.2 * f(i+1)
       fr4 = +0.2 * f(i-4) &
            -1.05 * f(i-3) &
            +2.28333333333333 * f(i-2) &
            -2.71666666666667 * f(i-1) &
            +2.28333333333333 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3 + &
            fr4 * omega4
    end do
  end subroutine reconstruct_right_k5
  subroutine smoothness_k6(f, n, sigma)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: sigma(n,0:5)
    integer :: i

    do i=6, n-6
       sigma(i,0) = +50.8449983465608 * f(i+0) * f(i+0) &
            -392.364947089947 * f(i+0) * f(i+1) &
            +630.016005291005 * f(i+0) * f(i+2) &
            -524.091633597884 * f(i+0) * f(i+3) &
            +223.711722883598 * f(i+0) * f(i+4) &
            -38.9611441798942 * f(i+0) * f(i+5) &
            +784.153745039683 * f(i+1) * f(i+1) &
            -2577.47390873016 * f(i+1) * f(i+2) &
            +2173.45958994709 * f(i+1) * f(i+3) &
            -935.902678571429 * f(i+1) * f(i+4) &
            +163.974454365079 * f(i+1) * f(i+5) &
            +2153.15287698413 * f(i+2) * f(i+2) &
            -3670.6671957672 * f(i+2) * f(i+3) &
            +1592.23273809524 * f(i+2) * f(i+4) &
            -280.413392857143 * f(i+2) * f(i+5) &
            +1577.03019179894 * f(i+3) * f(i+3) &
            -1376.16603835979 * f(i+3) * f(i+4) &
            +243.404894179894 * f(i+3) * f(i+5) &
            +301.592981150794 * f(i+4) * f(i+4) &
            -107.061706349206 * f(i+4) * f(i+5) &
            +9.52844742063492 * f(i+5) * f(i+5)
       sigma(i,1) = +9.52844742063492 * f(i-1) * f(i-1) &
            -75.3802248677249 * f(i-1) * f(i+0) &
            +121.878968253968 * f(i-1) * f(i+1) &
            -100.724503968254 * f(i-1) * f(i+2) &
            +42.4485284391534 * f(i-1) * f(i+3) &
            -7.2796626984127 * f(i-1) * f(i+4) &
            +160.102240410053 * f(i+0) * f(i+0) &
            -539.221593915344 * f(i+0) * f(i+1) &
            +455.140145502646 * f(i+0) * f(i+2) &
            -194.365641534392 * f(i+0) * f(i+3) &
            +33.622833994709 * f(i+0) * f(i+4) &
            +468.437599206349 * f(i+1) * f(i+1) &
            -808.852380952381 * f(i+1) * f(i+2) &
            +350.570701058201 * f(i+1) * f(i+3) &
            -61.2508928571429 * f(i+1) * f(i+4) &
            +356.263988095238 * f(i+2) * f(i+2) &
            -313.436871693122 * f(i+2) * f(i+3) &
            +55.3456349206349 * f(i+2) * f(i+4) &
            +69.8574487433862 * f(i+3) * f(i+3) &
            -24.9316137566138 * f(i+3) * f(i+4) &
            +2.2468501984127 * f(i+4) * f(i+4)
       sigma(i,2) = +2.2468501984127 * f(i-2) * f(i-2) &
            -19.6825396825397 * f(i-2) * f(i-1) &
            +33.782671957672 * f(i-2) * f(i+0) &
            -28.6231150793651 * f(i-2) * f(i+1) &
            +12.059871031746 * f(i-2) * f(i+2) &
            -2.03058862433862 * f(i-2) * f(i+3) &
            +46.7370783730159 * f(i-1) * f(i-1) &
            -168.881316137566 * f(i-1) * f(i+0) &
            +148.024404761905 * f(i-1) * f(i+1) &
            -63.8887896825397 * f(i-1) * f(i+2) &
            +10.954083994709 * f(i-1) * f(i+3) &
            +161.301025132275 * f(i+0) * f(i+0) &
            -296.11164021164 * f(i+0) * f(i+1) &
            +131.695701058201 * f(i+0) * f(i+2) &
            -23.0874669312169 * f(i+0) * f(i+3) &
            +142.159821428571 * f(i+1) * f(i+1) &
            -131.286408730159 * f(i+1) * f(i+2) &
            +23.6771164021164 * f(i+1) * f(i+3) &
            +31.6207589285714 * f(i+2) * f(i+2) &
            -11.8218915343915 * f(i+2) * f(i+3) &
            +1.15437334656085 * f(i+3) * f(i+3)
       sigma(i,3) = +1.15437334656085 * f(i-3) * f(i-3) &
            -11.8218915343915 * f(i-3) * f(i-2) &
            +23.6771164021164 * f(i-3) * f(i-1) &
            -23.0874669312169 * f(i-3) * f(i+0) &
            +10.954083994709 * f(i-3) * f(i+1) &
            -2.03058862433862 * f(i-3) * f(i+2) &
            +31.6207589285714 * f(i-2) * f(i-2) &
            -131.286408730159 * f(i-2) * f(i-1) &
            +131.695701058201 * f(i-2) * f(i+0) &
            -63.8887896825397 * f(i-2) * f(i+1) &
            +12.059871031746 * f(i-2) * f(i+2) &
            +142.159821428571 * f(i-1) * f(i-1) &
            -296.11164021164 * f(i-1) * f(i+0) &
            +148.024404761905 * f(i-1) * f(i+1) &
            -28.6231150793651 * f(i-1) * f(i+2) &
            +161.301025132275 * f(i+0) * f(i+0) &
            -168.881316137566 * f(i+0) * f(i+1) &
            +33.782671957672 * f(i+0) * f(i+2) &
            +46.7370783730159 * f(i+1) * f(i+1) &
            -19.6825396825397 * f(i+1) * f(i+2) &
            +2.2468501984127 * f(i+2) * f(i+2)
       sigma(i,4) = +2.2468501984127 * f(i-4) * f(i-4) &
            -24.9316137566138 * f(i-4) * f(i-3) &
            +55.3456349206349 * f(i-4) * f(i-2) &
            -61.2508928571429 * f(i-4) * f(i-1) &
            +33.622833994709 * f(i-4) * f(i+0) &
            -7.2796626984127 * f(i-4) * f(i+1) &
            +69.8574487433862 * f(i-3) * f(i-3) &
            -313.436871693122 * f(i-3) * f(i-2) &
            +350.570701058201 * f(i-3) * f(i-1) &
            -194.365641534392 * f(i-3) * f(i+0) &
            +42.4485284391534 * f(i-3) * f(i+1) &
            +356.263988095238 * f(i-2) * f(i-2) &
            -808.852380952381 * f(i-2) * f(i-1) &
            +455.140145502646 * f(i-2) * f(i+0) &
            -100.724503968254 * f(i-2) * f(i+1) &
            +468.437599206349 * f(i-1) * f(i-1) &
            -539.221593915344 * f(i-1) * f(i+0) &
            +121.878968253968 * f(i-1) * f(i+1) &
            +160.102240410053 * f(i+0) * f(i+0) &
            -75.3802248677249 * f(i+0) * f(i+1) &
            +9.52844742063492 * f(i+1) * f(i+1)
       sigma(i,5) = +9.52844742063492 * f(i-5) * f(i-5) &
            -107.061706349206 * f(i-5) * f(i-4) &
            +243.404894179894 * f(i-5) * f(i-3) &
            -280.413392857143 * f(i-5) * f(i-2) &
            +163.974454365079 * f(i-5) * f(i-1) &
            -38.9611441798942 * f(i-5) * f(i+0) &
            +301.592981150794 * f(i-4) * f(i-4) &
            -1376.16603835979 * f(i-4) * f(i-3) &
            +1592.23273809524 * f(i-4) * f(i-2) &
            -935.902678571429 * f(i-4) * f(i-1) &
            +223.711722883598 * f(i-4) * f(i+0) &
            +1577.03019179894 * f(i-3) * f(i-3) &
            -3670.6671957672 * f(i-3) * f(i-2) &
            +2173.45958994709 * f(i-3) * f(i-1) &
            -524.091633597884 * f(i-3) * f(i+0) &
            +2153.15287698413 * f(i-2) * f(i-2) &
            -2577.47390873016 * f(i-2) * f(i-1) &
            +630.016005291005 * f(i-2) * f(i+0) &
            +784.153745039683 * f(i-1) * f(i-1) &
            -392.364947089947 * f(i-1) * f(i+0) &
            +50.8449983465608 * f(i+0) * f(i+0)
    end do
  end subroutine smoothness_k6

  subroutine weights_left_k6(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:5)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:5,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5
    do i=6, n-6
       accumulator = 0.0
       omega0 = +0.00216450216450216 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.0649350649350649 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.324675324675325 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.432900432900433 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega4 = +0.162337662337662 / (1.e-36 + sigma(i,4)) / (1.e-36 + sigma(i,4))
       accumulator = accumulator + omega4
       omega5 = +0.012987012987013 / (1.e-36 + sigma(i,5)) / (1.e-36 + sigma(i,5))
       accumulator = accumulator + omega5
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega4 = omega4 / accumulator
       omega5 = omega5 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
       omega(i,1,4,1) = omega4
       omega(i,1,5,1) = omega5
    end do
  end subroutine weights_left_k6
  subroutine reconstruct_left_k6(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:5,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5
    real(kind=8) :: fr0, fr1, fr2, fr3, fr4, fr5
    do i=6, n-6
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       omega4 = omega(i,1,4,1)
       omega5 = omega(i,1,5,1)
       fr0 = +2.45 * f(i+0) &
            -3.55 * f(i+1) &
            +3.95 * f(i+2) &
            -2.71666666666667 * f(i+3) &
            +1.03333333333333 * f(i+4) &
            -0.166666666666667 * f(i+5)
       fr1 = +0.166666666666667 * f(i-1) &
            +1.45 * f(i+0) &
            -1.05 * f(i+1) &
            +0.616666666666667 * f(i+2) &
            -0.216666666666667 * f(i+3) &
            +0.0333333333333333 * f(i+4)
       fr2 = -0.0333333333333333 * f(i-2) &
            +0.366666666666667 * f(i-1) &
            +0.95 * f(i+0) &
            -0.383333333333333 * f(i+1) &
            +0.116666666666667 * f(i+2) &
            -0.0166666666666667 * f(i+3)
       fr3 = +0.0166666666666667 * f(i-3) &
            -0.133333333333333 * f(i-2) &
            +0.616666666666667 * f(i-1) &
            +0.616666666666667 * f(i+0) &
            -0.133333333333333 * f(i+1) &
            +0.0166666666666667 * f(i+2)
       fr4 = -0.0166666666666667 * f(i-4) &
            +0.116666666666667 * f(i-3) &
            -0.383333333333333 * f(i-2) &
            +0.95 * f(i-1) &
            +0.366666666666667 * f(i+0) &
            -0.0333333333333333 * f(i+1)
       fr5 = +0.0333333333333333 * f(i-5) &
            -0.216666666666667 * f(i-4) &
            +0.616666666666667 * f(i-3) &
            -1.05 * f(i-2) &
            +1.45 * f(i-1) &
            +0.166666666666667 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3 + &
            fr4 * omega4 + &
            fr5 * omega5
    end do
  end subroutine reconstruct_left_k6
  subroutine weights_right_k6(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:5)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:5,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5
    do i=6, n-6
       accumulator = 0.0
       omega0 = +0.012987012987013 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.162337662337662 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.432900432900433 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.324675324675325 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega4 = +0.0649350649350649 / (1.e-36 + sigma(i,4)) / (1.e-36 + sigma(i,4))
       accumulator = accumulator + omega4
       omega5 = +0.00216450216450216 / (1.e-36 + sigma(i,5)) / (1.e-36 + sigma(i,5))
       accumulator = accumulator + omega5
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega4 = omega4 / accumulator
       omega5 = omega5 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
       omega(i,1,4,1) = omega4
       omega(i,1,5,1) = omega5
    end do
  end subroutine weights_right_k6
  subroutine reconstruct_right_k6(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:5,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5
    real(kind=8) :: fr0, fr1, fr2, fr3, fr4, fr5
    do i=6, n-6
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       omega4 = omega(i,1,4,1)
       omega5 = omega(i,1,5,1)
       fr0 = +0.166666666666667 * f(i+0) &
            +1.45 * f(i+1) &
            -1.05 * f(i+2) &
            +0.616666666666667 * f(i+3) &
            -0.216666666666667 * f(i+4) &
            +0.0333333333333333 * f(i+5)
       fr1 = -0.0333333333333333 * f(i-1) &
            +0.366666666666667 * f(i+0) &
            +0.95 * f(i+1) &
            -0.383333333333333 * f(i+2) &
            +0.116666666666667 * f(i+3) &
            -0.0166666666666667 * f(i+4)
       fr2 = +0.0166666666666667 * f(i-2) &
            -0.133333333333333 * f(i-1) &
            +0.616666666666667 * f(i+0) &
            +0.616666666666667 * f(i+1) &
            -0.133333333333333 * f(i+2) &
            +0.0166666666666667 * f(i+3)
       fr3 = -0.0166666666666667 * f(i-3) &
            +0.116666666666667 * f(i-2) &
            -0.383333333333333 * f(i-1) &
            +0.95 * f(i+0) &
            +0.366666666666667 * f(i+1) &
            -0.0333333333333333 * f(i+2)
       fr4 = +0.0333333333333333 * f(i-4) &
            -0.216666666666667 * f(i-3) &
            +0.616666666666667 * f(i-2) &
            -1.05 * f(i-1) &
            +1.45 * f(i+0) &
            +0.166666666666667 * f(i+1)
       fr5 = -0.166666666666667 * f(i-5) &
            +1.03333333333333 * f(i-4) &
            -2.71666666666667 * f(i-3) &
            +3.95 * f(i-2) &
            -3.55 * f(i-1) &
            +2.45 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3 + &
            fr4 * omega4 + &
            fr5 * omega5
    end do
  end subroutine reconstruct_right_k6
  subroutine smoothness_k7(f, n, sigma)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: sigma(n,0:6)
    integer :: i

    do i=7, n-7
       sigma(i,0) = +119.876965822244 * f(i+0) * f(i+0) &
            -1140.52691383077 * f(i+0) * f(i+1) &
            +2345.30742098565 * f(i+0) * f(i+2) &
            -2648.53826913313 * f(i+0) * f(i+3) &
            +1719.4383816338 * f(i+0) * f(i+4) &
            -605.48066052389 * f(i+0) * f(i+5) &
            +90.0461092238523 * f(i+0) * f(i+6) &
            +2787.97471636003 * f(i+1) * f(i+1) &
            -11665.8977583874 * f(i+1) * f(i+2) &
            +13315.7065438111 * f(i+1) * f(i+3) &
            -8706.93798656205 * f(i+1) * f(i+4) &
            +3081.76169462482 * f(i+1) * f(i+5) &
            -460.055012375742 * f(i+1) * f(i+6) &
            +12350.3314303752 * f(i+2) * f(i+2) &
            -28424.0145572791 * f(i+2) * f(i+3) &
            +18693.1184907107 * f(i+2) * f(i+4) &
            -6644.20048656205 * f(i+2) * f(i+5) &
            +995.024029781946 * f(i+2) * f(i+6) &
            +16453.1759122308 * f(i+3) * f(i+3) &
            -21738.2182609828 * f(i+3) * f(i+4) &
            +7752.80284010742 * f(i+3) * f(i+5) &
            -1164.09012098498 * f(i+3) * f(i+6) &
            +7205.30018037518 * f(i+4) * f(i+4) &
            -5153.46025838745 * f(i+4) * f(i+5) &
            +775.459272837502 * f(i+4) * f(i+6) &
            +923.494716360029 * f(i+5) * f(i+5) &
            -278.412561978916 * f(i+5) * f(i+6) &
            +21.0141417481695 * f(i+6) * f(i+6)
       sigma(i,1) = +21.0141417481695 * f(i-1) * f(i-1) &
            -204.151875250521 * f(i-1) * f(i+0) &
            +422.538941047379 * f(i-1) * f(i+1) &
            -475.96589258992 * f(i-1) * f(i+2) &
            +306.899801386885 * f(i-1) * f(i+3) &
            -107.134680585618 * f(i-1) * f(i+4) &
            +15.7854224954572 * f(i-1) * f(i+5) &
            +519.247146915584 * f(i+0) * f(i+0) &
            -2207.33120746152 * f(i+0) * f(i+1) &
            +2525.45484628026 * f(i+0) * f(i+2) &
            -1645.22305600649 * f(i+0) * f(i+3) &
            +578.412852032227 * f(i+0) * f(i+4) &
            -85.6558534251243 * f(i+0) * f(i+5) &
            +2394.05596741222 * f(i+1) * f(i+1) &
            -5559.25606962482 * f(i+1) * f(i+2) &
            +3658.67693978475 * f(i+1) * f(i+3) &
            -1295.61101896946 * f(i+1) * f(i+4) &
            +192.87048039923 * f(i+1) * f(i+5) &
            +3266.81402951472 * f(i+2) * f(i+2) &
            -4339.66656345198 * f(i+2) * f(i+3) &
            +1547.32768578644 * f(i+2) * f(i+4) &
            -231.522065429427 * f(i+2) * f(i+5) &
            +1452.34531926407 * f(i+3) * f(i+3) &
            -1042.03954079485 * f(i+3) * f(i+4) &
            +156.661780553551 * f(i+3) * f(i+5) &
            +187.891961730399 * f(i+4) * f(i+4) &
            -56.7392209295334 * f(i+4) * f(i+5) &
            +4.29972816792261 * f(i+5) * f(i+5)
       sigma(i,2) = +4.29972816792261 * f(i-2) * f(i-2) &
            -44.4107718554594 * f(i-2) * f(i-1) &
            +94.9327296276255 * f(i-2) * f(i+0) &
            -108.110491355352 * f(i-2) * f(i+1) &
            +69.4589063251563 * f(i-2) * f(i+2) &
            -23.9268024991983 * f(i-2) * f(i+3) &
            +3.45697342138314 * f(i-2) * f(i+4) &
            +121.202864508177 * f(i-1) * f(i-1) &
            -537.187110239298 * f(i-1) * f(i+0) &
            +626.822593193843 * f(i-1) * f(i+1) &
            -409.688449525012 * f(i-1) * f(i+2) &
            +142.893546476671 * f(i-1) * f(i+3) &
            -20.8355370670996 * f(i-1) * f(i+4) &
            +616.654347041847 * f(i+0) * f(i+0) &
            -1479.69665604457 * f(i+0) * f(i+1) &
            +986.137009229197 * f(i+0) * f(i+2) &
            -348.912986562049 * f(i+0) * f(i+3) &
            +51.4183199054032 * f(i+0) * f(i+4) &
            +910.756159144354 * f(i+1) * f(i+1) &
            -1239.85097703223 * f(i+1) * f(i+2) &
            +445.834938872856 * f(i+1) * f(i+3) &
            -66.5117259232537 * f(i+1) * f(i+4) &
            +430.708745189995 * f(i+2) * f(i+2) &
            -315.141276905964 * f(i+2) * f(i+3) &
            +47.66729752886 * f(i+2) * f(i+4) &
            +58.6280496933622 * f(i+3) * f(i+3) &
            -18.0035187690396 * f(i+3) * f(i+4) &
            +1.40409545187323 * f(i+4) * f(i+4)
       sigma(i,3) = +1.40409545187323 * f(i-3) * f(i-3) &
            -16.2003629048421 * f(i-3) * f(i-2) &
            +38.1364719115761 * f(i-3) * f(i-1) &
            -46.8683617257228 * f(i-3) * f(i+0) &
            +31.7749557078724 * f(i-3) * f(i+1) &
            -11.3047114498156 * f(i-3) * f(i+2) &
            +1.65381755718561 * f(i-3) * f(i+3) &
            +48.9015913600289 * f(i-2) * f(i-2) &
            -238.769633387446 * f(i-2) * f(i-1) &
            +302.017191959275 * f(i-2) * f(i+0) &
            -209.541111562049 * f(i-2) * f(i+1) &
            +75.9954446248196 * f(i-2) * f(i+2) &
            -11.3047114498156 * f(i-2) * f(i+3) &
            +302.86268037518 * f(i-1) * f(i-1) &
            -792.178909130992 * f(i-1) * f(i+0) &
            +564.852865710678 * f(i-1) * f(i+1) &
            -209.541111562049 * f(i-1) * f(i+2) &
            +31.7749557078724 * f(i-1) * f(i+3) &
            +537.03007889744 * f(i+0) * f(i+0) &
            -792.178909130992 * f(i+0) * f(i+1) &
            +302.017191959275 * f(i+0) * f(i+2) &
            -46.8683617257228 * f(i+0) * f(i+3) &
            +302.86268037518 * f(i+1) * f(i+1) &
            -238.769633387446 * f(i+1) * f(i+2) &
            +38.1364719115761 * f(i+1) * f(i+3) &
            +48.9015913600289 * f(i+2) * f(i+2) &
            -16.2003629048421 * f(i+2) * f(i+3) &
            +1.40409545187323 * f(i+3) * f(i+3)
       sigma(i,4) = +1.40409545187323 * f(i-4) * f(i-4) &
            -18.0035187690396 * f(i-4) * f(i-3) &
            +47.66729752886 * f(i-4) * f(i-2) &
            -66.5117259232537 * f(i-4) * f(i-1) &
            +51.4183199054032 * f(i-4) * f(i+0) &
            -20.8355370670996 * f(i-4) * f(i+1) &
            +3.45697342138314 * f(i-4) * f(i+2) &
            +58.6280496933622 * f(i-3) * f(i-3) &
            -315.141276905964 * f(i-3) * f(i-2) &
            +445.834938872856 * f(i-3) * f(i-1) &
            -348.912986562049 * f(i-3) * f(i+0) &
            +142.893546476671 * f(i-3) * f(i+1) &
            -23.9268024991983 * f(i-3) * f(i+2) &
            +430.708745189995 * f(i-2) * f(i-2) &
            -1239.85097703223 * f(i-2) * f(i-1) &
            +986.137009229197 * f(i-2) * f(i+0) &
            -409.688449525012 * f(i-2) * f(i+1) &
            +69.4589063251563 * f(i-2) * f(i+2) &
            +910.756159144354 * f(i-1) * f(i-1) &
            -1479.69665604457 * f(i-1) * f(i+0) &
            +626.822593193843 * f(i-1) * f(i+1) &
            -108.110491355352 * f(i-1) * f(i+2) &
            +616.654347041847 * f(i+0) * f(i+0) &
            -537.187110239298 * f(i+0) * f(i+1) &
            +94.9327296276255 * f(i+0) * f(i+2) &
            +121.202864508177 * f(i+1) * f(i+1) &
            -44.4107718554594 * f(i+1) * f(i+2) &
            +4.29972816792261 * f(i+2) * f(i+2)
       sigma(i,5) = +4.29972816792261 * f(i-5) * f(i-5) &
            -56.7392209295334 * f(i-5) * f(i-4) &
            +156.661780553551 * f(i-5) * f(i-3) &
            -231.522065429427 * f(i-5) * f(i-2) &
            +192.87048039923 * f(i-5) * f(i-1) &
            -85.6558534251243 * f(i-5) * f(i+0) &
            +15.7854224954572 * f(i-5) * f(i+1) &
            +187.891961730399 * f(i-4) * f(i-4) &
            -1042.03954079485 * f(i-4) * f(i-3) &
            +1547.32768578644 * f(i-4) * f(i-2) &
            -1295.61101896946 * f(i-4) * f(i-1) &
            +578.412852032227 * f(i-4) * f(i+0) &
            -107.134680585618 * f(i-4) * f(i+1) &
            +1452.34531926407 * f(i-3) * f(i-3) &
            -4339.66656345198 * f(i-3) * f(i-2) &
            +3658.67693978475 * f(i-3) * f(i-1) &
            -1645.22305600649 * f(i-3) * f(i+0) &
            +306.899801386885 * f(i-3) * f(i+1) &
            +3266.81402951472 * f(i-2) * f(i-2) &
            -5559.25606962482 * f(i-2) * f(i-1) &
            +2525.45484628026 * f(i-2) * f(i+0) &
            -475.96589258992 * f(i-2) * f(i+1) &
            +2394.05596741222 * f(i-1) * f(i-1) &
            -2207.33120746152 * f(i-1) * f(i+0) &
            +422.538941047379 * f(i-1) * f(i+1) &
            +519.247146915584 * f(i+0) * f(i+0) &
            -204.151875250521 * f(i+0) * f(i+1) &
            +21.0141417481695 * f(i+1) * f(i+1)
       sigma(i,6) = +21.0141417481695 * f(i-6) * f(i-6) &
            -278.412561978916 * f(i-6) * f(i-5) &
            +775.459272837502 * f(i-6) * f(i-4) &
            -1164.09012098498 * f(i-6) * f(i-3) &
            +995.024029781946 * f(i-6) * f(i-2) &
            -460.055012375742 * f(i-6) * f(i-1) &
            +90.0461092238523 * f(i-6) * f(i+0) &
            +923.494716360029 * f(i-5) * f(i-5) &
            -5153.46025838745 * f(i-5) * f(i-4) &
            +7752.80284010742 * f(i-5) * f(i-3) &
            -6644.20048656205 * f(i-5) * f(i-2) &
            +3081.76169462482 * f(i-5) * f(i-1) &
            -605.48066052389 * f(i-5) * f(i+0) &
            +7205.30018037518 * f(i-4) * f(i-4) &
            -21738.2182609828 * f(i-4) * f(i-3) &
            +18693.1184907107 * f(i-4) * f(i-2) &
            -8706.93798656205 * f(i-4) * f(i-1) &
            +1719.4383816338 * f(i-4) * f(i+0) &
            +16453.1759122308 * f(i-3) * f(i-3) &
            -28424.0145572791 * f(i-3) * f(i-2) &
            +13315.7065438111 * f(i-3) * f(i-1) &
            -2648.53826913313 * f(i-3) * f(i+0) &
            +12350.3314303752 * f(i-2) * f(i-2) &
            -11665.8977583874 * f(i-2) * f(i-1) &
            +2345.30742098565 * f(i-2) * f(i+0) &
            +2787.97471636003 * f(i-1) * f(i-1) &
            -1140.52691383077 * f(i-1) * f(i+0) &
            +119.876965822244 * f(i+0) * f(i+0)
    end do
  end subroutine smoothness_k7

  subroutine weights_left_k7(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:6)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:6,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6
    do i=7, n-7
       accumulator = 0.0
       omega0 = +0.000582750582750583 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.0244755244755245 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.183566433566434 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.407925407925408 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega4 = +0.305944055944056 / (1.e-36 + sigma(i,4)) / (1.e-36 + sigma(i,4))
       accumulator = accumulator + omega4
       omega5 = +0.0734265734265734 / (1.e-36 + sigma(i,5)) / (1.e-36 + sigma(i,5))
       accumulator = accumulator + omega5
       omega6 = +0.00407925407925408 / (1.e-36 + sigma(i,6)) / (1.e-36 + sigma(i,6))
       accumulator = accumulator + omega6
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega4 = omega4 / accumulator
       omega5 = omega5 / accumulator
       omega6 = omega6 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
       omega(i,1,4,1) = omega4
       omega(i,1,5,1) = omega5
       omega(i,1,6,1) = omega6
    end do
  end subroutine weights_left_k7
  subroutine reconstruct_left_k7(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:6,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6
    real(kind=8) :: fr0, fr1, fr2, fr3, fr4, fr5, fr6
    do i=7, n-7
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       omega4 = omega(i,1,4,1)
       omega5 = omega(i,1,5,1)
       omega6 = omega(i,1,6,1)
       fr0 = +2.59285714285714 * f(i+0) &
            -4.40714285714286 * f(i+1) &
            +6.09285714285714 * f(i+2) &
            -5.57380952380952 * f(i+3) &
            +3.17619047619048 * f(i+4) &
            -1.02380952380952 * f(i+5) &
            +0.142857142857143 * f(i+6)
       fr1 = +0.142857142857143 * f(i-1) &
            +1.59285714285714 * f(i+0) &
            -1.40714285714286 * f(i+1) &
            +1.09285714285714 * f(i+2) &
            -0.573809523809524 * f(i+3) &
            +0.176190476190476 * f(i+4) &
            -0.0238095238095238 * f(i+5)
       fr2 = -0.0238095238095238 * f(i-2) &
            +0.30952380952381 * f(i-1) &
            +1.09285714285714 * f(i+0) &
            -0.573809523809524 * f(i+1) &
            +0.25952380952381 * f(i+2) &
            -0.0738095238095238 * f(i+3) &
            +0.00952380952380952 * f(i+4)
       fr3 = +0.00952380952380952 * f(i-3) &
            -0.0904761904761905 * f(i-2) &
            +0.509523809523809 * f(i-1) &
            +0.759523809523809 * f(i+0) &
            -0.24047619047619 * f(i+1) &
            +0.0595238095238095 * f(i+2) &
            -0.00714285714285714 * f(i+3)
       fr4 = -0.00714285714285714 * f(i-4) &
            +0.0595238095238095 * f(i-3) &
            -0.24047619047619 * f(i-2) &
            +0.759523809523809 * f(i-1) &
            +0.509523809523809 * f(i+0) &
            -0.0904761904761905 * f(i+1) &
            +0.00952380952380952 * f(i+2)
       fr5 = +0.00952380952380952 * f(i-5) &
            -0.0738095238095238 * f(i-4) &
            +0.25952380952381 * f(i-3) &
            -0.573809523809524 * f(i-2) &
            +1.09285714285714 * f(i-1) &
            +0.30952380952381 * f(i+0) &
            -0.0238095238095238 * f(i+1)
       fr6 = -0.0238095238095238 * f(i-6) &
            +0.176190476190476 * f(i-5) &
            -0.573809523809524 * f(i-4) &
            +1.09285714285714 * f(i-3) &
            -1.40714285714286 * f(i-2) &
            +1.59285714285714 * f(i-1) &
            +0.142857142857143 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3 + &
            fr4 * omega4 + &
            fr5 * omega5 + &
            fr6 * omega6
    end do
  end subroutine reconstruct_left_k7
  subroutine weights_right_k7(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:6)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:6,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6
    do i=7, n-7
       accumulator = 0.0
       omega0 = +0.00407925407925408 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.0734265734265734 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.305944055944056 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.407925407925408 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega4 = +0.183566433566434 / (1.e-36 + sigma(i,4)) / (1.e-36 + sigma(i,4))
       accumulator = accumulator + omega4
       omega5 = +0.0244755244755245 / (1.e-36 + sigma(i,5)) / (1.e-36 + sigma(i,5))
       accumulator = accumulator + omega5
       omega6 = +0.000582750582750583 / (1.e-36 + sigma(i,6)) / (1.e-36 + sigma(i,6))
       accumulator = accumulator + omega6
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega4 = omega4 / accumulator
       omega5 = omega5 / accumulator
       omega6 = omega6 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
       omega(i,1,4,1) = omega4
       omega(i,1,5,1) = omega5
       omega(i,1,6,1) = omega6
    end do
  end subroutine weights_right_k7
  subroutine reconstruct_right_k7(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:6,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6
    real(kind=8) :: fr0, fr1, fr2, fr3, fr4, fr5, fr6
    do i=7, n-7
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       omega4 = omega(i,1,4,1)
       omega5 = omega(i,1,5,1)
       omega6 = omega(i,1,6,1)
       fr0 = +0.142857142857143 * f(i+0) &
            +1.59285714285714 * f(i+1) &
            -1.40714285714286 * f(i+2) &
            +1.09285714285714 * f(i+3) &
            -0.573809523809524 * f(i+4) &
            +0.176190476190476 * f(i+5) &
            -0.0238095238095238 * f(i+6)
       fr1 = -0.0238095238095238 * f(i-1) &
            +0.30952380952381 * f(i+0) &
            +1.09285714285714 * f(i+1) &
            -0.573809523809524 * f(i+2) &
            +0.25952380952381 * f(i+3) &
            -0.0738095238095238 * f(i+4) &
            +0.00952380952380952 * f(i+5)
       fr2 = +0.00952380952380952 * f(i-2) &
            -0.0904761904761905 * f(i-1) &
            +0.509523809523809 * f(i+0) &
            +0.759523809523809 * f(i+1) &
            -0.24047619047619 * f(i+2) &
            +0.0595238095238095 * f(i+3) &
            -0.00714285714285714 * f(i+4)
       fr3 = -0.00714285714285714 * f(i-3) &
            +0.0595238095238095 * f(i-2) &
            -0.24047619047619 * f(i-1) &
            +0.759523809523809 * f(i+0) &
            +0.509523809523809 * f(i+1) &
            -0.0904761904761905 * f(i+2) &
            +0.00952380952380952 * f(i+3)
       fr4 = +0.00952380952380952 * f(i-4) &
            -0.0738095238095238 * f(i-3) &
            +0.25952380952381 * f(i-2) &
            -0.573809523809524 * f(i-1) &
            +1.09285714285714 * f(i+0) &
            +0.30952380952381 * f(i+1) &
            -0.0238095238095238 * f(i+2)
       fr5 = -0.0238095238095238 * f(i-5) &
            +0.176190476190476 * f(i-4) &
            -0.573809523809524 * f(i-3) &
            +1.09285714285714 * f(i-2) &
            -1.40714285714286 * f(i-1) &
            +1.59285714285714 * f(i+0) &
            +0.142857142857143 * f(i+1)
       fr6 = +0.142857142857143 * f(i-6) &
            -1.02380952380952 * f(i-5) &
            +3.17619047619048 * f(i-4) &
            -5.57380952380952 * f(i-3) &
            +6.09285714285714 * f(i-2) &
            -4.40714285714286 * f(i-1) &
            +2.59285714285714 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3 + &
            fr4 * omega4 + &
            fr5 * omega5 + &
            fr6 * omega6
    end do
  end subroutine reconstruct_right_k7
  subroutine smoothness_k8(f, n, sigma)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: sigma(n,0:7)
    integer :: i

    do i=8, n-8
       sigma(i,0) = +282.837600612977 * f(i+0) * f(i+0) &
            -3217.68658751845 * f(i+0) * f(i+1) &
            +8098.17149379591 * f(i+0) * f(i+2) &
            -11602.0447459035 * f(i+0) * f(i+3) &
            +10159.0931716562 * f(i+0) * f(i+4) &
            -5415.63416825253 * f(i+0) * f(i+5) &
            +1622.88916698582 * f(i+0) * f(i+6) &
            -210.463531989358 * f(i+0) * f(i+7) &
            +9343.02132742786 * f(i+1) * f(i+1) &
            -47645.872787025 * f(i+1) * f(i+2) &
            +68840.1294128134 * f(i+1) * f(i+3) &
            -60634.3990483282 * f(i+1) * f(i+4) &
            +32462.762767691 * f(i+1) * f(i+5) &
            -9759.93192303144 * f(i+1) * f(i+6) &
            +1268.95551054292 * f(i+1) * f(i+7) &
            +61294.8370166772 * f(i+2) * f(i+2) &
            -178245.759975439 * f(i+2) * f(i+3) &
            +157723.978487162 * f(i+2) * f(i+4) &
            -84736.2897924522 * f(i+2) * f(i+5) &
            +25544.3501239796 * f(i+2) * f(i+6) &
            -3328.25158337599 * f(i+2) * f(i+7) &
            +130199.124980546 * f(i+3) * f(i+3) &
            -231245.307361434 * f(i+3) * f(i+4) &
            +124579.67848041 * f(i+3) * f(i+5) &
            -37637.4314325875 * f(i+3) * f(i+6) &
            +4912.48566104661 * f(i+3) * f(i+7) &
            +102966.44021251 * f(i+4) * f(i+4) &
            -111189.450476982 * f(i+4) * f(i+5) &
            +33651.8387772039 * f(i+4) * f(i+6) &
            -4398.63397429859 * f(i+4) * f(i+7) &
            +30071.0784359481 * f(i+5) * f(i+5) &
            -18228.7647006052 * f(i+5) * f(i+6) &
            +2385.54101829437 * f(i+5) * f(i+7) &
            +2765.84444133604 * f(i+6) * f(i+6) &
            -724.638894617214 * f(i+6) * f(i+7) &
            +47.5028971986251 * f(i+7) * f(i+7)
       sigma(i,1) = +47.5028971986251 * f(i-1) * f(i-1) &
            -549.582823188643 * f(i-1) * f(i+0) &
            +1391.20673258008 * f(i-1) * f(i+1) &
            -1992.07290287002 * f(i-1) * f(i+2) &
            +1737.9199467609 * f(i-1) * f(i+3) &
            -921.690511947415 * f(i-1) * f(i+4) &
            +274.621224828637 * f(i-1) * f(i+5) &
            -35.407460560787 * f(i-1) * f(i+6) &
            +1639.31476541011 * f(i+0) * f(i+0) &
            -8454.36155245706 * f(i+0) * f(i+1) &
            +12248.796925352 * f(i+0) * f(i+2) &
            -10772.9570807356 * f(i+0) * f(i+3) &
            +5746.65947583143 * f(i+0) * f(i+4) &
            -1719.62507117959 * f(i+0) * f(i+5) &
            +222.440595557253 * f(i+0) * f(i+6) &
            +11054.5384359481 * f(i+1) * f(i+1) &
            -32362.4054769819 * f(i+1) * f(i+2) &
            +28675.0021841139 * f(i+1) * f(i+3) &
            -15380.2247924521 * f(i+1) * f(i+4) &
            +4621.402767691 * f(i+1) * f(i+5) &
            -599.696734390095 * f(i+1) * f(i+6) &
            +23881.8339625102 * f(i+2) * f(i+2) &
            -42591.6661577299 * f(i+2) * f(i+3) &
            +22956.5584871618 * f(i+2) * f(i+4) &
            -6924.03404832822 * f(i+2) * f(i+5) &
            +901.155248375756 * f(i+2) * f(i+6) &
            +19089.3249805465 * f(i+3) * f(i+3) &
            -20664.4461791424 * f(i+3) * f(i+4) &
            +6253.56570910972 * f(i+3) * f(i+5) &
            -816.068383469668 * f(i+3) * f(i+6) &
            +5612.02326667724 * f(i+4) * f(i+4) &
            -3406.48778702496 * f(i+4) * f(i+5) &
            +445.58477421919 * f(i+4) * f(i+6) &
            +518.201327427861 * f(i+5) * f(i+5) &
            -135.845449952311 * f(i+5) * f(i+6) &
            +8.91870511033141 * f(i+6) * f(i+6)
       sigma(i,2) = +8.91870511033141 * f(i-2) * f(i-2) &
            -107.291821204516 * f(i-2) * f(i-1) &
            +277.006890621306 * f(i-2) * f(i+0) &
            -399.198237967023 * f(i-2) * f(i+1) &
            +347.463467070642 * f(i-2) * f(i+2) &
            -182.82658888745 * f(i-2) * f(i+3) &
            +53.8627119593691 * f(i-2) * f(i+4) &
            -6.85383181299153 * f(i-2) * f(i+5) &
            +335.04033977354 * f(i-1) * f(i-1) &
            -1774.22905245706 * f(i-1) * f(i+0) &
            +2601.97484491219 * f(i-1) * f(i+1) &
            -2293.25840018007 * f(i-1) * f(i+2) &
            +1217.71486645643 * f(i-1) * f(i+3) &
            -361.183311920333 * f(i-1) * f(i+4) &
            +46.1921948462738 * f(i-1) * f(i+5) &
            +2403.24289630687 * f(i+0) * f(i+0) &
            -7175.23886432754 * f(i+0) * f(i+1) &
            +6406.93231432228 * f(i+0) * f(i+2) &
            -3435.42821837807 * f(i+0) * f(i+3) &
            +1026.47873509069 * f(i+0) * f(i+4) &
            -132.007597485334 * f(i+0) * f(i+5) &
            +5440.58053610203 * f(i+1) * f(i+1) &
            -9841.58822563111 * f(i+1) * f(i+2) &
            +5330.2740359658 * f(i+1) * f(i+3) &
            -1605.02809925414 * f(i+1) * f(i+4) &
            +207.643474097758 * f(i+1) * f(i+5) &
            +4502.6216168312 * f(i+2) * f(i+2) &
            -4924.83347080903 * f(i+2) * f(i+3) &
            +1494.60136979645 * f(i+2) * f(i+4) &
            -194.560288231573 * f(i+2) * f(i+5) &
            +1358.55473224437 * f(i+3) * f(i+3) &
            -830.843311716319 * f(i+3) * f(i+4) &
            +108.833222879904 * f(i+3) * f(i+5) &
            +127.914395039744 * f(i+4) * f(i+4) &
            -33.7168840352035 * f(i+4) * f(i+5) &
            +2.23485487058274 * f(i+5) * f(i+5)
       sigma(i,3) = +2.23485487058274 * f(i-3) * f(i-3) &
            -28.9038461163322 * f(i-3) * f(i-2) &
            +78.9596779063593 * f(i-3) * f(i-1) &
            -118.296148019933 * f(i-3) * f(i+0) &
            +105.236207783825 * f(i-3) * f(i+1) &
            -55.7434572736934 * f(i-3) * f(i+2) &
            +16.3186498727289 * f(i-3) * f(i+3) &
            -2.04079389412028 * f(i-3) * f(i+4) &
            +97.1187623236942 * f(i-2) * f(i-2) &
            -547.061953691627 * f(i-2) * f(i-1) &
            +839.561493253242 * f(i-2) * f(i+0) &
            -761.319673328215 * f(i-2) * f(i+1) &
            +409.596543732663 * f(i-2) * f(i+2) &
            -121.468497105518 * f(i-2) * f(i+3) &
            +15.3584086083991 * f(i-2) * f(i+4) &
            +793.785102614737 * f(i-1) * f(i-1) &
            -2499.75828562384 * f(i-1) * f(i+0) &
            +2315.13502362012 * f(i-1) * f(i+1) &
            -1267.31229245215 * f(i-1) * f(i+2) &
            +381.255607197169 * f(i-1) * f(i+3) &
            -48.78798218551 * f(i-1) * f(i+4) &
            +2019.32231127564 * f(i+0) * f(i+0) &
            -3827.93467624839 * f(i+0) * f(i+1) &
            +2136.14046247043 * f(i+0) * f(i+2) &
            -653.059881661548 * f(i+0) * f(i+3) &
            +84.7024132787544 * f(i+0) * f(i+4) &
            +1856.32621511438 * f(i+1) * f(i+1) &
            -2115.5956853152 * f(i+1) * f(i+2) &
            +658.5622523196 * f(i+1) * f(i+3) &
            -86.7358790604971 * f(i+1) * f(i+4) &
            +615.75035001057 * f(i+2) * f(i+2) &
            -390.9897931978 * f(i+2) * f(i+3) &
            +52.4035220146045 * f(i+2) * f(i+4) &
            +63.3507101439102 * f(i+3) * f(i+3) &
            -17.3197577124522 * f(i+3) * f(i+4) &
            +1.21003447541078 * f(i+4) * f(i+4)
       sigma(i,4) = +1.21003447541078 * f(i-4) * f(i-4) &
            -17.3197577124522 * f(i-4) * f(i-3) &
            +52.4035220146045 * f(i-4) * f(i-2) &
            -86.7358790604971 * f(i-4) * f(i-1) &
            +84.7024132787544 * f(i-4) * f(i+0) &
            -48.78798218551 * f(i-4) * f(i+1) &
            +15.3584086083991 * f(i-4) * f(i+2) &
            -2.04079389412028 * f(i-4) * f(i+3) &
            +63.3507101439102 * f(i-3) * f(i-3) &
            -390.9897931978 * f(i-3) * f(i-2) &
            +658.5622523196 * f(i-3) * f(i-1) &
            -653.059881661548 * f(i-3) * f(i+0) &
            +381.255607197169 * f(i-3) * f(i+1) &
            -121.468497105518 * f(i-3) * f(i+2) &
            +16.3186498727289 * f(i-3) * f(i+3) &
            +615.75035001057 * f(i-2) * f(i-2) &
            -2115.5956853152 * f(i-2) * f(i-1) &
            +2136.14046247043 * f(i-2) * f(i+0) &
            -1267.31229245215 * f(i-2) * f(i+1) &
            +409.596543732663 * f(i-2) * f(i+2) &
            -55.7434572736934 * f(i-2) * f(i+3) &
            +1856.32621511438 * f(i-1) * f(i-1) &
            -3827.93467624839 * f(i-1) * f(i+0) &
            +2315.13502362012 * f(i-1) * f(i+1) &
            -761.319673328215 * f(i-1) * f(i+2) &
            +105.236207783825 * f(i-1) * f(i+3) &
            +2019.32231127564 * f(i+0) * f(i+0) &
            -2499.75828562384 * f(i+0) * f(i+1) &
            +839.561493253242 * f(i+0) * f(i+2) &
            -118.296148019933 * f(i+0) * f(i+3) &
            +793.785102614737 * f(i+1) * f(i+1) &
            -547.061953691627 * f(i+1) * f(i+2) &
            +78.9596779063593 * f(i+1) * f(i+3) &
            +97.1187623236942 * f(i+2) * f(i+2) &
            -28.9038461163322 * f(i+2) * f(i+3) &
            +2.23485487058274 * f(i+3) * f(i+3)
       sigma(i,5) = +2.23485487058274 * f(i-5) * f(i-5) &
            -33.7168840352035 * f(i-5) * f(i-4) &
            +108.833222879904 * f(i-5) * f(i-3) &
            -194.560288231573 * f(i-5) * f(i-2) &
            +207.643474097758 * f(i-5) * f(i-1) &
            -132.007597485334 * f(i-5) * f(i+0) &
            +46.1921948462738 * f(i-5) * f(i+1) &
            -6.85383181299153 * f(i-5) * f(i+2) &
            +127.914395039744 * f(i-4) * f(i-4) &
            -830.843311716319 * f(i-4) * f(i-3) &
            +1494.60136979645 * f(i-4) * f(i-2) &
            -1605.02809925414 * f(i-4) * f(i-1) &
            +1026.47873509069 * f(i-4) * f(i+0) &
            -361.183311920333 * f(i-4) * f(i+1) &
            +53.8627119593691 * f(i-4) * f(i+2) &
            +1358.55473224437 * f(i-3) * f(i-3) &
            -4924.83347080903 * f(i-3) * f(i-2) &
            +5330.2740359658 * f(i-3) * f(i-1) &
            -3435.42821837807 * f(i-3) * f(i+0) &
            +1217.71486645643 * f(i-3) * f(i+1) &
            -182.82658888745 * f(i-3) * f(i+2) &
            +4502.6216168312 * f(i-2) * f(i-2) &
            -9841.58822563111 * f(i-2) * f(i-1) &
            +6406.93231432228 * f(i-2) * f(i+0) &
            -2293.25840018007 * f(i-2) * f(i+1) &
            +347.463467070642 * f(i-2) * f(i+2) &
            +5440.58053610203 * f(i-1) * f(i-1) &
            -7175.23886432754 * f(i-1) * f(i+0) &
            +2601.97484491219 * f(i-1) * f(i+1) &
            -399.198237967023 * f(i-1) * f(i+2) &
            +2403.24289630687 * f(i+0) * f(i+0) &
            -1774.22905245706 * f(i+0) * f(i+1) &
            +277.006890621306 * f(i+0) * f(i+2) &
            +335.04033977354 * f(i+1) * f(i+1) &
            -107.291821204516 * f(i+1) * f(i+2) &
            +8.91870511033141 * f(i+2) * f(i+2)
       sigma(i,6) = +8.91870511033141 * f(i-6) * f(i-6) &
            -135.845449952311 * f(i-6) * f(i-5) &
            +445.58477421919 * f(i-6) * f(i-4) &
            -816.068383469668 * f(i-6) * f(i-3) &
            +901.155248375756 * f(i-6) * f(i-2) &
            -599.696734390095 * f(i-6) * f(i-1) &
            +222.440595557253 * f(i-6) * f(i+0) &
            -35.407460560787 * f(i-6) * f(i+1) &
            +518.201327427861 * f(i-5) * f(i-5) &
            -3406.48778702496 * f(i-5) * f(i-4) &
            +6253.56570910972 * f(i-5) * f(i-3) &
            -6924.03404832822 * f(i-5) * f(i-2) &
            +4621.402767691 * f(i-5) * f(i-1) &
            -1719.62507117959 * f(i-5) * f(i+0) &
            +274.621224828637 * f(i-5) * f(i+1) &
            +5612.02326667724 * f(i-4) * f(i-4) &
            -20664.4461791424 * f(i-4) * f(i-3) &
            +22956.5584871618 * f(i-4) * f(i-2) &
            -15380.2247924521 * f(i-4) * f(i-1) &
            +5746.65947583143 * f(i-4) * f(i+0) &
            -921.690511947415 * f(i-4) * f(i+1) &
            +19089.3249805465 * f(i-3) * f(i-3) &
            -42591.6661577299 * f(i-3) * f(i-2) &
            +28675.0021841139 * f(i-3) * f(i-1) &
            -10772.9570807356 * f(i-3) * f(i+0) &
            +1737.9199467609 * f(i-3) * f(i+1) &
            +23881.8339625102 * f(i-2) * f(i-2) &
            -32362.4054769819 * f(i-2) * f(i-1) &
            +12248.796925352 * f(i-2) * f(i+0) &
            -1992.07290287002 * f(i-2) * f(i+1) &
            +11054.5384359481 * f(i-1) * f(i-1) &
            -8454.36155245706 * f(i-1) * f(i+0) &
            +1391.20673258008 * f(i-1) * f(i+1) &
            +1639.31476541011 * f(i+0) * f(i+0) &
            -549.582823188643 * f(i+0) * f(i+1) &
            +47.5028971986251 * f(i+1) * f(i+1)
       sigma(i,7) = +47.5028971986251 * f(i-7) * f(i-7) &
            -724.638894617214 * f(i-7) * f(i-6) &
            +2385.54101829437 * f(i-7) * f(i-5) &
            -4398.63397429859 * f(i-7) * f(i-4) &
            +4912.48566104661 * f(i-7) * f(i-3) &
            -3328.25158337599 * f(i-7) * f(i-2) &
            +1268.95551054292 * f(i-7) * f(i-1) &
            -210.463531989358 * f(i-7) * f(i+0) &
            +2765.84444133604 * f(i-6) * f(i-6) &
            -18228.7647006052 * f(i-6) * f(i-5) &
            +33651.8387772039 * f(i-6) * f(i-4) &
            -37637.4314325875 * f(i-6) * f(i-3) &
            +25544.3501239796 * f(i-6) * f(i-2) &
            -9759.93192303144 * f(i-6) * f(i-1) &
            +1622.88916698582 * f(i-6) * f(i+0) &
            +30071.0784359481 * f(i-5) * f(i-5) &
            -111189.450476982 * f(i-5) * f(i-4) &
            +124579.67848041 * f(i-5) * f(i-3) &
            -84736.2897924522 * f(i-5) * f(i-2) &
            +32462.762767691 * f(i-5) * f(i-1) &
            -5415.63416825253 * f(i-5) * f(i+0) &
            +102966.44021251 * f(i-4) * f(i-4) &
            -231245.307361434 * f(i-4) * f(i-3) &
            +157723.978487162 * f(i-4) * f(i-2) &
            -60634.3990483282 * f(i-4) * f(i-1) &
            +10159.0931716562 * f(i-4) * f(i+0) &
            +130199.124980546 * f(i-3) * f(i-3) &
            -178245.759975439 * f(i-3) * f(i-2) &
            +68840.1294128134 * f(i-3) * f(i-1) &
            -11602.0447459035 * f(i-3) * f(i+0) &
            +61294.8370166772 * f(i-2) * f(i-2) &
            -47645.872787025 * f(i-2) * f(i-1) &
            +8098.17149379591 * f(i-2) * f(i+0) &
            +9343.02132742786 * f(i-1) * f(i-1) &
            -3217.68658751845 * f(i-1) * f(i+0) &
            +282.837600612977 * f(i+0) * f(i+0)
    end do
  end subroutine smoothness_k8

  subroutine weights_left_k8(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:7)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:7,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7
    do i=8, n-8
       accumulator = 0.0
       omega0 = +0.000155400155400155 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.0087024087024087 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.0913752913752914 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.304584304584305 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega4 = +0.380730380730381 / (1.e-36 + sigma(i,4)) / (1.e-36 + sigma(i,4))
       accumulator = accumulator + omega4
       omega5 = +0.182750582750583 / (1.e-36 + sigma(i,5)) / (1.e-36 + sigma(i,5))
       accumulator = accumulator + omega5
       omega6 = +0.0304584304584305 / (1.e-36 + sigma(i,6)) / (1.e-36 + sigma(i,6))
       accumulator = accumulator + omega6
       omega7 = +0.00124320124320124 / (1.e-36 + sigma(i,7)) / (1.e-36 + sigma(i,7))
       accumulator = accumulator + omega7
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega4 = omega4 / accumulator
       omega5 = omega5 / accumulator
       omega6 = omega6 / accumulator
       omega7 = omega7 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
       omega(i,1,4,1) = omega4
       omega(i,1,5,1) = omega5
       omega(i,1,6,1) = omega6
       omega(i,1,7,1) = omega7
    end do
  end subroutine weights_left_k8
  subroutine reconstruct_left_k8(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:7,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7
    real(kind=8) :: fr0, fr1, fr2, fr3, fr4, fr5, fr6, fr7
    do i=8, n-8
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       omega4 = omega(i,1,4,1)
       omega5 = omega(i,1,5,1)
       omega6 = omega(i,1,6,1)
       omega7 = omega(i,1,7,1)
       fr0 = +2.71785714285714 * f(i+0) &
            -5.28214285714286 * f(i+1) &
            +8.71785714285714 * f(i+2) &
            -9.94880952380952 * f(i+3) &
            +7.55119047619048 * f(i+4) &
            -3.64880952380952 * f(i+5) &
            +1.01785714285714 * f(i+6) &
            -0.125 * f(i+7)
       fr1 = +0.125 * f(i-1) &
            +1.71785714285714 * f(i+0) &
            -1.78214285714286 * f(i+1) &
            +1.71785714285714 * f(i+2) &
            -1.19880952380952 * f(i+3) &
            +0.551190476190476 * f(i+4) &
            -0.148809523809524 * f(i+5) &
            +0.0178571428571429 * f(i+6)
       fr2 = -0.0178571428571429 * f(i-2) &
            +0.267857142857143 * f(i-1) &
            +1.21785714285714 * f(i+0) &
            -0.782142857142857 * f(i+1) &
            +0.467857142857143 * f(i+2) &
            -0.198809523809524 * f(i+3) &
            +0.0511904761904762 * f(i+4) &
            -0.00595238095238095 * f(i+5)
       fr3 = +0.00595238095238095 * f(i-3) &
            -0.0654761904761905 * f(i-2) &
            +0.43452380952381 * f(i-1) &
            +0.884523809523809 * f(i+0) &
            -0.36547619047619 * f(i+1) &
            +0.13452380952381 * f(i+2) &
            -0.0321428571428571 * f(i+3) &
            +0.00357142857142857 * f(i+4)
       fr4 = -0.00357142857142857 * f(i-4) &
            +0.0345238095238095 * f(i-3) &
            -0.16547619047619 * f(i-2) &
            +0.634523809523809 * f(i-1) &
            +0.634523809523809 * f(i+0) &
            -0.16547619047619 * f(i+1) &
            +0.0345238095238095 * f(i+2) &
            -0.00357142857142857 * f(i+3)
       fr5 = +0.00357142857142857 * f(i-5) &
            -0.0321428571428571 * f(i-4) &
            +0.13452380952381 * f(i-3) &
            -0.36547619047619 * f(i-2) &
            +0.884523809523809 * f(i-1) &
            +0.43452380952381 * f(i+0) &
            -0.0654761904761905 * f(i+1) &
            +0.00595238095238095 * f(i+2)
       fr6 = -0.00595238095238095 * f(i-6) &
            +0.0511904761904762 * f(i-5) &
            -0.198809523809524 * f(i-4) &
            +0.467857142857143 * f(i-3) &
            -0.782142857142857 * f(i-2) &
            +1.21785714285714 * f(i-1) &
            +0.267857142857143 * f(i+0) &
            -0.0178571428571429 * f(i+1)
       fr7 = +0.0178571428571429 * f(i-7) &
            -0.148809523809524 * f(i-6) &
            +0.551190476190476 * f(i-5) &
            -1.19880952380952 * f(i-4) &
            +1.71785714285714 * f(i-3) &
            -1.78214285714286 * f(i-2) &
            +1.71785714285714 * f(i-1) &
            +0.125 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3 + &
            fr4 * omega4 + &
            fr5 * omega5 + &
            fr6 * omega6 + &
            fr7 * omega7
    end do
  end subroutine reconstruct_left_k8
  subroutine weights_right_k8(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:7)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:7,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7
    do i=8, n-8
       accumulator = 0.0
       omega0 = +0.00124320124320124 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.0304584304584305 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.182750582750583 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.380730380730381 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega4 = +0.304584304584305 / (1.e-36 + sigma(i,4)) / (1.e-36 + sigma(i,4))
       accumulator = accumulator + omega4
       omega5 = +0.0913752913752914 / (1.e-36 + sigma(i,5)) / (1.e-36 + sigma(i,5))
       accumulator = accumulator + omega5
       omega6 = +0.0087024087024087 / (1.e-36 + sigma(i,6)) / (1.e-36 + sigma(i,6))
       accumulator = accumulator + omega6
       omega7 = +0.000155400155400155 / (1.e-36 + sigma(i,7)) / (1.e-36 + sigma(i,7))
       accumulator = accumulator + omega7
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega4 = omega4 / accumulator
       omega5 = omega5 / accumulator
       omega6 = omega6 / accumulator
       omega7 = omega7 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
       omega(i,1,4,1) = omega4
       omega(i,1,5,1) = omega5
       omega(i,1,6,1) = omega6
       omega(i,1,7,1) = omega7
    end do
  end subroutine weights_right_k8
  subroutine reconstruct_right_k8(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:7,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7
    real(kind=8) :: fr0, fr1, fr2, fr3, fr4, fr5, fr6, fr7
    do i=8, n-8
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       omega4 = omega(i,1,4,1)
       omega5 = omega(i,1,5,1)
       omega6 = omega(i,1,6,1)
       omega7 = omega(i,1,7,1)
       fr0 = +0.125 * f(i+0) &
            +1.71785714285714 * f(i+1) &
            -1.78214285714286 * f(i+2) &
            +1.71785714285714 * f(i+3) &
            -1.19880952380952 * f(i+4) &
            +0.551190476190476 * f(i+5) &
            -0.148809523809524 * f(i+6) &
            +0.0178571428571429 * f(i+7)
       fr1 = -0.0178571428571429 * f(i-1) &
            +0.267857142857143 * f(i+0) &
            +1.21785714285714 * f(i+1) &
            -0.782142857142857 * f(i+2) &
            +0.467857142857143 * f(i+3) &
            -0.198809523809524 * f(i+4) &
            +0.0511904761904762 * f(i+5) &
            -0.00595238095238095 * f(i+6)
       fr2 = +0.00595238095238095 * f(i-2) &
            -0.0654761904761905 * f(i-1) &
            +0.43452380952381 * f(i+0) &
            +0.884523809523809 * f(i+1) &
            -0.36547619047619 * f(i+2) &
            +0.13452380952381 * f(i+3) &
            -0.0321428571428571 * f(i+4) &
            +0.00357142857142857 * f(i+5)
       fr3 = -0.00357142857142857 * f(i-3) &
            +0.0345238095238095 * f(i-2) &
            -0.16547619047619 * f(i-1) &
            +0.634523809523809 * f(i+0) &
            +0.634523809523809 * f(i+1) &
            -0.16547619047619 * f(i+2) &
            +0.0345238095238095 * f(i+3) &
            -0.00357142857142857 * f(i+4)
       fr4 = +0.00357142857142857 * f(i-4) &
            -0.0321428571428571 * f(i-3) &
            +0.13452380952381 * f(i-2) &
            -0.36547619047619 * f(i-1) &
            +0.884523809523809 * f(i+0) &
            +0.43452380952381 * f(i+1) &
            -0.0654761904761905 * f(i+2) &
            +0.00595238095238095 * f(i+3)
       fr5 = -0.00595238095238095 * f(i-5) &
            +0.0511904761904762 * f(i-4) &
            -0.198809523809524 * f(i-3) &
            +0.467857142857143 * f(i-2) &
            -0.782142857142857 * f(i-1) &
            +1.21785714285714 * f(i+0) &
            +0.267857142857143 * f(i+1) &
            -0.0178571428571429 * f(i+2)
       fr6 = +0.0178571428571429 * f(i-6) &
            -0.148809523809524 * f(i-5) &
            +0.551190476190476 * f(i-4) &
            -1.19880952380952 * f(i-3) &
            +1.71785714285714 * f(i-2) &
            -1.78214285714286 * f(i-1) &
            +1.71785714285714 * f(i+0) &
            +0.125 * f(i+1)
       fr7 = -0.125 * f(i-7) &
            +1.01785714285714 * f(i-6) &
            -3.64880952380952 * f(i-5) &
            +7.55119047619048 * f(i-4) &
            -9.94880952380952 * f(i-3) &
            +8.71785714285714 * f(i-2) &
            -5.28214285714286 * f(i-1) &
            +2.71785714285714 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3 + &
            fr4 * omega4 + &
            fr5 * omega5 + &
            fr6 * omega6 + &
            fr7 * omega7
    end do
  end subroutine reconstruct_right_k8
  subroutine smoothness_k9(f, n, sigma)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: sigma(n,0:8)
    integer :: i

    do i=9, n-9
       sigma(i,0) = +669.714981108808 * f(i+0) * f(i+0) &
            -8893.78045641284 * f(i+0) * f(i+1) &
            +26542.9748247279 * f(i+0) * f(i+2) &
            -46202.7564018809 * f(i+0) * f(i+3) &
            +51067.48172473 * f(i+0) * f(i+4) &
            -36585.0576509742 * f(i+0) * f(i+5) &
            +16552.0014823105 * f(i+0) * f(i+6) &
            -4316.36463814691 * f(i+0) * f(i+7) &
            +496.071153428775 * f(i+0) * f(i+8) &
            +29991.6199268498 * f(i+1) * f(i+1) &
            -180813.861301392 * f(i+1) * f(i+2) &
            +316864.946394454 * f(i+1) * f(i+3) &
            -351925.412140196 * f(i+1) * f(i+4) &
            +253037.274363286 * f(i+1) * f(i+5) &
            -114802.39231254 * f(i+1) * f(i+6) &
            +30004.6106074921 * f(i+1) * f(i+7) &
            -3454.6250083913 * f(i+1) * f(i+8) &
            +274437.463974042 * f(i+2) * f(i+2) &
            -966727.207640071 * f(i+2) * f(i+3) &
            +1077735.75959481 * f(i+2) * f(i+4) &
            -777141.668605927 * f(i+2) * f(i+5) &
            +353390.255601705 * f(i+2) * f(i+6) &
            -92530.8427411113 * f(i+2) * f(i+7) &
            +10669.6623191712 * f(i+2) * f(i+8) &
            +854591.512480352 * f(i+3) * f(i+3) &
            -1911046.21916467 * f(i+3) * f(i+4) &
            +1381212.31577771 * f(i+3) * f(i+5) &
            -629240.182223211 * f(i+3) * f(i+6) &
            +165006.507929424 * f(i+3) * f(i+7) &
            -19050.4296324557 * f(i+3) * f(i+8) &
            +1070854.47449811 * f(i+4) * f(i+4) &
            -1550800.19015232 * f(i+4) * f(i+5) &
            +707565.248792344 * f(i+4) * f(i+6) &
            -185776.153674587 * f(i+4) * f(i+7) &
            +21470.5360236718 * f(i+4) * f(i+8) &
            +562311.328233438 * f(i+5) * f(i+5) &
            -513755.446640071 * f(i+5) * f(i+6) &
            +135029.257900627 * f(i+5) * f(i+7) &
            -15619.1414592001 * f(i+5) * f(i+8) &
            +117469.122961696 * f(i+6) * f(i+6) &
            -61801.7019274941 * f(i+6) * f(i+7) &
            +7153.9713035639 * f(i+6) * f(i+8) &
            +8134.55939472587 * f(i+7) * f(i+7) &
            -1884.43224565446 * f(i+7) * f(i+8) &
            +109.193772932945 * f(i+8) * f(i+8)
       sigma(i,1) = +109.193772932945 * f(i-1) * f(i-1) &
            -1469.41675936423 * f(i-1) * f(i+0) &
            +4407.32664278072 * f(i-1) * f(i+1) &
            -7674.89153356352 * f(i-1) * f(i+2) &
            +8466.40114664635 * f(i-1) * f(i+3) &
            -6046.29475543025 * f(i-1) * f(i+4) &
            +2725.41239353462 * f(i-1) * f(i+5) &
            -707.980347608125 * f(i-1) * f(i+6) &
            +81.0556671385418 * f(i-1) * f(i+7) &
            +5049.77020781835 * f(i+0) * f(i+0) &
            -30701.1587180034 * f(i+0) * f(i+1) &
            +53947.0217387824 * f(i+0) * f(i+2) &
            -59895.4013896723 * f(i+0) * f(i+3) &
            +42979.1691915766 * f(i+0) * f(i+4) &
            -19443.7923047687 * f(i+0) * f(i+5) &
            +5065.26308734773 * f(i+0) * f(i+6) &
            -581.225261534807 * f(i+0) * f(i+7) &
            +47140.2493458593 * f(i+1) * f(i+1) &
            -166921.455804809 * f(i+1) * f(i+2) &
            +186372.636616418 * f(i+1) * f(i+3) &
            -134309.272278381 * f(i+1) * f(i+4) &
            +60963.6198256632 * f(i+1) * f(i+5) &
            -15923.1845243454 * f(i+1) * f(i+6) &
            +1830.9895489579 * f(i+1) * f(i+7) &
            +148657.090978519 * f(i+2) * f(i+2) &
            -333527.451742794 * f(i+2) * f(i+3) &
            +241247.066835384 * f(i+2) * f(i+4) &
            -109824.674852455 * f(i+2) * f(i+5) &
            +28752.7613106242 * f(i+2) * f(i+6) &
            -3312.55790820803 * f(i+2) * f(i+7) &
            +187797.717874362 * f(i+3) * f(i+3) &
            -272525.224659461 * f(i+3) * f(i+4) &
            +124378.188236788 * f(i+3) * f(i+5) &
            -32630.2392534289 * f(i+3) * f(i+6) &
            +3765.65529677863 * f(i+3) * f(i+7) &
            +99127.2745988898 * f(i+4) * f(i+4) &
            -90677.1257492533 * f(i+4) * f(i+5) &
            +23831.4757387824 * f(i+4) * f(i+6) &
            -2754.34352099738 * f(i+4) * f(i+7) &
            +20774.7074754889 * f(i+5) * f(i+5) &
            -10936.7033079505 * f(i+5) * f(i+6) &
            +1265.66080746326 * f(i+5) * f(i+7) &
            +1441.28575449258 * f(i+6) * f(i+6) &
            -333.964212406558 * f(i+6) * f(i+7) &
            +19.3647914042221 * f(i+7) * f(i+7)
       sigma(i,2) = +19.3647914042221 * f(i-2) * f(i-2) &
            -267.510578137456 * f(i-2) * f(i-1) &
            +813.039719569186 * f(i-2) * f(i+0) &
            -1422.29540695141 * f(i-2) * f(i+1) &
            +1567.36952565594 * f(i-2) * f(i+2) &
            -1114.27213708535 * f(i-2) * f(i+3) &
            +498.941434911939 * f(i-2) * f(i+4) &
            -128.604173640737 * f(i-2) * f(i+5) &
            +14.6020328694406 * f(i-2) * f(i+6) &
            +948.240872428061 * f(i-1) * f(i-1) &
            -5868.7702184994 * f(i-1) * f(i+0) &
            +10399.3092657059 * f(i-1) * f(i+1) &
            -11568.2032050107 * f(i-1) * f(i+2) &
            +8281.83632095821 * f(i-1) * f(i+3) &
            -3728.0916300002 * f(i-1) * f(i+4) &
            +964.845939313753 * f(i-1) * f(i+5) &
            -109.897639186214 * f(i-1) * f(i+6) &
            +9222.43045243719 * f(i+0) * f(i+0) &
            -33080.8713993305 * f(i+0) * f(i+1) &
            +37137.9417090107 * f(i+0) * f(i+2) &
            -26774.8153713592 * f(i+0) * f(i+3) &
            +12118.1388794826 * f(i+0) * f(i+4) &
            -3149.43314058221 * f(i+0) * f(i+5) &
            +359.908916834338 * f(i+0) * f(i+6) &
            +29975.0953815866 * f(i+1) * f(i+1) &
            -67875.8127912121 * f(i+1) * f(i+2) &
            +49266.8129628913 * f(i+1) * f(i+3) &
            -22417.2306985197 * f(i+1) * f(i+4) &
            +5850.74664900054 * f(i+1) * f(i+5) &
            -670.849344757305 * f(i+1) * f(i+6) &
            +38710.2228777378 * f(i+2) * f(i+2) &
            -56543.4445813355 * f(i+2) * f(i+3) &
            +25858.8233448132 * f(i+2) * f(i+4) &
            -6776.71603569857 * f(i+2) * f(i+5) &
            +779.596278301446 * f(i+2) * f(i+6) &
            +20760.5788136853 * f(i+3) * f(i+3) &
            -19076.8005289601 * f(i+3) * f(i+4) &
            +5018.72343149075 * f(i+3) * f(i+5) &
            -579.197723970813 * f(i+3) * f(i+6) &
            +4400.38698330138 * f(i+4) * f(i+4) &
            -2323.5095791696 * f(i+4) * f(i+5) &
            +268.954810839028 * f(i+4) * f(i+6) &
            +307.688066683541 * f(i+5) * f(i+5) &
            -71.4292240810191 * f(i+5) * f(i+6) &
            +4.15594657554992 * f(i+6) * f(i+6)
       sigma(i,3) = +4.15594657554992 * f(i-3) * f(i-3) &
            -60.2050054904579 * f(i-3) * f(i-2) &
            +189.330514253379 * f(i-3) * f(i-1) &
            -338.290107858048 * f(i-3) * f(i+0) &
            +376.449192281274 * f(i-3) * f(i+1) &
            -267.702258737133 * f(i-3) * f(i+2) &
            +119.001300721573 * f(i-3) * f(i+3) &
            -30.2733426005663 * f(i-3) * f(i+4) &
            +3.3778142788794 * f(i-3) * f(i+5) &
            +224.5781681988 * f(i-2) * f(i-2) &
            -1445.81202311801 * f(i-2) * f(i-1) &
            +2631.07992925861 * f(i-2) * f(i+0) &
            -2970.48199593336 * f(i-2) * f(i+1) &
            +2136.83371274062 * f(i-2) * f(i+2) &
            -958.713082546492 * f(i-2) * f(i+3) &
            +245.728335017173 * f(i-2) * f(i+4) &
            -27.586206325686 * f(i-2) * f(i+5) &
            +2378.03262363703 * f(i-1) * f(i-1) &
            -8815.81240974712 * f(i-1) * f(i+0) &
            +10104.3776503687 * f(i-1) * f(i+1) &
            -7358.38198208451 * f(i-1) * f(i+2) &
            +3334.48145529283 * f(i-1) * f(i+3) &
            -861.616952916863 * f(i-1) * f(i+4) &
            +97.3685006774814 * f(i-1) * f(i+5) &
            +8314.44047543303 * f(i+0) * f(i+0) &
            -19354.0800298309 * f(i+0) * f(i+1) &
            +14276.4079218034 * f(i+0) * f(i+2) &
            -6538.57561788698 * f(i+0) * f(i+3) &
            +1704.37865189406 * f(i+0) * f(i+4) &
            -193.989288499038 * f(i+0) * f(i+5) &
            +11427.8857755966 * f(i+1) * f(i+1) &
            -17079.2799526704 * f(i+1) * f(i+2) &
            +7909.63189419591 * f(i+1) * f(i+3) &
            -2081.09545492807 * f(i+1) * f(i+4) &
            +238.707145323602 * f(i+1) * f(i+5) &
            +6460.89964518612 * f(i+2) * f(i+2) &
            -6051.52109493231 * f(i+2) * f(i+3) &
            +1607.7984936325 * f(i+2) * f(i+4) &
            -185.954130124362 * f(i+2) * f(i+5) &
            +1432.32903721728 * f(i+3) * f(i+3) &
            -768.643244458397 * f(i+3) * f(i+4) &
            +89.6803151793054 * f(i+3) * f(i+5) &
            +104.120555009079 * f(i+4) * f(i+4) &
            -24.5175956580063 * f(i+4) * f(i+5) &
            +1.45672257391222 * f(i+5) * f(i+5)
       sigma(i,4) = +1.45672257391222 * f(i-4) * f(i-4) &
            -22.8431920515406 * f(i-4) * f(i-3) &
            +77.2978189959941 * f(i-4) * f(i-2) &
            -147.360891739772 * f(i-4) * f(i-1) &
            +173.104800126842 * f(i-4) * f(i+0) &
            -128.386943302279 * f(i-4) * f(i+1) &
            +58.7752622928914 * f(i-4) * f(i+2) &
            -15.2037101423747 * f(i-4) * f(i+3) &
            +1.70341067241367 * f(i-4) * f(i+4) &
            +91.7501465525254 * f(i-3) * f(i-3) &
            -634.284062414746 * f(i-3) * f(i-2) &
            +1231.84214048546 * f(i-3) * f(i-1) &
            -1470.62870986083 * f(i-3) * f(i+0) &
            +1106.32708286298 * f(i-3) * f(i+1) &
            -512.943219947286 * f(i-3) * f(i+2) &
            +134.233377963287 * f(i-3) * f(i+3) &
            -15.2037101423747 * f(i-3) * f(i+4) &
            +1119.38719626435 * f(i-2) * f(i-2) &
            -4433.56279439219 * f(i-2) * f(i-1) &
            +5386.9907367885 * f(i-2) * f(i+0) &
            -4116.54995777895 * f(i-2) * f(i+1) &
            +1935.50182392709 * f(i-2) * f(i+2) &
            -512.943219947286 * f(i-2) * f(i+3) &
            +58.7752622928914 * f(i-2) * f(i+4) &
            +4477.71304825325 * f(i-1) * f(i-1) &
            -11088.1845350392 * f(i-1) * f(i+0) &
            +8620.4498023975 * f(i-1) * f(i+1) &
            -4116.54995777895 * f(i-1) * f(i+2) &
            +1106.32708286298 * f(i-1) * f(i+3) &
            -128.386943302279 * f(i-1) * f(i+4) &
            +6998.71770798472 * f(i+0) * f(i+0) &
            -11088.1845350392 * f(i+0) * f(i+1) &
            +5386.9907367885 * f(i+0) * f(i+2) &
            -1470.62870986083 * f(i+0) * f(i+3) &
            +173.104800126842 * f(i+0) * f(i+4) &
            +4477.71304825325 * f(i+1) * f(i+1) &
            -4433.56279439219 * f(i+1) * f(i+2) &
            +1231.84214048546 * f(i+1) * f(i+3) &
            -147.360891739772 * f(i+1) * f(i+4) &
            +1119.38719626435 * f(i+2) * f(i+2) &
            -634.284062414746 * f(i+2) * f(i+3) &
            +77.2978189959941 * f(i+2) * f(i+4) &
            +91.7501465525254 * f(i+3) * f(i+3) &
            -22.8431920515406 * f(i+3) * f(i+4) &
            +1.45672257391222 * f(i+4) * f(i+4)
       sigma(i,5) = +1.45672257391222 * f(i-5) * f(i-5) &
            -24.5175956580063 * f(i-5) * f(i-4) &
            +89.6803151793054 * f(i-5) * f(i-3) &
            -185.954130124362 * f(i-5) * f(i-2) &
            +238.707145323602 * f(i-5) * f(i-1) &
            -193.989288499038 * f(i-5) * f(i+0) &
            +97.3685006774814 * f(i-5) * f(i+1) &
            -27.586206325686 * f(i-5) * f(i+2) &
            +3.3778142788794 * f(i-5) * f(i+3) &
            +104.120555009079 * f(i-4) * f(i-4) &
            -768.643244458397 * f(i-4) * f(i-3) &
            +1607.7984936325 * f(i-4) * f(i-2) &
            -2081.09545492807 * f(i-4) * f(i-1) &
            +1704.37865189406 * f(i-4) * f(i+0) &
            -861.616952916863 * f(i-4) * f(i+1) &
            +245.728335017173 * f(i-4) * f(i+2) &
            -30.2733426005663 * f(i-4) * f(i+3) &
            +1432.32903721728 * f(i-3) * f(i-3) &
            -6051.52109493231 * f(i-3) * f(i-2) &
            +7909.63189419591 * f(i-3) * f(i-1) &
            -6538.57561788698 * f(i-3) * f(i+0) &
            +3334.48145529283 * f(i-3) * f(i+1) &
            -958.713082546492 * f(i-3) * f(i+2) &
            +119.001300721573 * f(i-3) * f(i+3) &
            +6460.89964518612 * f(i-2) * f(i-2) &
            -17079.2799526704 * f(i-2) * f(i-1) &
            +14276.4079218034 * f(i-2) * f(i+0) &
            -7358.38198208451 * f(i-2) * f(i+1) &
            +2136.83371274062 * f(i-2) * f(i+2) &
            -267.702258737133 * f(i-2) * f(i+3) &
            +11427.8857755966 * f(i-1) * f(i-1) &
            -19354.0800298309 * f(i-1) * f(i+0) &
            +10104.3776503687 * f(i-1) * f(i+1) &
            -2970.48199593336 * f(i-1) * f(i+2) &
            +376.449192281274 * f(i-1) * f(i+3) &
            +8314.44047543303 * f(i+0) * f(i+0) &
            -8815.81240974712 * f(i+0) * f(i+1) &
            +2631.07992925861 * f(i+0) * f(i+2) &
            -338.290107858048 * f(i+0) * f(i+3) &
            +2378.03262363703 * f(i+1) * f(i+1) &
            -1445.81202311801 * f(i+1) * f(i+2) &
            +189.330514253379 * f(i+1) * f(i+3) &
            +224.5781681988 * f(i+2) * f(i+2) &
            -60.2050054904579 * f(i+2) * f(i+3) &
            +4.15594657554992 * f(i+3) * f(i+3)
       sigma(i,6) = +4.15594657554992 * f(i-6) * f(i-6) &
            -71.4292240810191 * f(i-6) * f(i-5) &
            +268.954810839028 * f(i-6) * f(i-4) &
            -579.197723970813 * f(i-6) * f(i-3) &
            +779.596278301446 * f(i-6) * f(i-2) &
            -670.849344757305 * f(i-6) * f(i-1) &
            +359.908916834338 * f(i-6) * f(i+0) &
            -109.897639186214 * f(i-6) * f(i+1) &
            +14.6020328694406 * f(i-6) * f(i+2) &
            +307.688066683541 * f(i-5) * f(i-5) &
            -2323.5095791696 * f(i-5) * f(i-4) &
            +5018.72343149075 * f(i-5) * f(i-3) &
            -6776.71603569857 * f(i-5) * f(i-2) &
            +5850.74664900054 * f(i-5) * f(i-1) &
            -3149.43314058221 * f(i-5) * f(i+0) &
            +964.845939313753 * f(i-5) * f(i+1) &
            -128.604173640737 * f(i-5) * f(i+2) &
            +4400.38698330138 * f(i-4) * f(i-4) &
            -19076.8005289601 * f(i-4) * f(i-3) &
            +25858.8233448132 * f(i-4) * f(i-2) &
            -22417.2306985197 * f(i-4) * f(i-1) &
            +12118.1388794826 * f(i-4) * f(i+0) &
            -3728.0916300002 * f(i-4) * f(i+1) &
            +498.941434911939 * f(i-4) * f(i+2) &
            +20760.5788136853 * f(i-3) * f(i-3) &
            -56543.4445813355 * f(i-3) * f(i-2) &
            +49266.8129628913 * f(i-3) * f(i-1) &
            -26774.8153713592 * f(i-3) * f(i+0) &
            +8281.83632095821 * f(i-3) * f(i+1) &
            -1114.27213708535 * f(i-3) * f(i+2) &
            +38710.2228777378 * f(i-2) * f(i-2) &
            -67875.8127912121 * f(i-2) * f(i-1) &
            +37137.9417090107 * f(i-2) * f(i+0) &
            -11568.2032050107 * f(i-2) * f(i+1) &
            +1567.36952565594 * f(i-2) * f(i+2) &
            +29975.0953815866 * f(i-1) * f(i-1) &
            -33080.8713993305 * f(i-1) * f(i+0) &
            +10399.3092657059 * f(i-1) * f(i+1) &
            -1422.29540695141 * f(i-1) * f(i+2) &
            +9222.43045243719 * f(i+0) * f(i+0) &
            -5868.7702184994 * f(i+0) * f(i+1) &
            +813.039719569186 * f(i+0) * f(i+2) &
            +948.240872428061 * f(i+1) * f(i+1) &
            -267.510578137456 * f(i+1) * f(i+2) &
            +19.3647914042221 * f(i+2) * f(i+2)
       sigma(i,7) = +19.3647914042221 * f(i-7) * f(i-7) &
            -333.964212406558 * f(i-7) * f(i-6) &
            +1265.66080746326 * f(i-7) * f(i-5) &
            -2754.34352099738 * f(i-7) * f(i-4) &
            +3765.65529677863 * f(i-7) * f(i-3) &
            -3312.55790820803 * f(i-7) * f(i-2) &
            +1830.9895489579 * f(i-7) * f(i-1) &
            -581.225261534807 * f(i-7) * f(i+0) &
            +81.0556671385418 * f(i-7) * f(i+1) &
            +1441.28575449258 * f(i-6) * f(i-6) &
            -10936.7033079505 * f(i-6) * f(i-5) &
            +23831.4757387824 * f(i-6) * f(i-4) &
            -32630.2392534289 * f(i-6) * f(i-3) &
            +28752.7613106242 * f(i-6) * f(i-2) &
            -15923.1845243454 * f(i-6) * f(i-1) &
            +5065.26308734773 * f(i-6) * f(i+0) &
            -707.980347608125 * f(i-6) * f(i+1) &
            +20774.7074754889 * f(i-5) * f(i-5) &
            -90677.1257492533 * f(i-5) * f(i-4) &
            +124378.188236788 * f(i-5) * f(i-3) &
            -109824.674852455 * f(i-5) * f(i-2) &
            +60963.6198256632 * f(i-5) * f(i-1) &
            -19443.7923047687 * f(i-5) * f(i+0) &
            +2725.41239353462 * f(i-5) * f(i+1) &
            +99127.2745988898 * f(i-4) * f(i-4) &
            -272525.224659461 * f(i-4) * f(i-3) &
            +241247.066835384 * f(i-4) * f(i-2) &
            -134309.272278381 * f(i-4) * f(i-1) &
            +42979.1691915766 * f(i-4) * f(i+0) &
            -6046.29475543025 * f(i-4) * f(i+1) &
            +187797.717874362 * f(i-3) * f(i-3) &
            -333527.451742794 * f(i-3) * f(i-2) &
            +186372.636616418 * f(i-3) * f(i-1) &
            -59895.4013896723 * f(i-3) * f(i+0) &
            +8466.40114664635 * f(i-3) * f(i+1) &
            +148657.090978519 * f(i-2) * f(i-2) &
            -166921.455804809 * f(i-2) * f(i-1) &
            +53947.0217387824 * f(i-2) * f(i+0) &
            -7674.89153356352 * f(i-2) * f(i+1) &
            +47140.2493458593 * f(i-1) * f(i-1) &
            -30701.1587180034 * f(i-1) * f(i+0) &
            +4407.32664278072 * f(i-1) * f(i+1) &
            +5049.77020781835 * f(i+0) * f(i+0) &
            -1469.41675936423 * f(i+0) * f(i+1) &
            +109.193772932945 * f(i+1) * f(i+1)
       sigma(i,8) = +109.193772932945 * f(i-8) * f(i-8) &
            -1884.43224565446 * f(i-8) * f(i-7) &
            +7153.9713035639 * f(i-8) * f(i-6) &
            -15619.1414592001 * f(i-8) * f(i-5) &
            +21470.5360236718 * f(i-8) * f(i-4) &
            -19050.4296324557 * f(i-8) * f(i-3) &
            +10669.6623191712 * f(i-8) * f(i-2) &
            -3454.6250083913 * f(i-8) * f(i-1) &
            +496.071153428775 * f(i-8) * f(i+0) &
            +8134.55939472587 * f(i-7) * f(i-7) &
            -61801.7019274941 * f(i-7) * f(i-6) &
            +135029.257900627 * f(i-7) * f(i-5) &
            -185776.153674587 * f(i-7) * f(i-4) &
            +165006.507929424 * f(i-7) * f(i-3) &
            -92530.8427411113 * f(i-7) * f(i-2) &
            +30004.6106074921 * f(i-7) * f(i-1) &
            -4316.36463814691 * f(i-7) * f(i+0) &
            +117469.122961696 * f(i-6) * f(i-6) &
            -513755.446640071 * f(i-6) * f(i-5) &
            +707565.248792344 * f(i-6) * f(i-4) &
            -629240.182223211 * f(i-6) * f(i-3) &
            +353390.255601705 * f(i-6) * f(i-2) &
            -114802.39231254 * f(i-6) * f(i-1) &
            +16552.0014823105 * f(i-6) * f(i+0) &
            +562311.328233438 * f(i-5) * f(i-5) &
            -1550800.19015232 * f(i-5) * f(i-4) &
            +1381212.31577771 * f(i-5) * f(i-3) &
            -777141.668605927 * f(i-5) * f(i-2) &
            +253037.274363286 * f(i-5) * f(i-1) &
            -36585.0576509742 * f(i-5) * f(i+0) &
            +1070854.47449811 * f(i-4) * f(i-4) &
            -1911046.21916467 * f(i-4) * f(i-3) &
            +1077735.75959481 * f(i-4) * f(i-2) &
            -351925.412140196 * f(i-4) * f(i-1) &
            +51067.48172473 * f(i-4) * f(i+0) &
            +854591.512480352 * f(i-3) * f(i-3) &
            -966727.207640071 * f(i-3) * f(i-2) &
            +316864.946394454 * f(i-3) * f(i-1) &
            -46202.7564018809 * f(i-3) * f(i+0) &
            +274437.463974042 * f(i-2) * f(i-2) &
            -180813.861301392 * f(i-2) * f(i-1) &
            +26542.9748247279 * f(i-2) * f(i+0) &
            +29991.6199268498 * f(i-1) * f(i-1) &
            -8893.78045641284 * f(i-1) * f(i+0) &
            +669.714981108808 * f(i+0) * f(i+0)
    end do
  end subroutine smoothness_k9

  subroutine weights_left_k9(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:8)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:8,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7, omega8
    do i=9, n-9
       accumulator = 0.0
       omega0 = +4.11353352529823e-05 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.00296174413821473 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.0414644179350062 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.193500617030029 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega4 = +0.362813656931304 / (1.e-36 + sigma(i,4)) / (1.e-36 + sigma(i,4))
       accumulator = accumulator + omega4
       omega5 = +0.290250925545043 / (1.e-36 + sigma(i,5)) / (1.e-36 + sigma(i,5))
       accumulator = accumulator + omega5
       omega6 = +0.0967503085150144 / (1.e-36 + sigma(i,6)) / (1.e-36 + sigma(i,6))
       accumulator = accumulator + omega6
       omega7 = +0.0118469765528589 / (1.e-36 + sigma(i,7)) / (1.e-36 + sigma(i,7))
       accumulator = accumulator + omega7
       omega8 = +0.000370218017276841 / (1.e-36 + sigma(i,8)) / (1.e-36 + sigma(i,8))
       accumulator = accumulator + omega8
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega4 = omega4 / accumulator
       omega5 = omega5 / accumulator
       omega6 = omega6 / accumulator
       omega7 = omega7 / accumulator
       omega8 = omega8 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
       omega(i,1,4,1) = omega4
       omega(i,1,5,1) = omega5
       omega(i,1,6,1) = omega6
       omega(i,1,7,1) = omega7
       omega(i,1,8,1) = omega8
    end do
  end subroutine weights_left_k9
  subroutine reconstruct_left_k9(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:8,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7, omega8
    real(kind=8) :: fr0, fr1, fr2, fr3, fr4, fr5, fr6, fr7, fr8
    do i=9, n-9
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       omega4 = omega(i,1,4,1)
       omega5 = omega(i,1,5,1)
       omega6 = omega(i,1,6,1)
       omega7 = omega(i,1,7,1)
       omega8 = omega(i,1,8,1)
       fr0 = +2.82896825396825 * f(i+0) &
            -6.17103174603175 * f(i+1) &
            +11.8289682539683 * f(i+2) &
            -16.1710317460317 * f(i+3) &
            +15.3289682539683 * f(i+4) &
            -9.87103174603175 * f(i+5) &
            +4.12896825396825 * f(i+6) &
            -1.01388888888889 * f(i+7) &
            +0.111111111111111 * f(i+8)
       fr1 = +0.111111111111111 * f(i-1) &
            +1.82896825396825 * f(i+0) &
            -2.17103174603175 * f(i+1) &
            +2.49563492063492 * f(i+2) &
            -2.17103174603175 * f(i+3) &
            +1.32896825396825 * f(i+4) &
            -0.537698412698413 * f(i+5) &
            +0.128968253968254 * f(i+6) &
            -0.0138888888888889 * f(i+7)
       fr2 = -0.0138888888888889 * f(i-2) &
            +0.236111111111111 * f(i-1) &
            +1.32896825396825 * f(i+0) &
            -1.00436507936508 * f(i+1) &
            +0.745634920634921 * f(i+2) &
            -0.421031746031746 * f(i+3) &
            +0.162301587301587 * f(i+4) &
            -0.0376984126984127 * f(i+5) &
            +0.00396825396825397 * f(i+6)
       fr3 = +0.00396825396825397 * f(i-3) &
            -0.0496031746031746 * f(i-2) &
            +0.378968253968254 * f(i-1) &
            +0.995634920634921 * f(i+0) &
            -0.504365079365079 * f(i+1) &
            +0.245634920634921 * f(i+2) &
            -0.0876984126984127 * f(i+3) &
            +0.0194444444444444 * f(i+4) &
            -0.00198412698412698 * f(i+5)
       fr4 = -0.00198412698412698 * f(i-4) &
            +0.0218253968253968 * f(i-3) &
            -0.121031746031746 * f(i-2) &
            +0.545634920634921 * f(i-1) &
            +0.745634920634921 * f(i+0) &
            -0.254365079365079 * f(i+1) &
            +0.078968253968254 * f(i+2) &
            -0.0162698412698413 * f(i+3) &
            +0.00158730158730159 * f(i+4)
       fr5 = +0.00158730158730159 * f(i-5) &
            -0.0162698412698413 * f(i-4) &
            +0.078968253968254 * f(i-3) &
            -0.254365079365079 * f(i-2) &
            +0.745634920634921 * f(i-1) &
            +0.545634920634921 * f(i+0) &
            -0.121031746031746 * f(i+1) &
            +0.0218253968253968 * f(i+2) &
            -0.00198412698412698 * f(i+3)
       fr6 = -0.00198412698412698 * f(i-6) &
            +0.0194444444444444 * f(i-5) &
            -0.0876984126984127 * f(i-4) &
            +0.245634920634921 * f(i-3) &
            -0.504365079365079 * f(i-2) &
            +0.995634920634921 * f(i-1) &
            +0.378968253968254 * f(i+0) &
            -0.0496031746031746 * f(i+1) &
            +0.00396825396825397 * f(i+2)
       fr7 = +0.00396825396825397 * f(i-7) &
            -0.0376984126984127 * f(i-6) &
            +0.162301587301587 * f(i-5) &
            -0.421031746031746 * f(i-4) &
            +0.745634920634921 * f(i-3) &
            -1.00436507936508 * f(i-2) &
            +1.32896825396825 * f(i-1) &
            +0.236111111111111 * f(i+0) &
            -0.0138888888888889 * f(i+1)
       fr8 = -0.0138888888888889 * f(i-8) &
            +0.128968253968254 * f(i-7) &
            -0.537698412698413 * f(i-6) &
            +1.32896825396825 * f(i-5) &
            -2.17103174603175 * f(i-4) &
            +2.49563492063492 * f(i-3) &
            -2.17103174603175 * f(i-2) &
            +1.82896825396825 * f(i-1) &
            +0.111111111111111 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3 + &
            fr4 * omega4 + &
            fr5 * omega5 + &
            fr6 * omega6 + &
            fr7 * omega7 + &
            fr8 * omega8
    end do
  end subroutine reconstruct_left_k9
  subroutine weights_right_k9(sigma, n, omega)

    implicit none

    real(kind=8), intent(in) :: sigma(n,0:8)
    integer, intent(in) :: n
    real(kind=8), intent(out) :: omega(n,1,0:8,2)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7, omega8
    do i=9, n-9
       accumulator = 0.0
       omega0 = +0.000370218017276841 / (1.e-36 + sigma(i,0)) / (1.e-36 + sigma(i,0))
       accumulator = accumulator + omega0
       omega1 = +0.0118469765528589 / (1.e-36 + sigma(i,1)) / (1.e-36 + sigma(i,1))
       accumulator = accumulator + omega1
       omega2 = +0.0967503085150144 / (1.e-36 + sigma(i,2)) / (1.e-36 + sigma(i,2))
       accumulator = accumulator + omega2
       omega3 = +0.290250925545043 / (1.e-36 + sigma(i,3)) / (1.e-36 + sigma(i,3))
       accumulator = accumulator + omega3
       omega4 = +0.362813656931304 / (1.e-36 + sigma(i,4)) / (1.e-36 + sigma(i,4))
       accumulator = accumulator + omega4
       omega5 = +0.193500617030029 / (1.e-36 + sigma(i,5)) / (1.e-36 + sigma(i,5))
       accumulator = accumulator + omega5
       omega6 = +0.0414644179350062 / (1.e-36 + sigma(i,6)) / (1.e-36 + sigma(i,6))
       accumulator = accumulator + omega6
       omega7 = +0.00296174413821473 / (1.e-36 + sigma(i,7)) / (1.e-36 + sigma(i,7))
       accumulator = accumulator + omega7
       omega8 = +4.11353352529823e-05 / (1.e-36 + sigma(i,8)) / (1.e-36 + sigma(i,8))
       accumulator = accumulator + omega8
       omega0 = omega0 / accumulator
       omega1 = omega1 / accumulator
       omega2 = omega2 / accumulator
       omega3 = omega3 / accumulator
       omega4 = omega4 / accumulator
       omega5 = omega5 / accumulator
       omega6 = omega6 / accumulator
       omega7 = omega7 / accumulator
       omega8 = omega8 / accumulator
       omega(i,1,0,1) = omega0
       omega(i,1,1,1) = omega1
       omega(i,1,2,1) = omega2
       omega(i,1,3,1) = omega3
       omega(i,1,4,1) = omega4
       omega(i,1,5,1) = omega5
       omega(i,1,6,1) = omega6
       omega(i,1,7,1) = omega7
       omega(i,1,8,1) = omega8
    end do
  end subroutine weights_right_k9
  subroutine reconstruct_right_k9(f, n, omega, fr)

    implicit none

    real(kind=8), intent(in) :: f(n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: omega(n,1,0:8,2)
    real(kind=8), intent(out) :: fr(n,1)
    integer :: i
    real(kind=8) :: accumulator

    real(kind=8) :: omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7, omega8
    real(kind=8) :: fr0, fr1, fr2, fr3, fr4, fr5, fr6, fr7, fr8
    do i=9, n-9
       omega0 = omega(i,1,0,1)
       omega1 = omega(i,1,1,1)
       omega2 = omega(i,1,2,1)
       omega3 = omega(i,1,3,1)
       omega4 = omega(i,1,4,1)
       omega5 = omega(i,1,5,1)
       omega6 = omega(i,1,6,1)
       omega7 = omega(i,1,7,1)
       omega8 = omega(i,1,8,1)
       fr0 = +0.111111111111111 * f(i+0) &
            +1.82896825396825 * f(i+1) &
            -2.17103174603175 * f(i+2) &
            +2.49563492063492 * f(i+3) &
            -2.17103174603175 * f(i+4) &
            +1.32896825396825 * f(i+5) &
            -0.537698412698413 * f(i+6) &
            +0.128968253968254 * f(i+7) &
            -0.0138888888888889 * f(i+8)
       fr1 = -0.0138888888888889 * f(i-1) &
            +0.236111111111111 * f(i+0) &
            +1.32896825396825 * f(i+1) &
            -1.00436507936508 * f(i+2) &
            +0.745634920634921 * f(i+3) &
            -0.421031746031746 * f(i+4) &
            +0.162301587301587 * f(i+5) &
            -0.0376984126984127 * f(i+6) &
            +0.00396825396825397 * f(i+7)
       fr2 = +0.00396825396825397 * f(i-2) &
            -0.0496031746031746 * f(i-1) &
            +0.378968253968254 * f(i+0) &
            +0.995634920634921 * f(i+1) &
            -0.504365079365079 * f(i+2) &
            +0.245634920634921 * f(i+3) &
            -0.0876984126984127 * f(i+4) &
            +0.0194444444444444 * f(i+5) &
            -0.00198412698412698 * f(i+6)
       fr3 = -0.00198412698412698 * f(i-3) &
            +0.0218253968253968 * f(i-2) &
            -0.121031746031746 * f(i-1) &
            +0.545634920634921 * f(i+0) &
            +0.745634920634921 * f(i+1) &
            -0.254365079365079 * f(i+2) &
            +0.078968253968254 * f(i+3) &
            -0.0162698412698413 * f(i+4) &
            +0.00158730158730159 * f(i+5)
       fr4 = +0.00158730158730159 * f(i-4) &
            -0.0162698412698413 * f(i-3) &
            +0.078968253968254 * f(i-2) &
            -0.254365079365079 * f(i-1) &
            +0.745634920634921 * f(i+0) &
            +0.545634920634921 * f(i+1) &
            -0.121031746031746 * f(i+2) &
            +0.0218253968253968 * f(i+3) &
            -0.00198412698412698 * f(i+4)
       fr5 = -0.00198412698412698 * f(i-5) &
            +0.0194444444444444 * f(i-4) &
            -0.0876984126984127 * f(i-3) &
            +0.245634920634921 * f(i-2) &
            -0.504365079365079 * f(i-1) &
            +0.995634920634921 * f(i+0) &
            +0.378968253968254 * f(i+1) &
            -0.0496031746031746 * f(i+2) &
            +0.00396825396825397 * f(i+3)
       fr6 = +0.00396825396825397 * f(i-6) &
            -0.0376984126984127 * f(i-5) &
            +0.162301587301587 * f(i-4) &
            -0.421031746031746 * f(i-3) &
            +0.745634920634921 * f(i-2) &
            -1.00436507936508 * f(i-1) &
            +1.32896825396825 * f(i+0) &
            +0.236111111111111 * f(i+1) &
            -0.0138888888888889 * f(i+2)
       fr7 = -0.0138888888888889 * f(i-7) &
            +0.128968253968254 * f(i-6) &
            -0.537698412698413 * f(i-5) &
            +1.32896825396825 * f(i-4) &
            -2.17103174603175 * f(i-3) &
            +2.49563492063492 * f(i-2) &
            -2.17103174603175 * f(i-1) &
            +1.82896825396825 * f(i+0) &
            +0.111111111111111 * f(i+1)
       fr8 = +0.111111111111111 * f(i-8) &
            -1.01388888888889 * f(i-7) &
            +4.12896825396825 * f(i-6) &
            -9.87103174603175 * f(i-5) &
            +15.3289682539683 * f(i-4) &
            -16.1710317460317 * f(i-3) &
            +11.8289682539683 * f(i-2) &
            -6.17103174603175 * f(i-1) &
            +2.82896825396825 * f(i+0)
       fr(i,1) = fr0 * omega0 + &
            fr1 * omega1 + &
            fr2 * omega2 + &
            fr3 * omega3 + &
            fr4 * omega4 + &
            fr5 * omega5 + &
            fr6 * omega6 + &
            fr7 * omega7 + &
            fr8 * omega8
    end do
  end subroutine reconstruct_right_k9

end module reconstruct
