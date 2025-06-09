subroutine derivadas_x(nx, ny, hx, p, eta, dxp)
  implicit none
  integer, intent(in)    :: nx, ny
  real(kind=8), intent(in)  :: hx
  real(kind=8), intent(in)  :: p(nx, ny), eta(nx, ny)
  real(kind=8), intent(out) :: dxp(nx, ny)

  ! WENO parameters
  real(kind=8), parameter :: gamma0 = 1.0d0/6.0d0, gamma1 = 4.0d0/6.0d0, gamma2 = 1.0d0/6.0d0
  real(kind=8), parameter :: epsilon = 1.0d0

  real(kind=8) :: beta0, beta1, beta2, alpha0, alpha1, alpha2, sum_alpha
  real(kind=8) :: w0, w1, w2, D0, D1, D2
  integer :: i, j, b

    !-- Interior: 5th-order C-WENO derivatives in x --
    do i = 3, nx-2
      do j = 1, ny
        ! Candidate stencils
        D0 = (  p(i-2,j) - 4.d0*p(i-1,j) + 3.d0*p(i  ,j) ) / (2.d0*hx)
        D1 = (              p(i+1,j) -         p(i-1,j) ) / (2.d0*hx)
        D2 = ( -3.d0*p(i  ,j) + 4.d0*p(i+1,j) -    p(i+2,j) ) / (2.d0*hx)
        ! Smoothness indicators
        beta0 = (13.d0/12.d0)*(p(i-2,j)-2.d0*p(i-1,j)+p(i,j))**2 + &
                (1.d0/4.d0)*(p(i-2,j)-4.d0*p(i-1,j)+3.d0*p(i,j))**2
        beta1 = (13.d0/12.d0)*(p(i-1,j)-2.d0*p(i  ,j)+p(i+1,j))**2 + &
                (1.d0/4.d0)*(p(i+1,j)-p(i-1,j))**2
        beta2 = (13.d0/12.d0)*(p(i  ,j)-2.d0*p(i+1,j)+p(i+2,j))**2 + &
                (1.d0/4.d0)*(3.d0*p(i  ,j)-4.d0*p(i+1,j)+p(i+2,j))**2
        ! Nonlinear weights
        alpha0 = gamma0 / ( (epsilon+beta0)**2 )
        alpha1 = gamma1 / ( (epsilon+beta1)**2 )
        alpha2 = gamma2 / ( (epsilon+beta2)**2 )
        sum_alpha = alpha0 + alpha1 + alpha2
        w0 = alpha0 / sum_alpha;  w1 = alpha1 / sum_alpha;  w2 = alpha2 / sum_alpha
        ! Final derivative
        dxp(i,j) = w0*D0 + w1*D1 + w2*D2
      end do
    end do


b=0

    !-- Boundaries: depending on b flag --
    select case (b)
    case (0)  ! 2nd-order one-sided
      do j = 1, ny
        dxp(1  ,j) = (-3.d0*p(1,j) + 4.d0*p(2,j) -   p(3,j)) / (2.d0*hx)
        dxp(2  ,j) = (-3.d0*p(2,j) + 4.d0*p(3,j) -   p(4,j)) / (2.d0*hx)
        dxp(nx-1,j) = ( 3.d0*p(nx-1,j) - 4.d0*p(nx-2,j) + p(nx-3,j)) / (2.d0*hx)
        dxp(nx  ,j) = ( 3.d0*p(nx  ,j) - 4.d0*p(nx-1,j) + p(nx-2,j)) / (2.d0*hx)
      end do
    case (1)  ! periodic
      do j = 1, ny
        dxp(1  ,j) = (p(2,j) - p(nx,j))    / (2.d0*hx)
        dxp(2  ,j) = (p(3,j) - p(1 ,j))    / (2.d0*hx)
        dxp(nx-1,j) = (p(1,j) - p(nx-2,j)) / (2.d0*hx)
        dxp(nx  ,j) = (p(2,j) - p(nx-1,j)) / (2.d0*hx)
      end do
    case (2)  ! outflow (one-sided)
      do j = 1, ny
        dxp(1  ,j) = (-3.d0*p(1,j) + 4.d0*p(2,j) -   p(3,j)) / (2.d0*hx)
        dxp(2  ,j) = (-3.d0*p(2,j) + 4.d0*p(3,j) -   p(4,j)) / (2.d0*hx)
        dxp(nx-1,j) = ( 3.d0*p(nx-1,j) - 4.d0*p(nx-2,j) + p(nx-3,j)) / (2.d0*hx)
        dxp(nx  ,j) = ( 3.d0*p(nx  ,j) - 4.d0*p(nx-1,j) + p(nx-2,j)) / (2.d0*hx)
      end do
    case default
      do j = 1, ny
        dxp(:,j) = 0.d0
      end do
    end select


end subroutine derivadas_x
