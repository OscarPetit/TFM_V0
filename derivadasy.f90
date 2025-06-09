subroutine derivadas_y(nx, ny, hy, p, eta, dyp)
  implicit none
  integer, intent(in)    :: nx, ny
  real(kind=8), intent(in)  :: hy
  real(kind=8), intent(in)  :: p(nx, ny), eta(nx, ny)
  real(kind=8), intent(out) :: dyp(nx, ny)

  ! WENO parameters
  real(kind=8), parameter :: gamma0 = 1.0d0/6.0d0, gamma1 = 4.0d0/6.0d0, gamma2 = 1.0d0/6.0d0
  real(kind=8), parameter :: epsilon = 1.0d0

  real(kind=8) :: beta0, beta1, beta2, alpha0, alpha1, alpha2, sum_alpha
  real(kind=8) :: w0, w1, w2, D0, D1, D2
  integer :: i, j, b

    !-- Interior: 5th-order C-WENO derivatives in y --
    do i = 1, nx
      do j = 3, ny-2
        ! Candidate stencils
        D0 = (  p(i,j-2) - 4.d0*p(i,j-1) + 3.d0*p(i,j  )) / (2.d0*hy)
        D1 = (              p(i,j+1) -         p(i,j-1)) / (2.d0*hy)
        D2 = ( -3.d0*p(i,j  ) + 4.d0*p(i,j+1) -    p(i,j+2)) / (2.d0*hy)
        ! Smoothness indicators
        beta0 = (13.d0/12.d0)*(p(i,j-2)-2.d0*p(i,j-1)+p(i,j))**2 + &
                (1.d0/4.d0)*(p(i,j-2)-4.d0*p(i,j-1)+3.d0*p(i,j))**2
        beta1 = (13.d0/12.d0)*(p(i,j-1)-2.d0*p(i,j  )+p(i,j+1))**2 + &
                (1.d0/4.d0)*(p(i,j+1)-p(i,j-1))**2
        beta2 = (13.d0/12.d0)*(p(i,j  )-2.d0*p(i,j+1)+p(i,j+2))**2 + &
                (1.d0/4.d0)*(3.d0*p(i,j  )-4.d0*p(i,j+1)+p(i,j+2))**2
        ! Nonlinear weights
        alpha0 = gamma0 / ( (epsilon+beta0)**4 )
        alpha1 = gamma1 / ( (epsilon+beta1)**4 )
        alpha2 = gamma2 / ( (epsilon+beta2)**4 )
        sum_alpha = alpha0 + alpha1 + alpha2
        w0 = alpha0 / sum_alpha;  w1 = alpha1 / sum_alpha;  w2 = alpha2 / sum_alpha
        ! Final derivative
        dyp(i,j) = w0*D0 + w1*D1 + w2*D2
      end do
    end do


    b=0
    !-- Boundaries: depending on b flag --
    select case (b)
    case (0)  ! 2nd-order one-sided
      do i = 1, nx
        dyp(i,1  ) = (-3.d0*p(i,1) + 4.d0*p(i,2) -   p(i,3)) / (2.d0*hy)
        dyp(i,2  ) = (-3.d0*p(i,2) + 4.d0*p(i,3) -   p(i,4)) / (2.d0*hy)
        dyp(i,ny-1) = ( 3.d0*p(i,ny-1) - 4.d0*p(i,ny-2) + p(i,ny-3)) / (2.d0*hy)
        dyp(i,ny  ) = ( 3.d0*p(i,ny  ) - 4.d0*p(i,ny-1) + p(i,ny-2)) / (2.d0*hy)
      end do
    case (1)  ! periodic
      do i = 1, nx
        dyp(i,1  ) = (p(i,2) - p(i,ny))    / (2.d0*hy)
        dyp(i,2  ) = (p(i,3) - p(i,1 ))    / (2.d0*hy)
        dyp(i,ny-1) = (p(i,ny) - p(i,ny-2)) / (2.d0*hy)
        dyp(i,ny  ) = (p(i,1)  - p(i,ny-1)) / (2.d0*hy)
      end do
    case (2)  ! outflow (one-sided)
      do i = 1, nx
        dyp(i,1  ) = (-3.d0*p(i,1) + 4.d0*p(i,2) -   p(i,3)) / (2.d0*hy)
        dyp(i,2  ) = (-3.d0*p(i,2) + 4.d0*p(i,3) -   p(i,4)) / (2.d0*hy)
        dyp(i,ny-1) = ( 3.d0*p(i,ny-1) - 4.d0*p(i,ny-2) + p(i,ny-3)) / (2.d0*hy)
        dyp(i,ny  ) = ( 3.d0*p(i,ny  ) - 4.d0*p(i,ny-1) + p(i,ny-2)) / (2.d0*hy)
      end do
    case default
      do i = 1, nx
        dyp(i,:) = 0.d0
      end do
    end select


end subroutine derivadas_y
