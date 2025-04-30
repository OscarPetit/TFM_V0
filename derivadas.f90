subroutine derivadas(nx, ny, hx, hy, p, dxp, dyp)
  implicit none
  integer, intent(in) :: nx, ny           ! Number of grid points in x and y directions
  real*8, intent(in) :: hx, hy             ! Grid spacings in x and y directions
  real*8, intent(in) :: p(nx, ny)          ! Array of p values
  real*8, intent(out) :: dxp(nx, ny), dyp(nx, ny) ! Arrays of dp/dx and dp/dy values

  real*8 :: w0, w1, w2                      ! Weights
  real*8 :: beta0, beta1, beta2             ! Smoothness indicators
  real*8 :: gamma0, gamma1, gamma2          ! Linear weights
  real*8 :: alpha0, alpha1, alpha2          ! Nonlinear weights
  real*8 :: sum_alpha                       ! Sum of nonlinear weights
  real*8 :: epsilon                         ! Small parameter to avoid division by zero

  integer :: i, j



  ! Initialize linear weights (example values; adjust if needed) 
  gamma0 = 1.0/6.0
  gamma1 = 2.0/3.0
  gamma2 = 1.0/6.0

  ! Small parameter to prevent division by zero (MIRAR FIGURA 3 BNDK PAPER)
  epsilon =   10.0e-2

  ! Loop over interior points to compute weights and derivatives
  do j = 3, ny-2
    do i = 3, nx-2
      ! Calculate smoothness indicators (beta0, beta1, beta2) for x-direction
      beta2 = (1.0/4.0) * (3.0 * p(i, j) - 4.0 * p(i+1, j) + p(i+2, j))**2 + &
              (13.0/12.0) * (p(i, j) - 2.0 * p(i+1, j) + p(i+2, j))**2
      beta1 = (1.0/4.0) * (p(i+1, j) - p(i-1, j))**2 + &
              (13.0/12.0) * (p(i-1, j) - 2.0 * p(i, j) + p(i+1, j))**2
      beta0 = (1.0/4.0) * (p(i-2, j) - 4.0 * p(i-1, j) + 3.0 * p(i, j))**2 + &
              (13.0/12.0) * (p(i-2, j) - 2.0 * p(i-1, j) + p(i, j))**2

      ! Calculate nonlinear weights (alpha0, alpha1, alpha2)
      alpha0 = gamma0 / (epsilon + beta0)**2
      alpha1 = gamma1 / (epsilon + beta1)**2
      alpha2 = gamma2 / (epsilon + beta2)**2

      ! Normalize weights
      sum_alpha = alpha0 + alpha1 + alpha2
      w0 = alpha0 / sum_alpha
      w1 = alpha1 / sum_alpha
      w2 = alpha2 / sum_alpha

      ! Calculate dxp (x-derivative)
      dxp(i, j) = ((p(i-2, j) - 4.0 * p(i-1, j) + 3.0 * p(i, j)) / (2.0 * hx)) * w0 + &
                  ((p(i+1, j) - p(i-1, j)) / (2.0 * hx)) * w1 + &
                  ((-3.0 * p(i, j) + 4.0 * p(i+1, j) - p(i+2, j)) / (2.0 * hx)) * w2


        enddo
enddo

do i = 3, nx-2
        do j = 3, ny-2
      ! Calculate smoothness indicators (beta0, beta1, beta2) for y-direction
      beta2 = (1.0/4.0) * (3.0 * p(i, j) - 4.0 * p(i, j+1) + p(i, j+2))**2 + &
              (13.0/12.0) * (p(i, j) - 2.0 * p(i, j+1) + p(i, j+2))**2
      beta1 = (1.0/4.0) * (p(i, j+1) - p(i, j-1))**2 + &
              (13.0/12.0) * (p(i, j-1) - 2.0 * p(i, j) + p(i, j+1))**2
      beta0 = (1.0/4.0) * (p(i, j-2) - 4.0 * p(i, j-1) + 3.0 * p(i, j))**2 + &
              (13.0/12.0) * (p(i, j-2) - 2.0 * p(i, j-1) + p(i, j))**2

      ! Calculate nonlinear weights (alpha0, alpha1, alpha2) for y-direction
      alpha0 = gamma0 / (epsilon + beta0)**2
      alpha1 = gamma1 / (epsilon + beta1)**2
      alpha2 = gamma2 / (epsilon + beta2)**2

      ! Normalize weights
      sum_alpha = alpha0 + alpha1 + alpha2
      w0 = alpha0 / sum_alpha
      w1 = alpha1 / sum_alpha
      w2 = alpha2 / sum_alpha

      ! Calculate dyp (y-derivative)
      dyp(i, j) = ((p(i, j-2) - 4.0 * p(i, j-1) + 3.0 * p(i, j)) / (2.0 * hy)) * w0 + &
                  ((p(i, j+1) - p(i, j-1)) / (2.0 * hy)) * w1 + &
                  ((-3.0 * p(i, j) + 4.0 * p(i, j+1) - p(i, j+2)) / (2.0 * hy)) * w2
    end do
  end do





end subroutine derivadas
