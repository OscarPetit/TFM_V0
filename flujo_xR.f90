subroutine flujos_xR(nx, ny, psi_lx, psi_rx, ux_lx, ux_rx, uy_lx, uy_rx, ut_lx, ut_rx, A_lx, A_rx, Qx_lx, Qx_rx, Qy_lx, Qy_rx, Qt_lx, Qt_rx, sigmatx_lx, sigmatx_rx, sigmaxx_lx, sigmaxx_rx, sigmaxy_lx, sigmaxy_rx, dens_lx, dens_rx, fxt_lxR, fxt_rxR, fxx_lxR, fxx_rxR, fxy_lxR, fxy_rxR, f0xt_lxR, f0xt_rxR, f0xx_lxR, f0xx_rxR, f0xy_lxR, f0xy_rxR, fdensx_lxR, fdensx_rxR)
  implicit none
  
  !----------------------------------------------------------------------
  !  Input parameters
  !----------------------------------------------------------------------
  integer, intent(in) :: nx, ny
  
  real*8, intent(in) :: psi_lx(nx, ny), psi_rx(nx, ny)
  real*8, intent(in) :: ux_lx(nx, ny),  ux_rx(nx, ny)
  real*8, intent(in) :: uy_lx(nx, ny),  uy_rx(nx, ny)
  real*8, intent(in) :: ut_lx(nx, ny),  ut_rx(nx, ny)

  
  real*8, intent(in) :: A_lx(nx, ny),   A_rx(nx, ny)
  real*8, intent(in) :: Qx_lx(nx, ny),  Qx_rx(nx, ny)
  real*8, intent(in) :: Qy_lx(nx, ny),  Qy_rx(nx, ny)
  real*8, intent(in) :: Qt_lx(nx, ny),  Qt_rx(nx, ny)
  
  real*8, intent(in) :: sigmatx_lx(nx, ny), sigmatx_rx(nx, ny)
  real*8, intent(in) :: sigmaxx_lx(nx, ny), sigmaxx_rx(nx, ny)
  real*8, intent(in) :: sigmaxy_lx(nx, ny), sigmaxy_rx(nx, ny)
  
  real*8, intent(in) :: dens_lx(nx, ny), dens_rx(nx, ny)
  
  !----------------------------------------------------------------------
  !  Output parameters
  !----------------------------------------------------------------------
  real*8, intent(out) :: fxt_lxR(nx, ny), fxt_rxR(nx, ny)
  real*8, intent(out) :: fxx_lxR(nx, ny), fxx_rxR(nx, ny)
  real*8, intent(out) :: fxy_lxR(nx, ny), fxy_rxR(nx, ny)
  real*8, intent(out) :: f0xt_lxR(nx, ny), f0xt_rxR(nx, ny)
  real*8, intent(out) :: f0xx_lxR(nx, ny), f0xx_rxR(nx, ny)
  real*8, intent(out) :: f0xy_lxR(nx, ny), f0xy_rxR(nx, ny)
  real*8, intent(out) :: fdensx_lxR(nx, ny), fdensx_rxR(nx, ny)
    integer :: i, j
  
  
  
  
    do j = 3, ny-3
      do i = 3, nx-3
        
        fxt_rxR(i,j) = exp(psi_rx(i,j))*(4.0/3.0*ux_rx(i,j)*ut_rx(i,j)) + A_rx(i,j)*(4.0/3.0*ux_rx(i,j)*ut_rx(i,j)) + ut_rx(i,j)*Qx_rx(i,j) + ux_rx(i,j)*Qt_rx(i,j) + sigmatx_rx(i,j)
        fxt_lxR(i,j) = exp(psi_lx(i,j))*(4.0/3.0*ux_lx(i,j)*ut_lx(i,j)) + A_lx(i,j)*(4.0/3.0*ux_lx(i,j)*ut_lx(i,j)) + ut_lx(i,j)*Qx_lx(i,j) + ux_lx(i,j)*Qt_lx(i,j) + sigmatx_lx(i,j)
  
        fxx_rxR(i,j) = (exp(psi_rx(i,j)) + A_rx(i,j))*(4.0/3.0*ux_rx(i,j)**2 + 1.0/3.0) + 2.0*ux_rx(i,j)*Qx_rx(i,j) + sigmaxx_rx(i,j)
        fxx_lxR(i,j) = (exp(psi_lx(i,j)) + A_lx(i,j))*(4.0/3.0*ux_lx(i,j)**2 + 1.0/3.0) + 2.0*ux_lx(i,j)*Qx_lx(i,j) + sigmaxx_lx(i,j)
  
        fxy_rxR(i,j) = (exp(psi_rx(i,j)) + A_rx(i,j))*(4.0/3.0*ux_rx(i,j)*uy_rx(i,j)) + ux_rx(i,j)*Qy_rx(i,j) + uy_rx(i,j)*Qx_rx(i,j) + sigmaxy_rx(i,j)
        fxy_lxR(i,j) = (exp(psi_lx(i,j)) + A_lx(i,j))*(4.0/3.0*ux_lx(i,j)*uy_lx(i,j)) + ux_lx(i,j)*Qy_lx(i,j) + uy_lx(i,j)*Qx_lx(i,j) + sigmaxy_lx(i,j)
  
        fdensx_rxR(i,j) = dens_rx(i,j)*ux_rx(i,j)
        fdensx_lxR(i,j) = dens_lx(i,j)*ux_lx(i,j)
  
        f0xt_rxR(i,j) = exp(psi_rx(i,j))*(4.0/3.0*ux_rx(i,j)*ut_rx(i,j))
        f0xt_lxR(i,j) = exp(psi_lx(i,j))*(4.0/3.0*ux_lx(i,j)*ut_lx(i,j))
  
        f0xx_rxR(i,j) = exp(psi_rx(i,j))*(4.0/3.0*ux_rx(i,j)*ux_rx(i,j) + 1.0/3.0)
        f0xx_lxR(i,j) = exp(psi_lx(i,j))*(4.0/3.0*ux_lx(i,j)*ux_lx(i,j) + 1.0/3.0)
    
        f0xy_rxR(i,j) = exp(psi_rx(i,j))*(4.0/3.0*ux_rx(i,j)*uy_rx(i,j))
        f0xy_lxR(i,j) = exp(psi_lx(i,j))*(4.0/3.0*ux_lx(i,j)*uy_lx(i,j))
  
    
  
      enddo
    enddo
  
  
  
  
  
  
  
  
  
  end subroutine flujos_xR
  