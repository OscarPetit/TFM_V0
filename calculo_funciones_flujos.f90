subroutine calculo_funciones_flujos(nx, ny, psi_lx, psi_rx, psi_ly, psi_ry, ux_lx, ux_rx, ux_ly, ux_ry, uy_lx, uy_rx, uy_ly, uy_ry, ut_lx, ut_rx, ut_ly, ut_ry, dtpsi_lx, dtpsi_rx, dtpsi_ly, dtpsi_ry, dtux_lx, dtux_rx, dtux_ly, dtux_ry, dtuy_lx, dtuy_rx, dtuy_ly, dtuy_ry,dxpsi_lx, dxpsi_rx, dxpsi_ly, dxpsi_ry, dxux_lx, dxux_rx, dxux_ly, dxux_ry, dxuy_lx, dxuy_rx, dxuy_ly, dxuy_ry, dypsi_lx, dypsi_rx, dypsi_ly, dypsi_ry, dyux_lx, dyux_rx, dyux_ly, dyux_ry, dyuy_lx, dyuy_rx, dyuy_ly, dyuy_ry, eta, A_lx, A_rx, A_ly, A_ry, Qx_lx, Qx_rx, Qx_ly, Qx_ry, Qy_lx, Qy_rx, Qy_ly, Qy_ry, Qt_lx, Qt_rx, Qt_ly, Qt_ry, sigmaxx_lx, sigmaxx_rx, sigmaxx_ly, sigmaxx_ry, sigmaxy_lx, sigmaxy_rx, sigmaxy_ly, sigmaxy_ry, sigmayy_lx, sigmayy_rx, sigmayy_ly, sigmayy_ry, sigmatx_lx, sigmatx_rx, sigmatx_ly, sigmatx_ry, sigmaty_lx, sigmaty_rx, sigmaty_ly, sigmaty_ry)
  implicit none
  integer, intent(in) :: nx, ny           ! Number of grid points in x and y directions
  real(kind=8), intent(in) :: dtpsi_lx(nx, ny), dtpsi_rx(nx, ny), dtpsi_ly(nx, ny), dtpsi_ry(nx, ny), dtux_lx(nx, ny), dtux_rx(nx, ny), dtux_ly(nx, ny), dtux_ry(nx, ny), dtuy_lx(nx, ny), dtuy_rx(nx, ny), dtuy_ly(nx, ny), dtuy_ry(nx, ny), dxpsi_lx(nx, ny), dxpsi_rx(nx, ny), dxpsi_ly(nx, ny), dxpsi_ry(nx, ny), dxux_lx(nx, ny), dxux_rx(nx, ny), dxux_ly(nx, ny), dxux_ry(nx, ny), dxuy_lx(nx, ny), dxuy_rx(nx, ny), dxuy_ly(nx, ny), dxuy_ry(nx, ny), dypsi_lx(nx, ny), dypsi_rx(nx, ny), dypsi_ly(nx, ny), dypsi_ry(nx, ny), dyux_lx(nx, ny), dyux_rx(nx, ny), dyux_ly(nx, ny), dyux_ry(nx, ny), dyuy_lx(nx, ny), dyuy_rx(nx, ny), dyuy_ly(nx, ny), dyuy_ry(nx, ny), psi_lx(nx, ny), psi_rx(nx, ny), psi_ly(nx, ny), psi_ry(nx, ny), ux_lx(nx, ny), ux_rx(nx, ny), ux_ly(nx, ny), ux_ry(nx, ny), uy_lx(nx, ny), uy_rx(nx, ny), uy_ly(nx, ny), uy_ry(nx, ny), ut_lx(nx, ny), ut_rx(nx, ny), ut_ly(nx, ny), ut_ry(nx, ny), eta(nx,ny)            ! Grid spacings in x and y directions

  real(kind=8), intent(out) :: A_lx(nx, ny), A_rx(nx, ny), A_ly(nx, ny), A_ry(nx, ny), Qx_lx(nx, ny), Qx_rx(nx, ny), Qx_ly(nx, ny), Qx_ry(nx, ny), Qy_lx(nx, ny), Qy_rx(nx, ny), Qy_ly(nx, ny), Qy_ry(nx, ny), Qt_lx(nx, ny), Qt_rx(nx, ny), Qt_ly(nx, ny), Qt_ry(nx, ny), sigmaxx_lx(nx, ny), sigmaxx_rx(nx, ny), sigmaxx_ly(nx, ny), sigmaxx_ry(nx, ny), sigmaxy_lx(nx, ny), sigmaxy_rx(nx, ny), sigmaxy_ly(nx, ny), sigmaxy_ry(nx, ny), sigmayy_lx(nx, ny), sigmayy_rx(nx, ny), sigmayy_ly(nx, ny), sigmayy_ry(nx, ny), sigmatx_lx(nx, ny), sigmatx_rx(nx, ny), sigmatx_ly(nx, ny), sigmatx_ry(nx, ny), sigmaty_lx(nx, ny), sigmaty_rx(nx, ny), sigmaty_ly(nx, ny), sigmaty_ry(nx, ny)


  integer :: i, j






  do i = 3, nx-3 
    do j = 3, ny-3
  !Auxiliares para calcular flujos 
  
  ! if (eta(i,j)>0.00001) then


  A_rx(i,j) = 25.0d0/16.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_rx(i,j))*(ux_rx(i,j)*(4.0d0*dtux_rx(i,j)/ut_rx(i,j) + 3.0d0*dxpsi_rx(i,j)) + uy_rx(i,j)*(4.0d0*dtuy_rx(i,j)/ut_rx(i,j) + 3.0d0*dypsi_rx(i,j)) + 4.0d0*dxux_rx(i,j) + 4.0d0*dyuy_rx(i,j) + 3.0d0*dtpsi_rx(i,j)*ut_rx(i,j))
  A_lx(i,j) = 25.0d0/16.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_lx(i,j))*(ux_lx(i,j)*(4.0d0*dtux_lx(i,j)/ut_lx(i,j) + 3.0d0*dxpsi_lx(i,j)) + uy_lx(i,j)*(4.0d0*dtuy_lx(i,j)/ut_lx(i,j) + 3.0d0*dypsi_lx(i,j)) + 4.0d0*dxux_lx(i,j) + 4.0d0*dyuy_lx(i,j) + 3.0d0*dtpsi_lx(i,j)*ut_lx(i,j))
  A_ry(i,j) = 25.0d0/16.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_ry(i,j))*(ux_ry(i,j)*(4.0d0*dtux_ry(i,j)/ut_ry(i,j) + 3.0d0*dxpsi_ry(i,j)) + uy_ry(i,j)*(4.0d0*dtuy_ry(i,j)/ut_ry(i,j) + 3.0d0*dypsi_ry(i,j)) + 4.0d0*dxux_ry(i,j) + 4.0d0*dyuy_ry(i,j) + 3.0d0*dtpsi_ry(i,j)*ut_ry(i,j))
  A_ly(i,j) = 25.0d0/16.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_ly(i,j))*(ux_ly(i,j)*(4.0d0*dtux_ly(i,j)/ut_ly(i,j) + 3.0d0*dxpsi_ly(i,j)) + uy_ly(i,j)*(4.0d0*dtuy_ly(i,j)/ut_ly(i,j) + 3.0d0*dypsi_ly(i,j)) + 4.0d0*dxux_ly(i,j) + 4.0d0*dyuy_ly(i,j) + 3.0d0*dtpsi_ly(i,j)*ut_ly(i,j))

  Qx_rx(i,j) = 25.0d0/28.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_rx(i,j))*(4.0d0*dtux_rx(i,j)*ut_rx(i,j) + ux_rx(i,j)*dtpsi_rx(i,j)*ut_rx(i,j) + (ux_rx(i,j)**2.0d0 + 1.0d0)*dxpsi_rx(i,j) + 4.0d0*ux_rx(i,j)*dxux_rx(i,j) + uy_rx(i,j)*(ux_rx(i,j)*dypsi_rx(i,j) + 4.0d0*dyux_rx(i,j)))
  Qx_lx(i,j) = 25.0d0/28.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_lx(i,j))*(4.0d0*dtux_lx(i,j)*ut_lx(i,j) + ux_lx(i,j)*dtpsi_lx(i,j)*ut_lx(i,j) + (ux_lx(i,j)**2.0d0 + 1.0d0)*dxpsi_lx(i,j) + 4.0d0*ux_lx(i,j)*dxux_lx(i,j) + uy_lx(i,j)*(ux_lx(i,j)*dypsi_lx(i,j) + 4.0d0*dyux_lx(i,j)))
  Qx_ry(i,j) = 25.0d0/28.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_ry(i,j))*(4.0d0*dtux_ry(i,j)*ut_ry(i,j) + ux_ry(i,j)*dtpsi_ry(i,j)*ut_ry(i,j) + (ux_ry(i,j)**2.0d0 + 1.0d0)*dxpsi_ry(i,j) + 4.0d0*ux_ry(i,j)*dxux_ry(i,j) + uy_ry(i,j)*(ux_ry(i,j)*dypsi_ry(i,j) + 4.0d0*dyux_ry(i,j)))
  Qx_ly(i,j) = 25.0d0/28.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_ly(i,j))*(4.0d0*dtux_ly(i,j)*ut_ly(i,j) + ux_ly(i,j)*dtpsi_ly(i,j)*ut_ly(i,j) + (ux_ly(i,j)**2.0d0 + 1.0d0)*dxpsi_ly(i,j) + 4.0d0*ux_ly(i,j)*dxux_ly(i,j) + uy_ly(i,j)*(ux_ly(i,j)*dypsi_ly(i,j) + 4.0d0*dyux_ly(i,j)))

  Qy_rx(i,j) = 25.0d0/28.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_rx(i,j))*(uy_rx(i,j)*(ut_rx(i,j)*dtpsi_rx(i,j) + ux_rx(i,j)*dxpsi_rx(i,j) + uy_rx(i,j)*dypsi_rx(i,j) + 4.0d0*dyuy_rx(i,j)) + 4.0d0*ut_rx(i,j)*dtuy_rx(i,j) + 4.0d0*ux_rx(i,j)*dxuy_rx(i,j) + dypsi_rx(i,j))
  Qy_lx(i,j) = 25.0d0/28.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_lx(i,j))*(uy_lx(i,j)*(ut_lx(i,j)*dtpsi_lx(i,j) + ux_lx(i,j)*dxpsi_lx(i,j) + uy_lx(i,j)*dypsi_lx(i,j) + 4.0d0*dyuy_lx(i,j)) + 4.0d0*ut_lx(i,j)*dtuy_lx(i,j) + 4.0d0*ux_lx(i,j)*dxuy_lx(i,j) + dypsi_lx(i,j))
  Qy_ry(i,j) = 25.0d0/28.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_ry(i,j))*(uy_ry(i,j)*(ut_ry(i,j)*dtpsi_ry(i,j) + ux_ry(i,j)*dxpsi_ry(i,j) + uy_ry(i,j)*dypsi_ry(i,j) + 4.0d0*dyuy_ry(i,j)) + 4.0d0*ut_ry(i,j)*dtuy_ry(i,j) + 4.0d0*ux_ry(i,j)*dxuy_ry(i,j) + dypsi_ry(i,j))
  Qy_ly(i,j) = 25.0d0/28.0d0*eta(i,j)*exp(3.0d0/4.0d0*psi_ly(i,j))*(uy_ly(i,j)*(ut_ly(i,j)*dtpsi_ly(i,j) + ux_ly(i,j)*dxpsi_ly(i,j) + uy_ly(i,j)*dypsi_ly(i,j) + 4.0d0*dyuy_ly(i,j)) + 4.0d0*ut_ly(i,j)*dtuy_ly(i,j) + 4.0d0*ux_ly(i,j)*dxuy_ly(i,j) + dypsi_ly(i,j))

  Qt_rx(i,j) = 1.0d0/ut_rx(i,j)*(ux_rx(i,j)*Qx_rx(i,j) + uy_rx(i,j)*Qy_rx(i,j))
  Qt_lx(i,j) = 1.0d0/ut_lx(i,j)*(ux_lx(i,j)*Qx_lx(i,j) + uy_lx(i,j)*Qy_lx(i,j))
  Qt_ry(i,j) = 1.0d0/ut_ry(i,j)*(ux_ry(i,j)*Qx_ry(i,j) + uy_ry(i,j)*Qy_ry(i,j))
  Qt_ly(i,j) = 1.0d0/ut_ly(i,j)*(ux_ly(i,j)*Qx_ly(i,j) + uy_ly(i,j)*Qy_ly(i,j))

  !esto es tecnicamente -2eta*sigma, pero escribo solo sigma

  sigmaxx_rx(i,j) = 2.0d0/3.0d0/ut_rx(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_rx(i,j))*(ut_rx(i,j)*(-2.0d0*(ux_rx(i,j)**2.0d0 + 1.0d0)*dxux_rx(i,j) + (ux_rx(i,j)**2.0d0)*dyuy_rx(i,j) - 3.0d0*ux_rx(i,j)*dyux_rx(i,j)*uy_rx(i,j) + dyuy_rx(i,j)) - 2.0d0*(ux_rx(i,j)**3)*dtux_rx(i,j) + (ux_rx(i,j)**2.0d0)*uy_rx(i,j)*dtuy_rx(i,j) - ux_rx(i,j)*dtux_rx(i,j)*(3.0d0*(uy_rx(i,j)**2.0d0)+ 2.0d0) + uy_rx(i,j)*dtuy_rx(i,j))
  sigmaxx_lx(i,j) = 2.0d0/3.0d0/ut_lx(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_lx(i,j))*(ut_lx(i,j)*(-2.0d0*(ux_lx(i,j)**2.0d0 + 1.0d0)*dxux_lx(i,j) + (ux_lx(i,j)**2.0d0)*dyuy_lx(i,j) - 3.0d0*ux_lx(i,j)*dyux_lx(i,j)*uy_lx(i,j) + dyuy_lx(i,j)) - 2.0d0*(ux_lx(i,j)**3)*dtux_lx(i,j) + (ux_lx(i,j)**2.0d0)*uy_lx(i,j)*dtuy_lx(i,j) - ux_lx(i,j)*dtux_lx(i,j)*(3.0d0*(uy_lx(i,j)**2.0d0)+ 2.0d0) + uy_lx(i,j)*dtuy_lx(i,j))  
  sigmaxx_ry(i,j) = 2.0d0/3.0d0/ut_ry(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_ry(i,j))*(ut_ry(i,j)*(-2.0d0*(ux_ry(i,j)**2.0d0 + 1.0d0)*dxux_ry(i,j) + (ux_ry(i,j)**2.0d0)*dyuy_ry(i,j) - 3.0d0*ux_ry(i,j)*dyux_ry(i,j)*uy_ry(i,j) + dyuy_ry(i,j)) - 2.0d0*(ux_ry(i,j)**3)*dtux_ry(i,j) + (ux_ry(i,j)**2.0d0)*uy_ry(i,j)*dtuy_ry(i,j) - ux_ry(i,j)*dtux_ry(i,j)*(3.0d0*(uy_ry(i,j)**2.0d0)+ 2.0d0) + uy_ry(i,j)*dtuy_ry(i,j))
  sigmaxx_ly(i,j) = 2.0d0/3.0d0/ut_ly(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_ly(i,j))*(ut_ly(i,j)*(-2.0d0*(ux_ly(i,j)**2.0d0 + 1.0d0)*dxux_ly(i,j) + (ux_ly(i,j)**2.0d0)*dyuy_ly(i,j) - 3.0d0*ux_ly(i,j)*dyux_ly(i,j)*uy_ly(i,j) + dyuy_ly(i,j)) - 2.0d0*(ux_ly(i,j)**3)*dtux_ly(i,j) + (ux_ly(i,j)**2.0d0)*uy_ly(i,j)*dtuy_ly(i,j) - ux_ly(i,j)*dtux_ly(i,j)*(3.0d0*(uy_ly(i,j)**2.0d0)+ 2.0d0) + uy_ly(i,j)*dtuy_ly(i,j))

  sigmaxy_rx(i,j) = -1.0d0/3.0d0/ut_rx(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_rx(i,j))*(ut_rx(i,j)*(3.0d0*(ux_rx(i,j)**2.0d0)*dxuy_rx(i,j) + ux_rx(i,j)*uy_rx(i,j)*(dxux_rx(i,j)+dyuy_rx(i,j)) + 3.0d0*dyux_rx(i,j)*(uy_rx(i,j)**2.0d0 + 1.0d0) + 3.0d0*dxuy_rx(i,j)) + dtux_rx(i,j)*uy_rx(i,j)*(ux_rx(i,j)**2.0d0 + 3.0d0*uy_rx(i,j)**2.0d0 + 3.0d0) + ux_rx(i,j)*dtuy_rx(i,j)*(3.0d0*(ux_rx(i,j)**2.0d0) + uy_rx(i,j)**2.0d0 + 3.0d0))
  sigmaxy_lx(i,j) = -1.0d0/3.0d0/ut_lx(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_lx(i,j))*(ut_lx(i,j)*(3.0d0*(ux_lx(i,j)**2.0d0)*dxuy_lx(i,j) + ux_lx(i,j)*uy_lx(i,j)*(dxux_lx(i,j)+dyuy_lx(i,j)) + 3.0d0*dyux_lx(i,j)*(uy_lx(i,j)**2.0d0 + 1.0d0) + 3.0d0*dxuy_lx(i,j)) + dtux_lx(i,j)*uy_lx(i,j)*(ux_lx(i,j)**2.0d0 + 3.0d0*uy_lx(i,j)**2.0d0 + 3.0d0) + ux_lx(i,j)*dtuy_lx(i,j)*(3.0d0*(ux_lx(i,j)**2.0d0) + uy_lx(i,j)**2.0d0 + 3.0d0))
  sigmaxy_ry(i,j) = -1.0d0/3.0d0/ut_ry(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_ry(i,j))*(ut_ry(i,j)*(3.0d0*(ux_ry(i,j)**2.0d0)*dxuy_ry(i,j) + ux_ry(i,j)*uy_ry(i,j)*(dxux_ry(i,j)+dyuy_ry(i,j)) + 3.0d0*dyux_ry(i,j)*(uy_ry(i,j)**2.0d0 + 1.0d0) + 3.0d0*dxuy_ry(i,j)) + dtux_ry(i,j)*uy_ry(i,j)*(ux_ry(i,j)**2.0d0 + 3.0d0*uy_ry(i,j)**2.0d0 + 3.0d0) + ux_ry(i,j)*dtuy_ry(i,j)*(3.0d0*(ux_ry(i,j)**2.0d0) + uy_ry(i,j)**2.0d0 + 3.0d0))
  sigmaxy_ly(i,j) = -1.0d0/3.0d0/ut_ly(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_ly(i,j))*(ut_ly(i,j)*(3.0d0*(ux_ly(i,j)**2.0d0)*dxuy_ly(i,j) + ux_ly(i,j)*uy_ly(i,j)*(dxux_ly(i,j)+dyuy_ly(i,j)) + 3.0d0*dyux_ly(i,j)*(uy_ly(i,j)**2.0d0 + 1.0d0) + 3.0d0*dxuy_ly(i,j)) + dtux_ly(i,j)*uy_ly(i,j)*(ux_ly(i,j)**2.0d0 + 3.0d0*uy_ly(i,j)**2.0d0 + 3.0d0) + ux_ly(i,j)*dtuy_ly(i,j)*(3.0d0*(ux_ly(i,j)**2.0d0) + uy_ly(i,j)**2.0d0 + 3.0d0))

  sigmayy_rx(i,j) = 2.0d0/3.0d0/ut_rx(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_rx(i,j))*(-3.0d0*ut_rx(i,j)*ux_rx(i,j)*uy_rx(i,j)*dxuy_rx(i,j) + ut_rx(i,j)*dxux_rx(i,j)*(uy_rx(i,j)**2.0d0 + 1.0d0) - 2.0d0*ut_rx(i,j)*(uy_rx(i,j)**2.0d0 + 1.0d0)*dyuy_rx(i,j) - 3.0d0*(ux_rx(i,j)**2.0d0)*uy_rx(i,j)*dtuy_rx(i,j) + ux_rx(i,j)*dtux_rx(i,j)*(uy_rx(i,j)**2.0d0 + 1.0d0) - 2.0d0*(uy_rx(i,j)**3 + uy_rx(i,j))*dtuy_rx(i,j))
  sigmayy_lx(i,j) = 2.0d0/3.0d0/ut_lx(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_lx(i,j))*(-3.0d0*ut_lx(i,j)*ux_lx(i,j)*uy_lx(i,j)*dxuy_lx(i,j) + ut_lx(i,j)*dxux_lx(i,j)*(uy_lx(i,j)**2.0d0 + 1.0d0) - 2.0d0*ut_lx(i,j)*(uy_lx(i,j)**2.0d0 + 1.0d0)*dyuy_lx(i,j) - 3.0d0*(ux_lx(i,j)**2.0d0)*uy_lx(i,j)*dtuy_lx(i,j) + ux_lx(i,j)*dtux_lx(i,j)*(uy_lx(i,j)**2.0d0 + 1.0d0) - 2.0d0*(uy_lx(i,j)**3 + uy_lx(i,j))*dtuy_lx(i,j))
  sigmayy_ry(i,j) = 2.0d0/3.0d0/ut_ry(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_ry(i,j))*(-3.0d0*ut_ry(i,j)*ux_ry(i,j)*uy_ry(i,j)*dxuy_ry(i,j) + ut_ry(i,j)*dxux_ry(i,j)*(uy_ry(i,j)**2.0d0 + 1.0d0) - 2.0d0*ut_ry(i,j)*(uy_ry(i,j)**2.0d0 + 1.0d0)*dyuy_ry(i,j) - 3.0d0*(ux_ry(i,j)**2.0d0)*uy_ry(i,j)*dtuy_ry(i,j) + ux_ry(i,j)*dtux_ry(i,j)*(uy_ry(i,j)**2.0d0 + 1.0d0) - 2.0d0*(uy_ry(i,j)**3 + uy_ry(i,j))*dtuy_ry(i,j))
  sigmayy_ly(i,j) = 2.0d0/3.0d0/ut_ly(i,j)*eta(i,j)*exp(3.0d0/4.0d0*psi_ly(i,j))*(-3.0d0*ut_ly(i,j)*ux_ly(i,j)*uy_ly(i,j)*dxuy_ly(i,j) + ut_ly(i,j)*dxux_ly(i,j)*(uy_ly(i,j)**2.0d0 + 1.0d0) - 2.0d0*ut_ly(i,j)*(uy_ly(i,j)**2.0d0 + 1.0d0)*dyuy_ly(i,j) - 3.0d0*(ux_ly(i,j)**2.0d0)*uy_ly(i,j)*dtuy_ly(i,j) + ux_ly(i,j)*dtux_ly(i,j)*(uy_ly(i,j)**2.0d0 + 1.0d0) - 2.0d0*(uy_ly(i,j)**3 + uy_ly(i,j))*dtuy_ly(i,j))

  sigmatx_rx(i,j) = 1.0d0/ut_rx(i,j)*(ux_rx(i,j)*sigmaxx_rx(i,j) + uy_rx(i,j)*sigmaxy_rx(i,j))
  sigmatx_lx(i,j) = 1.0d0/ut_lx(i,j)*(ux_lx(i,j)*sigmaxx_lx(i,j) + uy_lx(i,j)*sigmaxy_lx(i,j))
  sigmatx_ry(i,j) = 1.0d0/ut_ry(i,j)*(ux_ry(i,j)*sigmaxx_ry(i,j) + uy_ry(i,j)*sigmaxy_ry(i,j))
  sigmatx_ly(i,j) = 1.0d0/ut_ly(i,j)*(ux_ly(i,j)*sigmaxx_ly(i,j) + uy_ly(i,j)*sigmaxy_ly(i,j))

  sigmaty_rx(i,j) = 1.0d0/ut_rx(i,j)*(ux_rx(i,j)*sigmaxy_rx(i,j) + uy_rx(i,j)*sigmayy_rx(i,j))
  sigmaty_lx(i,j) = 1.0d0/ut_lx(i,j)*(ux_lx(i,j)*sigmaxy_lx(i,j) + uy_lx(i,j)*sigmayy_lx(i,j))
  sigmaty_ry(i,j) = 1.0d0/ut_ry(i,j)*(ux_ry(i,j)*sigmaxy_ry(i,j) + uy_ry(i,j)*sigmayy_ry(i,j))
  sigmaty_ly(i,j) = 1.0d0/ut_ly(i,j)*(ux_ly(i,j)*sigmaxy_ly(i,j) + uy_ly(i,j)*sigmayy_ly(i,j))


  ! else

  ! A_lx(i,j) = 0.0d0
  ! A_rx(i,j) = 0.0d0
  ! A_ly(i,j) = 0.0d0
  ! A_ry(i,j) = 0.0d0

  ! Qx_lx(i,j) = 0.0d0
  ! Qx_rx(i,j) = 0.0d0
  ! Qx_ly(i,j) = 0.0d0
  ! Qx_ry(i,j) = 0.0d0

  ! Qy_lx(i,j) = 0.0d0
  ! Qy_rx(i,j) = 0.0d0
  ! Qy_ly(i,j) = 0.0d0
  ! Qy_ry(i,j) = 0.0d0

  ! Qt_lx(i,j) = 0.0d0
  ! Qt_rx(i,j) = 0.0d0
  ! Qt_ly(i,j) = 0.0d0
  ! Qt_ry(i,j) = 0.0d0

  ! sigmaxx_lx(i,j) = 0.0d0
  ! sigmaxx_rx(i,j) = 0.0d0
  ! sigmaxx_ly(i,j) = 0.0d0
  ! sigmaxx_ry(i,j) = 0.0d0


  ! sigmaxy_lx(i,j) = 0.0d0
  ! sigmaxy_rx(i,j) = 0.0d0
  ! sigmaxy_ly(i,j) = 0.0d0
  ! sigmaxy_ry(i,j) = 0.0d0


  ! sigmayy_lx(i,j) = 0.0d0
  ! sigmayy_rx(i,j) = 0.0d0
  ! sigmayy_ly(i,j) = 0.0d0
  ! sigmayy_ry(i,j) = 0.0d0


  ! sigmatx_lx(i,j) = 0.0d0
  ! sigmatx_rx(i,j) = 0.0d0
  ! sigmatx_ly(i,j) = 0.0d0
  ! sigmatx_ry(i,j) = 0.0d0

  ! sigmaty_lx(i,j) = 0.0d0
  ! sigmaty_rx(i,j) = 0.0d0
  ! sigmaty_ly(i,j) = 0.0d0
  ! sigmaty_ry(i,j) = 0.0d0

  ! end if



    enddo
  enddo











end subroutine calculo_funciones_flujos
