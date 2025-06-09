subroutine flujos_yR(eta, nx, ny, psi_ly, psi_ry, ux_ly, ux_ry, uy_ly, uy_ry, ut_ly, ut_ry, A_ly, A_ry, Qx_ly, Qx_ry, Qy_ly, Qy_ry, Qt_ly, Qt_ry, sigmaty_ly, sigmaty_ry, sigmaxy_ly, sigmaxy_ry, sigmayy_ly, sigmayy_ry, dens_ly, dens_ry, fyt_lyR, fyt_ryR, fyx_lyR, fyx_ryR, fyy_lyR, fyy_ryR, fdensy_lyR, fdensy_ryR)

  implicit none
  
  !----------------------------------------------------------------------
  !  Input parameters
  !----------------------------------------------------------------------
  integer, intent(in) :: nx, ny
  
  ! Psi fields
  real(kind=8), intent(in) :: psi_ly(nx, ny), psi_ry(nx, ny), eta(nx,ny)

  
  ! Velocities for "current" cell
  real(kind=8), intent(in) :: ux_ly(nx, ny), uy_ly(nx, ny), ut_ly(nx, ny)
  real(kind=8), intent(in) :: ux_ry(nx, ny), uy_ry(nx, ny), ut_ry(nx, ny)
  
  
  
  ! Amplitude fields
  real(kind=8), intent(in) :: A_ly(nx, ny), A_ry(nx, ny)
  
  ! Q fields
  real(kind=8), intent(in) :: Qx_ly(nx, ny), Qx_ry(nx, ny)
  real(kind=8), intent(in) :: Qy_ly(nx, ny), Qy_ry(nx, ny)
  
  ! Time-direction Q (Qt) fields
  real(kind=8), intent(in) :: Qt_ly(nx, ny), Qt_ry(nx, ny)
  
  ! Stress tensors
  real(kind=8), intent(in) :: sigmaty_ly(nx, ny), sigmaty_ry(nx, ny)
  real(kind=8), intent(in) :: sigmaxy_ly(nx, ny), sigmaxy_ry(nx, ny)
  real(kind=8), intent(in) :: sigmayy_ly(nx, ny), sigmayy_ry(nx, ny)
  
  ! Density fields
  real(kind=8), intent(in) :: dens_ly(nx, ny), dens_ry(nx, ny)
  
  !----------------------------------------------------------------------
  !  Output flux arrays
  !----------------------------------------------------------------------
  real(kind=8), intent(out) :: fyt_lyR(nx, ny), fyt_ryR(nx, ny)
  real(kind=8), intent(out) :: fyx_lyR(nx, ny), fyx_ryR(nx, ny)
  real(kind=8), intent(out) :: fyy_lyR(nx, ny), fyy_ryR(nx, ny)
  real(kind=8), intent(out) :: fdensy_lyR(nx, ny), fdensy_ryR(nx, ny)
  
  !----------------------------------------------------------------------
  !  Local variables
  !----------------------------------------------------------------------
  integer :: i, j
  
  
  do i = 3, nx-3 
    do j = 3, ny-3

      !if(eta(i,j) > 0.0001d0) then
  
        fyt_ryR(i,j) = exp(psi_ry(i,j))*(4.0d0/3.0d0*uy_ry(i,j)*ut_ry(i,j)) + A_ry(i,j)*(4.0d0/3.0d0*uy_ry(i,j)*ut_ry(i,j)) + ut_ry(i,j)*Qy_ry(i,j) + uy_ry(i,j)*Qt_ry(i,j) + sigmaty_ry(i,j)
        fyt_lyR(i,j) = exp(psi_ly(i,j))*(4.0d0/3.0d0*uy_ly(i,j)*ut_ly(i,j)) + A_ly(i,j)*(4.0d0/3.0d0*uy_ly(i,j)*ut_ly(i,j)) + ut_ly(i,j)*Qy_ly(i,j) + uy_ly(i,j)*Qt_ly(i,j) + sigmaty_ly(i,j)
    
        fyx_ryR(i,j) = (exp(psi_ry(i,j)) + A_ry(i,j))*(4.0d0/3.0d0*uy_ry(i,j)*ux_ry(i,j)) + uy_ry(i,j)*Qx_ry(i,j) + ux_ry(i,j)*Qy_ry(i,j) + sigmaxy_ry(i,j)
        fyx_lyR(i,j) = (exp(psi_ly(i,j)) + A_ly(i,j))*(4.0d0/3.0d0*uy_ly(i,j)*ux_ly(i,j)) + uy_ly(i,j)*Qx_ly(i,j) + ux_ly(i,j)*Qy_ly(i,j) + sigmaxy_ly(i,j)
    
        fyy_ryR(i,j) = (exp(psi_ry(i,j)) + A_ry(i,j))*(4.0d0/3.0d0*uy_ry(i,j)**2.0d0 + 1.0d0/3.0d0) + 2.0d0*uy_ry(i,j)*Qy_ry(i,j) + sigmayy_ry(i,j)
        fyy_lyR(i,j) = (exp(psi_ly(i,j)) + A_ly(i,j))*(4.0d0/3.0d0*uy_ly(i,j)**2.0d0 + 1.0d0/3.0d0) + 2.0d0*uy_ly(i,j)*Qy_ly(i,j) + sigmayy_ly(i,j)
    
        fdensy_ryR(i,j) = dens_ry(i,j)*uy_ry(i,j)
        fdensy_lyR(i,j) = dens_ly(i,j)*uy_ly(i,j)
    
  
      ! else
      !   fyt_ryR(i,j) = exp(psi_ry(i,j))*(4.0d0/3.0d0*uy_ry(i,j)*ut_ry(i,j))
      !   fyt_lyR(i,j) = exp(psi_ly(i,j))*(4.0d0/3.0d0*uy_ly(i,j)*ut_ly(i,j))

      !   fyx_ryR(i,j) = exp(psi_ry(i,j))*(4.0d0/3.0d0*uy_ry(i,j)*ux_ry(i,j))
      !   fyx_lyR(i,j) = exp(psi_ly(i,j))*(4.0d0/3.0d0*uy_ly(i,j)*ux_ly(i,j))

      !   fyy_ryR(i,j) = exp(psi_ry(i,j))*(4.0d0/3.0d0*uy_ry(i,j)*uy_ry(i,j) + 1.0/3.0)
      !   fyy_lyR(i,j) = exp(psi_ly(i,j))*(4.0d0/3.0d0*uy_ly(i,j)*uy_ly(i,j) + 1.0/3.0)

      !   fdensy_ryR(i,j) = dens_ry(i,j)*uy_ry(i,j)
      !   fdensy_lyR(i,j) = dens_ly(i,j)*uy_ly(i,j)

      ! endif
      enddo
    enddo
  
  
  
  
  
  
  
  
  
  end subroutine flujos_yR
  