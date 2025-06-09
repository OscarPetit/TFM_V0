subroutine compute_dtux(n, m, psi, ux, uy, dxpsi, dypsi, dxux, dyux, dxuy, dyuy, dtux, eta, q11, q12, q13, deltaeta)
  implicit none
  integer, intent(in) :: n, m
  integer :: i, j

  real(kind=8), intent(in)  :: psi(n,m), ux(n,m), uy(n,m)
  real(kind=8), intent(in)  :: dxpsi(n,m), dypsi(n,m)
  real(kind=8), intent(in)  :: dxux(n,m), dyux(n,m), dxuy(n,m), dyuy(n,m)
  real(kind=8), intent(inout)  :: eta(n,m)
  real(kind=8), intent(in) :: q11(n,m), q12(n,m), q13(n,m)
  real(kind=8), intent(out) :: dtux(n,m),deltaeta

  real(kind=8) :: chi0, lambda0
  real(kind=8) :: XIG, UXG, UYG, ch, l

  real(kind=8) :: TttPF, TtxPF, TtyPF
  real(kind=8) :: delta_tt, delta_tx, delta_ty
  real(kind=8) :: uxDpf, common, denomPF,eta0

  !-- constants
  chi0    = 25.0d0/4.0d0
  lambda0 = 25.0d0/7.0d0

  do i = 4, n-3
    do j = 4, m-3
      eta0=eta(i,j)
      ! map local notation with "G" suffix
      XIG = psi(i,j)
      UXG = ux(i,j)
      UYG = uy(i,j)
      if (eta0 >0.00001) then
        ch = chi0 / eta0
        l  = lambda0 / eta0
      else
        ch = chi0 / (eta0 + 1.0)
        l  = lambda0 / (eta0 +1.0)
      endif


      ! compute ideal‐fluid moments q01,q02,q03
      common = sqrt(1.0d0 + UXG**2 + UYG**2)
      ! q01    = exp(XIG)*(4.0d0/3.0d0*common**2 - 1.0d0/3.0d0)
      ! q02    = exp(XIG)*(4.0d0/3.0d0*common*UXG)
      ! q03    = exp(XIG)*(4.0d0/3.0d0*common*UYG)

      ! compute PF momenta TttPF, TtxPF, TtyPF
      denomPF = 9.0d0 + 6.0d0*UXG**2 + 6.0d0*UYG**2
      TttPF = ( &
           9.0d0*exp(XIG) &
         + 8.0d0*exp(XIG)*UXG**4 &
         + 8.0d0*exp(XIG)*UYG**4 &
         + 6.0d0*UYG**2*(3.0d0*exp(XIG) + exp(XIG)**0.75d0*eta0*dxux(i,j) - 2.0d0*exp(XIG)**0.75d0*eta0*dyuy(i,j)) &
         + 3.0d0*exp(XIG)**0.75d0*eta0*UXG**3*dxpsi(i,j) &
         + 3.0d0*exp(XIG)**0.75d0*eta0*UXG*UYG*(-6.0d0*dyux(i,j) - 6.0d0*dxuy(i,j) + UYG*dxpsi(i,j)) &
         + 3.0d0*exp(XIG)**0.75d0*eta0*UYG**3*dypsi(i,j) &
         + UXG**2*(16.0d0*exp(XIG)*UYG**2 &
                  + 6.0d0*(3.0d0*exp(XIG) - 2.0d0*exp(XIG)**0.75d0*eta0*dxux(i,j) + exp(XIG)**0.75d0*eta0*dyuy(i,j)) &
                  + 3.0d0*exp(XIG)**0.75d0*eta0*UYG*dypsi(i,j)) &
        ) / denomPF

      TtxPF = ( &
           32.0d0*exp(XIG)*UXG**5 &
         + 12.0d0*exp(XIG)**0.75d0*eta0*UXG**4*dxpsi(i,j) &
         + 3.0d0*exp(XIG)**0.75d0*eta0*UYG*(3.0d0+2.0d0*UYG**2)*(-4.0d0*dyux(i,j)-4.0d0*dxuy(i,j)+UYG*dxpsi(i,j)) &
         + 6.0d0*exp(XIG)**0.75d0*eta0*UXG**2*(-8.0d0*UYG*(dyux(i,j)+dxuy(i,j)) + 2.0d0*dxpsi(i,j) + 3.0d0*UYG**2*dxpsi(i,j)) &
         + UXG**3*(80.0d0*exp(XIG) - 48.0d0*exp(XIG)**0.75d0*eta0*dxux(i,j) + 64.0d0*exp(XIG)*UYG**2 + 24.0d0*exp(XIG)**0.75d0*eta0*dyuy(i,j) + 6.0d0*exp(XIG)**0.75d0*eta0*UYG*dypsi(i,j)) &
         + UXG*(8.0d0*(10.0d0*exp(XIG) - 3.0d0*exp(XIG)**0.75d0*eta0*dxux(i,j))*UYG**2 + 32.0d0*exp(XIG)*UYG**4 + 24.0d0*(2.0d0*exp(XIG) - 2.0d0*exp(XIG)**0.75d0*eta0*dxux(i,j) + exp(XIG)**0.75d0*eta0*dyuy(i,j)) + 3.0d0*exp(XIG)**0.75d0*eta0*UYG*dypsi(i,j) + 6.0d0*exp(XIG)**0.75d0*eta0*UYG**3*dypsi(i,j)) &
        ) / (12.0d0*sqrt(1.0d0+UXG**2+UYG**2)*(3.0d0+2.0d0*UXG**2+2.0d0*UYG**2))

      TtyPF = ( &
           32.0d0*exp(XIG)*UYG**5 &
         + UYG**3*(80.0d0*exp(XIG) + 64.0d0*exp(XIG)*UXG**2 + 24.0d0*exp(XIG)**0.75d0*eta0*dxux(i,j) - 48.0d0*exp(XIG)**0.75d0*eta0*dyuy(i,j) + 6.0d0*exp(XIG)**0.75d0*eta0*UXG*dxpsi(i,j)) &
         + UYG*(32.0d0*exp(XIG)*UXG**4 + 8.0d0*UXG**2*(10.0d0*exp(XIG) - 3.0d0*exp(XIG)**0.75d0*eta0*dyuy(i,j)) + 24.0d0*(2.0d0*exp(XIG) + exp(XIG)**0.75d0*eta0*dxux(i,j) - 2.0d0*exp(XIG)**0.75d0*eta0*dyuy(i,j)) + 3.0d0*exp(XIG)**0.75d0*eta0*UXG*dxpsi(i,j) + 6.0d0*exp(XIG)**0.75d0*eta0*UXG**3*dxpsi(i,j)) &
         + 12.0d0*exp(XIG)**0.75d0*eta0*UYG**4*dypsi(i,j) &
         + 3.0d0*exp(XIG)**0.75d0*eta0*UXG*(3.0d0+2.0d0*UXG**2)*(-4.0d0*dyux(i,j)-4.0d0*dxuy(i,j)+UXG*dypsi(i,j)) &
         + 6.0d0*exp(XIG)**0.75d0*eta0*UYG**2*(-8.0d0*UXG*(dyux(i,j)+dxuy(i,j)) + 2.0d0*dypsi(i,j) + 3.0d0*UXG**2*dypsi(i,j)) &
        ) / (12.0d0*sqrt(1.0d0+UXG**2+UYG**2)*(3.0d0+2.0d0*UXG**2+2.0d0*UYG**2))

      ! compute deltas (floor applied earlier in compute_dtpsi)
        if (eta0 >0.00001) then
          delta_tt = (TttPF - q11(i,j)) / eta0
          delta_tx = (TtxPF - q12(i,j)) / eta0
          delta_ty = (TtyPF - q13(i,j)) / eta0
        else
          delta_tt = 0.d0
          delta_tx = 0.d0
          delta_ty = 0.d0
        endif
  
        if(  (eta0*abs(delta_tt)<deltaeta) .or. (eta0*abs(delta_tx)<deltaeta) .or. (eta0*abs(delta_ty)<deltaeta) .or. (eta0==0.0)) then
  
  
          delta_tt = 0.d0
          delta_tx = 0.d0
          delta_ty = 0.d0
  
          !eta(i,j) = 0.d0
  
          !dtux(i,j) = 0.d0
        endif


      ! PF correction for ux
      uxDpf = ( &
         -8.0d0*dyux(i,j)*UYG**3 &
         + (1.0d0 + UXG**2)*(4.0d0*UXG*(-2.0d0*dxux(i,j) + dyuy(i,j)) - 3.0d0*dxpsi(i,j)) &
         - 2.0d0*UYG**2*(2.0d0*UXG*dxux(i,j) + dxpsi(i,j)) &
         - UYG*(12.0d0*(1.0d0 + UXG**2)*dyux(i,j) + UXG*(4.0d0*UXG*dxuy(i,j) + dypsi(i,j))) &
        ) / (4.0d0*sqrt(1.0d0+UXG**2+UYG**2)*(3.0d0+2.0d0*UXG**2+2.0d0*UYG**2))

      ! assemble final dtux
      dtux(i,j) = uxDpf &
        + (3.0d0*UXG*sqrt(1.0d0+UXG**2+UYG**2)*(4.0d0*ch + l + 2.0d0*(2.0d0*ch + l)*(UXG**2+UYG**2))*delta_tt) &
            / (exp(XIG)**0.75d0*(9.0d0*ch*l + 4.0d0*(ch*(-3.0d0+l)-l)*UXG**4 + 12.0d0*ch*(-1.0d0+l)*UYG**2 + 4.0d0*(ch*(-3.0d0+l)-l)*UYG**4 + 4.0d0*UXG**2*(3.0d0*ch*(-1.0d0+l)+2.0d0*(ch*(-3.0d0+l)-l)*UYG**2))) &
        + ((-9.0d0*ch*l - 6.0d0*(-1.0d0+l)*(2.0d0*ch + l)*UXG**6 - 12.0d0*ch*(-1.0d0+l)*UYG**2 + 4.0d0*(-(ch*(-3.0d0+l)) + l)*UYG**4 - 3.0d0*UXG**4*(2.0d0*l*(-1.0d0+2.0d0*l) + ch*(-7.0d0+11.0d0*l) + 4.0d0*(-1.0d0+l)*(2.0d0*ch + l)*UYG**2) + UXG**2*(ch*(9.0d0-30.0d0*l) - 6.0d0*l**2 + (ch*(33.0d0-37.0d0*l) + 2.0d0*(5.0d0-6.0d0*l)*l)*UYG**2 - 6.0d0*(-1.0d0+l)*(2.0d0*ch + l)*UYG**4)) * delta_tx) &
            / (exp(XIG)**0.75d0*(l + (-1.0d0+l)*(UXG**2+UYG**2))*(9.0d0*ch*l + 4.0d0*(ch*(-3.0d0+l)-l)*UXG**4 + 12.0d0*ch*(-1.0d0+l)*UYG**2 + 4.0d0*(ch*(-3.0d0+l)-l)*UYG**4 + 4.0d0*UXG**2*(3.0d0*ch*(-1.0d0+l)+2.0d0*(ch*(-3.0d0+l)-l)*UYG**2))) &
        - (UXG*UYG*(3.0d0*(ch + 6.0d0*ch*l + 2.0d0*l**2) + 6.0d0*(-1.0d0+l)*(2.0d0*ch + l)*UXG**4 + (2.0d0*l*(-1.0d0+6.0d0*l) + ch*(-9.0d0+29.0d0*l))*UYG**2 + 6.0d0*(-1.0d0+l)*(2.0d0*ch + l)*UYG**4 + UXG**2*(2.0d0*l*(-1.0d0+6.0d0*l) + ch*(-9.0d0+29.0d0*l) + 12.0d0*(-1.0d0+l)*(2.0d0*ch + l)*UYG**2)) * delta_ty) &
            / (exp(XIG)**0.75d0*(l + (-1.0d0+l)*(UXG**2+UYG**2))*(9.0d0*ch*l + 4.0d0*(ch*(-3.0d0+l)-l)*UXG**4 + 12.0d0*ch*(-1.0d0+l)*UYG**2 + 4.0d0*(ch*(-3.0d0+l)-l)*UYG**4 + 4.0d0*UXG**2*(3.0d0*ch*(-1.0d0+l)+2.0d0*(ch*(-3.0d0+l)-l)*UYG**2)))

    
          end do
  end do
!write(*,*) "dtux", maxval(dtux)
  ! write(*,*) "dtux", dtux(4:10,4:10)
  ! write(*,*) "dtux", dtux(4:10,4:10) - dtux(4:10,4:10)
  ! write(*,*) "dtux", dtux(4:10,4:10) - dtux(4:10,4:10)

end subroutine compute_dtux
