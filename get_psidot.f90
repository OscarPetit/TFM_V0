subroutine compute_dtpsi(n, m, psi, ux, uy, dxpsi, dypsi, dxux, dyux, dxuy, dyuy, dtpsi,eta0,q11,q12,q13,deltaeta)
  implicit none
  integer, intent(in) :: n, m
  integer :: i, j
  real(kind=8), intent(inout)  :: psi(n,m), ux(n,m), uy(n,m)
  real(kind=8), intent(in)  :: dxpsi(n,m), dypsi(n,m)
  real(kind=8), intent(in)  :: dxux(n,m), dyux(n,m), dxuy(n,m), dyuy(n,m)
  real(kind=8), intent(inout)  :: eta0(n,m)
  real(kind=8), intent(in) :: q11(n,m), q12(n,m), q13(n,m)
  real(kind=8), intent(out) :: dtpsi(n,m),deltaeta
  real(kind=8) :: XIG, UXG, UYG, ut
  real(kind=8) :: ch, l
  real(kind=8) :: epspf(n,m), psipf(n,m), absvpf(n,m)
  real(kind=8) :: utpf(n,m), uxpf(n,m), uypf(n,m)


  real(kind=8) :: TttPF, TtxPF, TtyPF
  real(kind=8) :: delta_tt, delta_tx, delta_ty
  real(kind=8) :: xiDpf, chi0, lambda0,eta

  chi0 = 25.0d0/4.0d0
  lambda0 = 25.0d0/7.0d0


  ! normalize couplings (chi0, eta, lambda0 are module‐level constants)
  do i = 4, n-3
    do j = 4, m-3
      eta=eta0(i,j)
      XIG = psi(i,j)
      UXG = ux(i,j)
      UYG = uy(i,j)
      if (eta >0.00001d0) then
        ch = chi0 / eta
        l  = lambda0 / eta
      else
      ch = chi0 / (eta + 1.0d0)
      l  = lambda0 / (eta +1.0d0)
      endif

      ! ideal‐fluid moments q01, q02, q03
      ut  = sqrt(1.d0 + UXG**2d0 + UYG**2d0)
      ! q01 = exp(XIG)*(4.d0/3.d0*ut**2 - 1.d0/3.d0)
      ! q02 = exp(XIG)*(4.d0/3.d0*ut*UXG)
      ! q03 = exp(XIG)*(4.d0/3.d0*ut*UYG)

      ! PF “momenta” TttPF, TtxPF, TtyPF (notation mapped: 
      !   uxcx→dxux, uxcy→dyux, uycx→dxuy, uycy→dyuy)
      TttPF = ( &
           9.d0*exp(XIG) &
         + 8.d0*exp(XIG)*UXG**4d0 &
         + 8.d0*exp(XIG)*UYG**4d0 &
         + 6.d0*UYG**2.0d0*(3.d0*exp(XIG) + exp(XIG)**0.75d0*eta*dxux(i,j) - 2.d0*exp(XIG)**0.75d0*eta*dyuy(i,j)) &
         + 3.d0*exp(XIG)**0.75d0*eta*UXG**3.d0*dxpsi(i,j) &
         + 3.d0*exp(XIG)**0.75d0*eta*UXG*UYG*(-6.d0*dyux(i,j) - 6.d0*dxuy(i,j) + UYG*dxpsi(i,j)) &
         + 3.d0*exp(XIG)**0.75d0*eta*UYG**3.d0*dypsi(i,j) &
         + UXG**2.d0*(16.d0*exp(XIG)*UYG**2.d0 &
                  + 6.d0*(3.d0*exp(XIG) - 2.d0*exp(XIG)**0.75d0*eta*dxux(i,j) + exp(XIG)**0.75d0*eta*dyuy(i,j)) &
                  + 3.d0*exp(XIG)**0.75d0*eta*UYG*dypsi(i,j)) &
         ) / (9.d0 + 6.d0*UXG**2 + 6.d0*UYG**2.d0)

      TtxPF = ( &
           32.d0*exp(XIG)*UXG**5 &
         + 12.d0*exp(XIG)**0.75d0*eta*UXG**4*dxpsi(i,j) &
         + 3.d0*exp(XIG)**0.75d0*eta*UYG*(3.d0+2.d0*UYG**2)*(-4.d0*dyux(i,j)-4.d0*dxuy(i,j)+UYG*dxpsi(i,j)) &
         + 6.d0*exp(XIG)**0.75d0*eta*UXG**2*(-8.d0*UYG*(dyux(i,j)+dxuy(i,j)) + 2.d0*dxpsi(i,j) + 3.d0*UYG**2*dxpsi(i,j)) &
         + UXG**3*(80.d0*exp(XIG) - 48.d0*exp(XIG)**0.75d0*eta*dxux(i,j) + 64.d0*exp(XIG)*UYG**2 + 24.d0*exp(XIG)**0.75d0*eta*dyuy(i,j) + 6.d0*exp(XIG)**0.75d0*eta*UYG*dypsi(i,j)) &
         + UXG*(8.d0*(10.d0*exp(XIG) - 3.d0*exp(XIG)**0.75d0*eta*dxux(i,j))*UYG**2 + 32.d0*exp(XIG)*UYG**4 + 24.d0*(2.d0*exp(XIG) - 2.d0*exp(XIG)**0.75d0*eta*dxux(i,j) + exp(XIG)**0.75d0*eta*dyuy(i,j)) + 3.d0*exp(XIG)**0.75d0*eta*UYG*dypsi(i,j) + 6.d0*exp(XIG)**0.75d0*eta*UYG**3.0d0*dypsi(i,j)) &
         ) / (12.d0*sqrt(1.d0+UXG**2+UYG**2)*(3.d0+2.d0*UXG**2+2.d0*UYG**2))

      TtyPF = ( &
           32.d0*exp(XIG)*UYG**5 &
         + UYG**3*(80.d0*exp(XIG) + 64.d0*exp(XIG)*UXG**2 + 24.d0*exp(XIG)**0.75d0*eta*dxux(i,j) - 48.d0*exp(XIG)**0.75d0*eta*dyuy(i,j) + 6.d0*exp(XIG)**0.75d0*eta*UXG*dxpsi(i,j)) &
         + UYG*(32.d0*exp(XIG)*UXG**4 + 8.d0*UXG**2*(10.d0*exp(XIG) - 3.d0*exp(XIG)**0.75d0*eta*dyuy(i,j)) + 24.d0*(2.d0*exp(XIG) + exp(XIG)**0.75d0*eta*dxux(i,j) - 2.d0*exp(XIG)**0.75d0*eta*dyuy(i,j)) + 3.d0*exp(XIG)**0.75d0*eta*UXG*dxpsi(i,j) + 6.d0*exp(XIG)**0.75d0*eta*UXG**3*dxpsi(i,j)) &
         + 12.d0*exp(XIG)**0.75d0*eta*UYG**4*dypsi(i,j) &
         + 3.d0*exp(XIG)**0.75d0*eta*UXG*(3.d0+2.d0*UXG**2)*(-4.d0*dyux(i,j)-4.d0*dxuy(i,j)+UXG*dypsi(i,j)) &
         + 6.d0*exp(XIG)**0.75d0*eta*UYG**2*(-8.d0*UXG*(dyux(i,j)+dxuy(i,j)) + 2.d0*dypsi(i,j) + 3.d0*UXG**2*dypsi(i,j)) &
         ) / (12.d0*sqrt(1.d0+UXG**2+UYG**2)*(3.d0+2.d0*UXG**2+2.d0*UYG**2))

      ! compute deltas
      if (eta >0.00001d0) then
        delta_tt = (TttPF - q11(i,j)) / eta
        delta_tx = (TtxPF - q12(i,j)) / eta
        delta_ty = (TtyPF - q13(i,j)) / eta
      else
        delta_tt = 0.d0
        delta_tx = 0.d0
        delta_ty = 0.d0
      endif

      if(  (eta*abs(delta_tt)<deltaeta) .or. (eta*abs(delta_tx)<deltaeta) .or. (eta*abs(delta_ty)<deltaeta) .or. (eta==0.0d0)) then



        epspf(i,j) = 3.0d0*(-2.0d0/6.0d0*q11(i,j)+sqrt(4.0d0/36.0d0*q11(i,j)**2 + 1.0d0/3.0d0*(q11(i,j)**2.0d0-q12(i,j)**2.0d0-q13(i,j)**2.0d0)))!-q11(i,j) + 3.0d0*sqrt(4.0d0/12.0d0*q11(i,j)**2.0d0 + 1.0d0/3.0d0*(q11(i,j)**2.d0 -q12(i,j)**2.d0-q13(i,j)**2.d0))!!!!!!!!!!!!!!!!!
        psipf(i,j) = log(epspf(i,j))
        absvpf(i,j) = (sqrt(q12(i,j)**2.d0 + q13(i,j)**2.d0))/(q11(i,j) + epspf(i,j)/3.0d0)
        utpf(i,j)=1.0d0/sqrt(1.0d0 - absvpf(i,j)**2.d0)
        uxpf(i,j) = (utpf(i,j)*q12(i,j))/(q11(i,j) + epspf(i,j)/3.0d0)
        uypf(i,j) = (utpf(i,j)*q13(i,j))/(q11(i,j) + epspf(i,j)/3.0d0)!(3.0d0*utpf(i,j)*q13(i,j))/(3.0d0*q11(i,j) + epspf(i,j))
  
        psi(i,j) = psipf(i,j)
        ux(i,j) = uxpf(i,j)
        uy(i,j) = uypf(i,j)

        delta_tt = 0.d0
        delta_tx = 0.d0
        delta_ty = 0.d0

        !eta0(i,j) = 0.d0

        !dtpsi(i,j) = 0.d0
        
      endif


      ! correction term
      xiDpf = -2.d0 * ( &
               2.d0*dyuy(i,j) &
             + UXG**3*dxpsi(i,j) &
             + UXG*(-2.d0*UYG*(dyux(i,j)+dxuy(i,j)) + dxpsi(i,j) + UYG**2*dxpsi(i,j)) &
             + (1.d0+UYG**2.d0)*(2.d0*dxux(i,j) + UYG*dypsi(i,j)) &
             + UXG**2.d0*(2.d0*dyuy(i,j) + UYG*dypsi(i,j)) &
             ) / ( sqrt(1.d0+UXG**2.d0+UYG**2.d0) * (3.d0+2.d0*UXG**2.d0+2.d0*UYG**2.d0) )

      ! final dtpsi
      dtpsi(i,j) = xiDpf &
        - 4.d0*sqrt(1.d0+UXG**2.d0+UYG**2.d0)*(3.d0*l + (-4.d0+4.d0*ch+6.d0*l)*(UXG**2.d0+UYG**2.d0))*delta_tt &
            / ( exp(XIG)**0.75d0*(9.d0*ch*l + 4.d0*(ch*(-3.d0+l)-l)*UXG**4.d0 + 12.d0*ch*(-1.d0+l)*UYG**2.d0 + 4.d0*(ch*(-3.d0+l)-l)*UYG**4.d0 + 4.d0*UXG**2.d0*(3.d0*ch*(-1.d0+l)+2.d0*(ch*(-3.d0+l)-l)*UYG**2.d0)) ) &
        + 4.d0*UXG*(3.d0*(ch+2.d0*l) + (-4.d0+4.d0*ch+6.d0*l)*(UXG**2.d0+UYG**2.d0))*delta_tx &
            / ( exp(XIG)**0.75d0*(9.d0*ch*l + 4.d0*(ch*(-3.d0+l)-l)*UXG**4.d0 + 12.d0*ch*(-1.d0+l)*UYG**2.d0 + 4.d0*(ch*(-3.d0+l)-l)*UYG**4.d0 + 4.d0*UXG**2.d0*(3.d0*ch*(-1.d0+l)+2.d0*(ch*(-3.d0+l)-l)*UYG**2.d0)) ) &
        + 4.d0*UYG*(3.d0*(ch+2.d0*l) + (-4.d0+4.d0*ch+6.d0*l)*(UXG**2.d0+UYG**2.d0))*delta_ty &
            / ( exp(XIG)**0.75d0*(9.d0*ch*l + 4.d0*(ch*(-3.d0+l)-l)*UXG**4.d0 + 12.d0*ch*(-1.d0+l)*UYG**2.d0 + 4.d0*(ch*(-3.d0+l)-l)*UYG**4.d0 + 4.d0*UXG**2.d0*(3.d0*ch*(-1.d0+l)+2.d0*(ch*(-3.d0+l)-l)*UYG**2.d0)) )
    
    

    
          end do
  end do



  !write(*,*) 'dtpsi', maxval(dtpsi)
end subroutine compute_dtpsi
