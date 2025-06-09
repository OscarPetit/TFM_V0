subroutine weno_2d_xL(Nx, Ny, p, pLx)
        implicit none
        integer, intent(in) :: Nx, Ny
        real(kind=8), intent(in) :: p(1:Nx,1:Ny)
        real(kind=8), intent(out):: pLx(1:Nx,1:Ny)
        integer :: i,j
        real(kind=8), parameter :: epsW=1.0d0
        
        real(kind=8), parameter :: d2=3.d0/10.d0, d1=3.d0/5.d0, d0=1.d0/10.d0
        real(kind=8) :: u0,u1,u2, beta0,beta1,beta2
        real(kind=8) :: a0,a1,a2,sumA, w0,w1,w2
      
        do i=3, Nx-3
          do j=3, Ny-3
            ! 3rd-order left-biased
            u0 = (2.d0*p(i-2,j)-7.d0*p(i-1,j)+11.d0*p(i,j))/6.d0
            u1 = (-1.d0*p(i-1,j)+5.d0*p(i,j)+2.d0*p(i+1,j))/6.d0
            u2 = (2.d0*p(i,j)+5.d0*p(i+1,j)-1.d0*p(i+2,j))/6.d0
            ! smoothness
            beta0=(13.d0/12.d0)*(p(i-2,j)-2.d0*p(i-1,j)+p(i,j))**2 + (1.d0/4.d0)*(p(i-2,j)-4.d0*p(i-1,j)+3.d0*p(i,j))**2
            beta1=(13.d0/12.d0)*(p(i-1,j)-2.d0*p(i,j)+p(i+1,j))**2 + (1.d0/4.d0)*(p(i+1,j)-p(i-1,j))**2
            beta2=(13.d0/12.d0)*(p(i,j)-2.d0*p(i+1,j)+p(i+2,j))**2 + (1.d0/4.d0)*(3.d0*p(i,j)-4.d0*p(i+1,j)+p(i+2,j))**2
            a0=d0/((epsW+beta0)**2); a1=d1/((epsW+beta1)**2); a2=d2/((epsW+beta2)**2)
            sumA=a0+a1+a2; w0=a0/sumA; w1=a1/sumA; w2=a2/sumA
            pLx(i,j)=w0*u0 + w1*u1 + w2*u2
          end do
        end do
      end subroutine weno_2d_xL