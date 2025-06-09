subroutine weno_2d_yR(Nx, Ny, p, pRy)
        implicit none
        integer, intent(in) :: Nx, Ny
        real(kind=8), intent(in) :: p(1:Nx,1:Ny)
        real(kind=8), intent(out):: pRy(1:Nx,1:Ny)
        integer :: i,j

        real(kind=8), parameter :: epsW=1.0d0
        
        real(kind=8), parameter :: d0=1.d0/10.d0, d1=3.d0/5.d0, d2=3.d0/10.d0
        real(kind=8) :: v0,v1,v2, beta0,beta1,beta2
        real(kind=8) :: a0,a1,a2,sumA, w0,w1,w2
      
        do i=3, Nx-3
          do j=3, Ny-3
            ! 3rd-order right-biased in y
            v0 = (2.d0*p(i,j+3)-7.d0*p(i,j+2)+11.d0*p(i,j+1))/6.d0
            v1 = (-1.d0*p(i,j+2)+5.d0*p(i,j+1)+2.d0*p(i,j))/6.d0
            v2 = (2.d0*p(i,j+1)+5.d0*p(i,j)-1.d0*p(i,j-1))/6.d0
            ! smoothness
            beta0=(13.d0/12.d0)*(p(i,j+1)-2.d0*p(i,j+2)+p(i,j+3))**2 + (1.d0/4.d0)*(3.d0*p(i,j+1)-4.d0*p(i,j+2)+p(i,j+3))**2
            beta1=(13.d0/12.d0)*(p(i,j)-2.d0*p(i,j+1)+p(i,j+2))**2 + (1.d0/4.d0)*(p(i,j+2)-p(i,j))**2
            beta2=(13.d0/12.d0)*(p(i,j-1)-2.d0*p(i,j)+p(i,j+1))**2 + (1.d0/4.d0)*(p(i,j-1)-4.d0*p(i,j)+3.d0*p(i,j+1))**2
            a0=d0/((epsW+beta0)**2); a1=d1/((epsW+beta1)**2); a2=d2/((epsW+beta2)**2)
            sumA=a0+a1+a2; w0=a0/sumA; w1=a1/sumA; w2=a2/sumA
            pRy(i,j)=w0*v0 + w1*v1 + w2*v2
          end do
        end do
      end subroutine weno_2d_yR



! subroutine weno5_y(nx, ny, p, pL, pR)
  
!   implicit none

!   integer, intent(in)    :: nx, ny
!   real(kind=8), intent(in)  :: p   (1:nx,1:ny)
!   real(kind=8), intent(out) :: pL  (1:nx,1:ny)
!   real(kind=8), intent(out) :: pR  (1:nx,1:ny)

!   integer :: i, j
!   real(kind=8), parameter :: eps = 1.0d-3
!   real(kind=8), parameter :: d0 = 0.1d0, d1 = 0.6d0, d2 = 0.3d0
!   real(kind=8) :: h0, h1, h2
!   real(kind=8) :: beta0, beta1, beta2
!   real(kind=8) :: alpha0, alpha1, alpha2, wsum
!   real(kind=8) :: w0, w1, w2
!   real(kind=8) :: wm0, wm1, wm2

!   do i = 1, nx
!     do j = 3, ny-3
!       ! Left state p^- at (i, j+1/2) using WENO-M
!       h0 = (  2.0d0*p(i,j-2) -  7.0d0*p(i,j-1) + 11.0d0*p(i,j  ) ) / 6.0d0
!       h1 = ( -1.0d0*p(i,j-1) +  5.0d0*p(i,j  ) +  2.0d0*p(i,j+1) ) / 6.0d0
!       h2 = (  2.0d0*p(i,j  ) +  5.0d0*p(i,j+1) -  1.0d0*p(i,j+2) ) / 6.0d0

!       beta0 = (13.0d0/12.0d0)*(p(i,j-2)-2.0d0*p(i,j-1)+p(i,j  ))**2 + &
!               0.25d0*(p(i,j-2)-4.0d0*p(i,j-1)+3.0d0*p(i,j  ))**2
!       beta1 = (13.0d0/12.0d0)*(p(i,j-1)-2.0d0*p(i,j  )+p(i,j+1))**2 + &
!               0.25d0*(p(i,j-1)           -   p(i,j+1))**2
!       beta2 = (13.0d0/12.0d0)*(p(i,j  )-2.0d0*p(i,j+1)+p(i,j+2))**2 + &
!               0.25d0*(3.0d0*p(i,j  )-4.0d0*p(i,j+1)+p(i,j+2))**2

!       alpha0 = d2 / ( (eps + beta0)**2 )
!       alpha1 = d1 / ( (eps + beta1)**2 )
!       alpha2 = d0 / ( (eps + beta2)**2 )
!       wsum   = alpha0 + alpha1 + alpha2
!       w0     = alpha0 / wsum
!       w1     = alpha1 / wsum
!       w2     = alpha2 / wsum

!       wm0 = ( w0 * (d0 + d0**2 - 3.0d0*d0*w0 + w0**2) ) / ( d0**2 + w0*(1.0d0 - 2.0d0*d0) )
!       wm1 = ( w1 * (d1 + d1**2 - 3.0d0*d1*w1 + w1**2) ) / ( d1**2 + w1*(1.0d0 - 2.0d0*d1) )
!       wm2 = ( w2 * (d2 + d2**2 - 3.0d0*d2*w2 + w2**2) ) / ( d2**2 + w2*(1.0d0 - 2.0d0*d2) )
!       wsum = wm0 + wm1 + wm2
!       w0   = wm0 / wsum
!       w1   = wm1 / wsum
!       w2   = wm2 / wsum

!       pL(i,j) = w0*h0 + w1*h1 + w2*h2

!       ! Right state p^+ at (i, j+1/2)
!       h0 = (  2.0d0*p(i,j+3) -  7.0d0*p(i,j+2) + 11.0d0*p(i,j+1) ) / 6.0d0
!       h1 = ( -1.0d0*p(i,j+2) +  5.0d0*p(i,j+1) +  2.0d0*p(i,j  ) ) / 6.0d0
!       h2 = (  2.0d0*p(i,j+1) +  5.0d0*p(i,j  ) -  1.0d0*p(i,j-1) ) / 6.0d0

!       beta0 = (13.0d0/12.0d0)*(p(i,j+3)-2.0d0*p(i,j+2)+p(i,j+1))**2 + &
!               0.25d0*(p(i,j+3)-4.0d0*p(i,j+2)+3.0d0*p(i,j+1))**2
!       beta1 = (13.0d0/12.0d0)*(p(i,j+2)-2.0d0*p(i,j+1)+p(i,j  ))**2 + &
!               0.25d0*(p(i,j+2)           -   p(i,j  ))**2
!       beta2 = (13.0d0/12.0d0)*(p(i,j+1)-2.0d0*p(i,j  )+p(i,j-1))**2 + &
!               0.25d0*(3.0d0*p(i,j+1)-4.0d0*p(i,j  )+p(i,j-1))**2

!       alpha0 = d0 / ( (eps + beta0)**2 )
!       alpha1 = d1 / ( (eps + beta1)**2 )
!       alpha2 = d2 / ( (eps + beta2)**2 )
!       wsum   = alpha0 + alpha1 + alpha2
!       w0     = alpha0 / wsum
!       w1     = alpha1 / wsum
!       w2     = alpha2 / wsum

!       wm0 = ( w0 * (d0 + d0**2 - 3.0d0*d0*w0 + w0**2) ) / ( d0**2 + w0*(1.0d0 - 2.0d0*d0) )
!       wm1 = ( w1 * (d1 + d1**2 - 3.0d0*d1*w1 + w1**2) ) / ( d1**2 + w1*(1.0d0 - 2.0d0*d1) )
!       wm2 = ( w2 * (d2 + d2**2 - 3.0d0*d2*w2 + w2**2) ) / ( d2**2 + w2*(1.0d0 - 2.0d0*d2) )
!       wsum = wm0 + wm1 + wm2
!       w0   = wm0 / wsum
!       w1   = wm1 / wsum
!       w2   = wm2 / wsum

!       pR(i,j) = w0*h0 + w1*h1 + w2*h2
!     end do
!   end do
  
!   ! implicit none

!   ! integer, intent(in)    :: nx, ny
!   ! real(kind=8), intent(in)  :: p   (1:nx,1:ny)
!   ! real(kind=8), intent(out) :: pL  (1:nx,1:ny)
!   ! real(kind=8), intent(out) :: pR  (1:nx,1:ny)

!   ! integer :: i, j
!   ! real(kind=8), parameter :: eps = 1.0d-10
!   ! real(kind=8), parameter :: d0 = 0.1d0, d1 = 0.6d0, d2 = 0.3d0

!   ! real(kind=8) :: h0, h1, h2
!   ! real(kind=8) :: beta0, beta1, beta2, tau5
!   ! real(kind=8) :: alpha0, alpha1, alpha2, wsum
!   ! real(kind=8) :: w0, w1, w2

!   ! do i = 1, nx
!   !   do j = 3, ny-3

!   !     !------------------------------------------------------
!   !     ! 1) LEFT STATE p^- at y-interface (i, j+1/2) from {j-2..j+2}
!   !     !------------------------------------------------------
!   !     h0 = (  2.0d0*p(i,j-2) -  7.0d0*p(i,j-1) + 11.0d0*p(i,j  ) ) / 6.0d0
!   !     h1 = ( -1.0d0*p(i,j-1) +  5.0d0*p(i,j  ) +  2.0d0*p(i,j+1) ) / 6.0d0
!   !     h2 = (  2.0d0*p(i,j  ) +  5.0d0*p(i,j+1) -  1.0d0*p(i,j+2) ) / 6.0d0

!   !     beta0 = (13.0d0/12.0d0)*(p(i,j-2)-2.0d0*p(i,j-1)+p(i,j  ))**2 &
!   !           +  0.25d0*(p(i,j-2)-4.0d0*p(i,j-1)+3.0d0*p(i,j  ))**2
!   !     beta1 = (13.0d0/12.0d0)*(p(i,j-1)-2.0d0*p(i,j  )+p(i,j+1))**2 &
!   !           +  0.25d0*(p(i,j-1)           -   p(i,j+1)        )**2
!   !     beta2 = (13.0d0/12.0d0)*(p(i,j  )-2.0d0*p(i,j+1)+p(i,j+2))**2 &
!   !           +  0.25d0*(3.0d0*p(i,j  )-4.0d0*p(i,j+1)+p(i,j+2))**2

!   !     tau5 = abs(beta0 - beta2)

!   !     alpha0 = d0 * (1.0d0 + tau5 / (beta0 + eps))
!   !     alpha1 = d1 * (1.0d0 + tau5 / (beta1 + eps))
!   !     alpha2 = d2 * (1.0d0 + tau5 / (beta2 + eps))
!   !     wsum   = alpha0 + alpha1 + alpha2

!   !     w0 = alpha0 / wsum
!   !     w1 = alpha1 / wsum
!   !     w2 = alpha2 / wsum

!   !     pL(i,j) = w0*h0 + w1*h1 + w2*h2


!   !     !------------------------------------------------------
!   !     ! 2) RIGHT STATE p^+ at y-interface (i, j+1/2) from {j-1..j+3}
!   !     !------------------------------------------------------
!   !     h0 = (  2.0d0*p(i,j+3) -  7.0d0*p(i,j+2) + 11.0d0*p(i,j+1) ) / 6.0d0
!   !     h1 = ( -1.0d0*p(i,j+2) +  5.0d0*p(i,j+1) +  2.0d0*p(i,j  ) ) / 6.0d0
!   !     h2 = (  2.0d0*p(i,j+1) +  5.0d0*p(i,j  ) -  1.0d0*p(i,j-1) ) / 6.0d0

!   !     beta0 = (13.0d0/12.0d0)*(p(i,j+3)-2.0d0*p(i,j+2)+p(i,j+1))**2 &
!   !           +  0.25d0*(p(i,j+3)-4.0d0*p(i,j+2)+3.0d0*p(i,j+1))**2
!   !     beta1 = (13.0d0/12.0d0)*(p(i,j+2)-2.0d0*p(i,j+1)+p(i,j  ))**2 &
!   !           +  0.25d0*(p(i,j+2)           -   p(i,j  )        )**2
!   !     beta2 = (13.0d0/12.0d0)*(p(i,j+1)-2.0d0*p(i,j  )+p(i,j-1))**2 &
!   !           +  0.25d0*(3.0d0*p(i,j+1)-4.0d0*p(i,j  )+p(i,j-1))**2

!   !     tau5 = abs(beta0 - beta2)

!   !     alpha0 = d0 * (1.0d0 + tau5 / (beta0 + eps))
!   !     alpha1 = d1 * (1.0d0 + tau5 / (beta1 + eps))
!   !     alpha2 = d2 * (1.0d0 + tau5 / (beta2 + eps))
!   !     wsum   = alpha0 + alpha1 + alpha2

!   !     w0 = alpha0 / wsum
!   !     w1 = alpha1 / wsum
!   !     w2 = alpha2 / wsum

!   !     pR(i,j) = w0*h0 + w1*h1 + w2*h2

!   !   end do
!   ! end do

! end subroutine weno5_y



! ! subroutine weno_2d_yR(Nx, Ny, p, pRy)
! !         implicit none
! !         integer, intent(in) :: Nx, Ny
! !         real(kind=8), intent(in) :: p(1:Nx,1:Ny)
! !         real(kind=8), intent(out):: pRy(1:Nx,1:Ny)
! !         integer :: i,j

! !         real(kind=8), parameter :: epsW=1.d-3
        
! !         real(kind=8), parameter :: d0=1.d0/10.d0, d1=3.d0/5.d0, d2=3.d0/10.d0
! !         real(kind=8) :: v0,v1,v2, beta0,beta1,beta2
! !         real(kind=8) :: a0,a1,a2,sumA, w0,w1,w2
      
! !         do i=3, Nx-3
! !           do j=3, Ny-3
! !             ! 3rd-order right-biased in y
! !             v0 = (2.d0*p(i,j+3)-7.d0*p(i,j+2)+11.d0*p(i,j+1))/6.d0
! !             v1 = (-1.d0*p(i,j+2)+5.d0*p(i,j+1)+2.d0*p(i,j))/6.d0
! !             v2 = (2.d0*p(i,j+1)+5.d0*p(i,j)-1.d0*p(i,j-1))/6.d0
! !             ! smoothness
! !             beta0=(13.d0/12.d0)*(p(i,j+1)-2.d0*p(i,j+2)+p(i,j+3))**2 + (1.d0/4.d0)*(3.d0*p(i,j+1)-4.d0*p(i,j+2)+p(i,j+3))**2
! !             beta1=(13.d0/12.d0)*(p(i,j)-2.d0*p(i,j+1)+p(i,j+2))**2 + (1.d0/4.d0)*(p(i,j+2)-p(i,j))**2
! !             beta2=(13.d0/12.d0)*(p(i,j-1)-2.d0*p(i,j)+p(i,j+1))**2 + (1.d0/4.d0)*(p(i,j-1)-4.d0*p(i,j)+3.d0*p(i,j+1))**2
! !             a0=d0/((epsW+beta0)**2); a1=d1/((epsW+beta1)**2); a2=d2/((epsW+beta2)**2)
! !             sumA=a0+a1+a2; w0=a0/sumA; w1=a1/sumA; w2=a2/sumA
! !             pRy(i,j)=w0*v0 + w1*v1 + w2*v2
! !           end do
! !         end do
! !       end subroutine weno_2d_yR