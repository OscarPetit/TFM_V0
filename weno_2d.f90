subroutine WENO_2D( Nx, Ny, p, p_left_x, p_right_x, p_left_y, p_right_y)
        implicit none
      
      
        integer, intent(in) :: Nx, Ny
        real*8, intent(in)  :: p(Nx,Ny)      ! cell-centered data (with ghost cells)
        real*8 :: epsW      ! small epsilon 
        real*8, intent(out) :: p_left_x(Nx,Ny), p_right_x(Nx,Ny)
        real*8, intent(out) :: p_left_y(Nx,Ny), p_right_y(Nx,Ny)
      
        integer :: i, j
        real*8 :: v0, v1, v2
        real*8 :: beta0, beta1, beta2
        real*8 :: alpha0, alpha1, alpha2, w0, w1, w2
        real*8 :: u0, u1, u2
      
        real*8 :: d0 ! 3/10  (modificados ahora de acuerdo a github alex pandya)
        real*8 :: d1  ! 3/5
        real*8 :: d2  ! 1/10
      

        epsW =  10.0e-2


        !RECONSTRUYE  pR en i+1/2
        do j = 3, Ny-3
            do i =3, Nx-3

            !i = i+1 (todos los i tienen sumados 1)

            d0 = 1.0/10.0
            d1 = 3.0/5.0
            d2 = 3.0/10.0   
      

              v2 = -1.0/6.0 * p(i-1,j)  +  5.0/6.0 * p(i,j)  +  1.0/3.0 * p(i+1  ,j)
              v1 =  1.0/3.0 * p(i,j)  +  5.0/6.0 * p(i+1  ,j)  -  1.0/6.0 * p(i+2,j)
              v0 = 11.0/6.0 * p(i+1  ,j)  -  7.0/6.0 * p(i+2,j)  +  1.0/3.0 * p(i+3,j)
      

              beta0 = 0.25* ( 3.0*p(i+1  ,j) - 4.0*p(i+2,j) + p(i+3,j) )**2   &
                    + (13.0/12.0) * ( p(i+1  ,j) - 2.0*p(i+2,j) + p(i+3,j) )**2
      
              beta1 = 0.25 * ( p(i+2,j) - p(i,j) )**2                     &
                    + (13.0/12.0) * ( p(i,j) - 2.0*p(i+1,j) + p(i+2,j) )**2
      
              beta2 = 0.25 * ( p(i-1,j) - 4.0*p(i,j) + 3.0*p(i+1,j) )**2     &
                    + (13.0/12.0) * ( p(i-1,j) - 2.0*p(i,j) + p(i+1,j) )**2
      

              alpha0 = d0 / ( (epsW + beta0)**2 )  ! uses beta^0
              alpha1 = d1 / ( (epsW + beta1)**2 )  ! uses beta^1
              alpha2 = d2 / ( (epsW + beta2)**2 )  ! uses beta^2
      
              w0 = alpha0 / ( alpha0 + alpha1 + alpha2 )
              w1 = alpha1 / ( alpha0 + alpha1 + alpha2 )
              w2 = alpha2 / ( alpha0 + alpha1 + alpha2 )
      

              p_right_x(i,j) = w0*v0 + w1*v1 + w2*v2
      
            enddo
      enddo


      !qL en i+1/2
      do j = 3, Ny-2
            do i =3, Nx-2


            u0 =  1.0/3.0*p(i  ,j) +  5.0/6.0*p(i+1,j) - 1.0/6.0*p(i+2,j)
            u1 = -1.0/6.0*p(i-1,j) +  5.0/6.0*p(i  ,j) + 1.0/3.0*p(i+1,j)
            u2 =  1.0/3.0*p(i-2,j) -  7.0/6.0*p(i-1,j) + 11.0/6.0*p(i  ,j)


              d0 = 3.0/10.0
              d1 = 3.0/5.0
              d2 = 1.0/10.0


              beta0 = 0.25* ( 3.0*p(i  ,j) - 4.0*p(i+1,j) + p(i+2,j) )**2   &
                    + (13.0/12.0) * ( p(i  ,j) - 2.0*p(i+1,j) + p(i+2,j) )**2
      
              beta1 = 0.25 * ( p(i+1,j) - p(i-1,j) )**2                     &
                    + (13.0/12.0) * ( p(i-1,j) - 2.0*p(i,j) + p(i+1,j) )**2
      
              beta2 = 0.25 * ( p(i-2,j) - 4.0*p(i-1,j) + 3.0*p(i,j) )**2     &
                    + (13.0/12.0) * ( p(i-2,j) - 2.0*p(i-1,j) + p(i,j) )**2




              alpha0 = d0 / ( (epsW + beta0)**2 )  ! uses beta^0
              alpha1 = d1 / ( (epsW + beta1)**2 )  ! uses beta^1
              alpha2 = d2 / ( (epsW + beta2)**2 )  ! uses beta^2
      
              w0 = alpha0 / ( alpha0 + alpha1 + alpha2 )
              w1 = alpha1 / ( alpha0 + alpha1 + alpha2 )
              w2 = alpha2 / ( alpha0 + alpha1 + alpha2 )
      


      
              ! same w0, w1, w2 as above
              p_left_x(i,j) = w0*u0 + w1*u1 + w2*u2


      
           end do
        end do
      
        !======================================================================
        ! 2) Y-DIRECTION RECONSTRUCTION
        !    Similarly, p_right_y(i,j) = p^+_{i,j+1/2}
        !                p_left_y(i,j)  = p^-_{i,j+1/2}
        !======================================================================



        !reconstruye pR en j+1/2
        do j = 3, Ny-3
            do i =3, Nx-3
      
                  !j=j+1 mismo que en xr

            d0 = 1.0/10.0
            d1 = 3.0/5.0
            d2 = 3.0/10.0 
      
              !----------------------------
              ! 2A) Polynomials v^0,v^1,v^2 (in y-direction)
              !     v^0 => {j-2,j-1,j}, etc.
              !----------------------------
              v2 = -1.0/6.0*p(i,j-1) +  5.0/6.0*p(i,j) +  1.0/3.0*p(i,j+1  )
              v1 =  1.0/3.0*p(i,j) +  5.0/6.0*p(i,j+1  ) -  1.0/6.0*p(i,j+2)
              v0 = 11.0/6.0*p(i,j +1 ) -  7.0/6.0*p(i,j+2) +  1.0/3.0*p(i,j+3)
      
              !----------------------------
              ! 2B) Beta^0 => {j,j+1,j+2}, etc.
              !----------------------------
              beta0 = 0.25 * ( 3.0*p(i,j +1 ) - 4.0*p(i,j+2) + p(i,j+3) )**2   &
                    + (13.0/12.0)*( p(i,j +1 ) - 2.0*p(i,j+2) + p(i,j+3) )**2
      
              beta1 = 0.25 * ( p(i,j+2) - p(i,j) )**2                     &
                    + (13.0/12.0)*( p(i,j) - 2.0*p(i,j +1 ) + p(i,j+2) )**2
      
              beta2 = 0.25 * ( p(i,j-1) - 4.0*p(i,j) + 3.0*p(i,j+1) )**2     &
                    + (13.0/12.0)*( p(i,j-1) - 2.0*p(i,j) + p(i,j +1 ) )**2
      
              alpha0 = d0 / ( (epsW + beta0)**2 )
              alpha1 = d1 / ( (epsW + beta1)**2 )
              alpha2 = d2 / ( (epsW + beta2)**2 )
      
              w0 = alpha0 / ( alpha0 + alpha1 + alpha2 )
              w1 = alpha1 / ( alpha0 + alpha1 + alpha2 )
              w2 = alpha2 / ( alpha0 + alpha1 + alpha2 )
      
              p_right_y(i,j) = w0*v0 + w1*v1 + w2*v2
      
            enddo
      enddo



              !reconstruye pL en j+1/2
      do j = 3, Ny-2
            do i =3, Nx-2

                  u0 =  1.0/3.0*p(i,j  ) +  5.0/6.0*p(i,j+1) - 1.0/6.0*p(i,j+2)
                  u1 = -1.0/6.0*p(i,j-1) +  5.0/6.0*p(i,j  ) + 1.0/3.0*p(i,j+1)
                  u2 =  1.0/3.0*p(i,j-2) -  7.0/6.0*p(i,j-1) + 11.0/6.0*p(i,j  )


              d0 = 3.0/10.0
              d1 = 3.0/5.0
              d2 = 1.0/10.0




            beta0 = 0.25 * ( 3.0*p(i,j  ) - 4.0*p(i,j+1) + p(i,j+2) )**2   &
              + (13.0/12.0)*( p(i,j  ) - 2.0*p(i,j+1) + p(i,j+2) )**2

            beta1 = 0.25 * ( p(i,j+1) - p(i,j-1) )**2                     &
              + (13.0/12.0)*( p(i,j-1) - 2.0*p(i,j  ) + p(i,j+1) )**2

            beta2 = 0.25 * ( p(i,j-2) - 4.0*p(i,j-1) + 3.0*p(i,j) )**2     &
              + (13.0/12.0)*( p(i,j-2) - 2.0*p(i,j-1) + p(i,j  ) )**2


              alpha0 = d0 / ( (epsW + beta0)**2 )  ! uses beta^0
              alpha1 = d1 / ( (epsW + beta1)**2 )  ! uses beta^1
              alpha2 = d2 / ( (epsW + beta2)**2 )  ! uses beta^2
      
              w0 = alpha0 / ( alpha0 + alpha1 + alpha2 )
              w1 = alpha1 / ( alpha0 + alpha1 + alpha2 )
              w2 = alpha2 / ( alpha0 + alpha1 + alpha2 )




      
              p_left_y(i,j)  = w0*u0 + w1*u1 + w2*u2

              
      
           end do
        end do


end subroutine WENO_2D
      