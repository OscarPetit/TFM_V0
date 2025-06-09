subroutine initial_conditions(test, nx, ny, x, y, psi, dens, ux, uy, ut, dtpsi, dxpsi, dypsi, dtux, dxux, dyux, dtuy,dxuy, dyuy, q11, q12, q13, jt)
   implicit none
   integer, intent(in) :: nx, ny,  test
   double precision, intent(inout) :: psi(nx, ny), dtpsi(nx, ny), dxpsi(nx, ny), dypsi(nx, ny)
   double precision, intent(inout) :: ut(nx, ny)
  !  double precision, intent(inout) :: q01(nx, ny), q02(nx, ny), q03(nx, ny)
   double precision, intent(inout) :: jt(nx, ny)
   double precision :: epspf(nx, ny)
   double precision, intent(inout) :: ux(nx, ny), dtux(nx, ny), dxux(nx, ny), dyux(nx, ny)
   double precision, intent(inout) :: uy(nx, ny), dtuy(nx, ny), dxuy(nx, ny), dyuy(nx, ny)
   double precision, intent(inout) :: dens(nx, ny)
   double precision, intent(inout) :: q11(nx, ny), q12(nx, ny), q13(nx, ny)
   integer :: i, j
   real(kind=8) :: delta, L, xval, yval, dn, D, B, theta, lorentz_factor, o, k, sigma, y1, y2, uflow, vx, vy, gamma, de, dvx, dvy,w
   real(kind=8), intent(in) :: x(nx), y(ny)
   


   if (test==0) then


      do j = 1, ny
        do i=1 ,nx
    
              ! Parameters
           delta = 0.05d0      ! Smoothness of the transitions    
           L = 3.0d0 
    
          ! Compute grid positions
          xval = x(i)
          yval = y(j)
        
            if ((sqrt(xval**2+yval**2) <= 0.5d0) .and. (abs(yval)<0.1d0)) then
              B=0.1d0
            else
              B=0.d0
            end if
          
          theta = atan2(yval,xval)
          dn = 0.5d0 - sqrt(xval**2+yval**2)
          ! Compute D(d, delta)
          D = 0.5d0 * (1.0d0 + tanh(dn/delta))
        
          ! Set initial conditions
          dens(i, j) = 0.5d0 * (1.0d0 + D) + B  ! n(x, y)
          !dens(i,j) = 10*exp(-(x(i))**2/(5**2))+0.5
          psi(i, j) = 0.0d0
        
          ux(i, j) = -1.0d0*sqrt(xval**2+yval**2)*sin(theta)*D           ! Initial vx
          uy(i, j) = 1.0d0*sqrt(xval**2+yval**2)*cos(theta)*D           ! Initial vy
          !ux(i,j) = 0.0d0
          !uy(i,j) = 0.0d0
          ! Compute the Lorentz factor
          lorentz_factor = sqrt(1.0d0/(1.0d0-(ux(i,j)**2+uy(i,j)**2)))
    
          ! Compute vx and vy
          
          ux(i, j) = lorentz_factor * ux(i, j)
          uy(i, j) = lorentz_factor * uy(i, j)
    
          ! Recalculate ut
          ut(i, j) = sqrt(1.0d0 + ux(i, j)**2 + uy(i, j)**2)
    
          dtpsi(i,j)=0.0d0 !Inicializamos las derivadas igual a 0 porque entran en q11 q12 y q13
          !dxpsi(i,j)=0.0d0
          !dypsi(i,j)=0.0d0
    
          dtux(i,j)=0.0d0
          !dxux(i,j)=0.0d0
          !dyux(i,j)=0.0d0
    
          dtuy(i,j)=0.0d0
          !dxuy(i,j)=0.0d0
          !dyuy(i,j)=0.0d0
    
          q11(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ut(i,j) - 1.0d0/3.0d0)
          q12(i,j) = (exp(psi(i,j)))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
          q13(i,j) = (exp(psi(i,j)))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
    
          ! q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
          ! q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
          ! q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
    
    
          jt(i,j) = dens(i,j)*ut(i,j)
    
    
        end do
      end do
    
    
    elseif (test==1) then
    
    
      do j = 4, ny-3
        do i=4 ,nx-3
    
    
          o=0.01d0
          k=0.05d0
          sigma=0.2d0
          y1=-0.5d0
          y2=0.5d0
          uflow=1.0d0/4.0d0/sqrt(3.0d0)
    
           xval = x(i)
           yval = y(j)
    
           ! n(x,y)
           dens(i,j) = 1.0d0 + 0.5d0 * (  &
         &    tanh((yval - y1)/k) - tanh((yval - y2)/k) )
    
           ! vx(x,y)
           vx = uflow * ( tanh((yval - y1)/k) - tanh((yval - y2)/k) - 1.0d0 )
    
           ! vy(x,y)
           vy = o * sin(2.0d0*acos(-1.0d0)*xval) * (  &
         &       exp(-((yval - y1)/sigma)**2)  &
         &     + exp(-((yval - y2)/sigma)**2) )
    
         lorentz_factor=sqrt(1.0d0/(1.0d0 - (vx**2 + vy**2)))
    
           ux(i,j) = vx*lorentz_factor
           uy(i,j) = vy*lorentz_factor
    
           ! You can keep epspf = 1.0d0 and psi = 0.0d0, as in your example:
           epspf(i,j) = 1.0d0
           psi(i,j)   = 0.0d0
    
           ! -- Relativistic corrections: compute Lorentz factor --
           ut(i,j) = sqrt(1.0d0 + ux(i,j)**2 + uy(i,j)**2)
    
           ! -- Zero out derivatives, as in your rotor code --
           dtpsi(i,j) = 0.0d0; dxpsi(i,j) = 0.0d0; dypsi(i,j) = 0.0d0
           dtux(i,j)  = 0.0d0; dxux(i,j)  = 0.0d0; dyux(i,j)  = 0.0d0
           dtuy(i,j)  = 0.0d0; dxuy(i,j)  = 0.0d0; dyuy(i,j)  = 0.0d0
    
           ! -- Define the q-tensors (your viscous variables) similarly --
           q11(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ut(i,j) - 1.0d0/3.0d0)
           q12(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
           q13(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
    
          !  q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
          !  q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
          !  q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
    
           ! -- Charge/current 4-vector, e.g. --
           jt(i,j) = dens(i,j)*ut(i,j)
    
        end do
      end do
    
    
    elseif(test==2) then
    
    
      ! do j = 4, ny-3
      !   do i=4, nx-3
    
    
      !     ! Parameters
      !    delta = 10.0d0       ! 10 by default Smoothness of the transitions
      !    gamma = 10.0d0       ! Controls the squareness of each quadrant
      !    L = 200.0d0            ! Size parameter
    
      !     ! Compute grid positions
      !     xval = x(i)
      !     yval = y(j)
        
      !     ! Compute dn, de, dvx, and dvy
      !     dn = L - ((xval + L)**gamma + (yval + L)**gamma)**(1.0d0 / gamma)
      !     de = L - ((xval - L)**gamma + (yval - L)**gamma)**(1.0d0 / gamma)
      !     dvx = L - ((xval + L)**gamma + (yval - L)**gamma)**(1.0d0 / gamma)
      !     dvy = L - ((xval - L)**gamma + (yval + L)**gamma)**(1.0d0 / gamma)
        
      !     ! Compute D(d, delta)
      !      D = 0.5d0 * (1.0d0 + TANH(dn / delta))
        
      !     ! Set initial conditions
      !     dens(i, j) = 0.4d0 * (1.0d0 + TANH(dn / delta)) + 0.1d0  ! n(x, y)
      !     !dens(i,j) = 10*exp(-(x(i))**2/(5**2))+0.5
      !     psi(i, j) = log(3.0d0 - 2.97d0 * 0.5d0 * (1.0d0 + TANH(de / delta)))  ! epsilon(x, y)
        
      !     ux(i, j) = 0.95d0 * 0.5d0 * (1.0d0 + TANH(dvx / delta))           ! Initial ux
      !     uy(i, j) = 0.95d0 * 0.5d0 * (1.0d0 + TANH(dvy / delta))          ! Initial uy
      !     !ux(i,j) = 0.0d0
      !     !uy(i,j) = 0.0d0
      !     ! Compute the Lorentz factor
      !     lorentz_factor = sqrt(1.0d0 + ux(i, j)**2 + uy(i, j)**2)
    
      !     ! Compute vx and vy
      !     ux(i, j) = lorentz_factor * ux(i, j)
      !     uy(i, j) = lorentz_factor * uy(i, j)
    
      !      ! -- Relativistic corrections: compute Lorentz factor --
      !      ut(i,j) = sqrt(1.0d0 + ux(i,j)**2 + uy(i,j)**2)
    
      !      ! -- Zero out derivatives, as in your rotor code --
      !      dtpsi(i,j) = 0.0d0; dxpsi(i,j) = 0.0d0; dypsi(i,j) = 0.0d0
      !      dtux(i,j)  = 0.0d0; dxux(i,j)  = 0.0d0; dyux(i,j)  = 0.0d0
      !      dtuy(i,j)  = 0.0d0; dxuy(i,j)  = 0.0d0; dyuy(i,j)  = 0.0d0
    
      !      ! -- Define the q-tensors (your viscous variables) similarly --
      !      q11(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ut(i,j) - 1.0d0/3.0d0)
      !      q12(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
      !      q13(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
    
      !     !  q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
      !     !  q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
      !     !  q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
    
      !      ! -- Charge/current 4-vector, e.g. --
      !      jt(i,j) = dens(i,j)*ut(i,j)
    
      !   end do
      ! end do
    
    
      elseif(test==3) then
    
    
        do j = 1, ny
          do i =1, nx
      
      
            ! Parameters
           delta = 10.0d0! 10 by default Smoothness of the transitions
           gamma = 10.0d0       ! Controls the squareness of each quadrant
           L = 200.0d0            ! Size parameter
      
            ! Compute grid positions
            xval = x(i)
            yval = y(j)
          
            ! Compute dn, de, dvx, and dvy
            dn = L - ((xval + L)**gamma + (yval + L)**gamma)**(1.0d0 / gamma)
            de = L - ((xval - L)**gamma + (yval - L)**gamma)**(1.0d0 / gamma)
            dvx = L - ((xval + L)**gamma + (yval - L)**gamma)**(1.0d0 / gamma)
            dvy = L - ((xval - L)**gamma + (yval + L)**gamma)**(1.0d0 / gamma)
          
            ! Compute D(d, delta)
             D = 0.5d0 * (1.0d0 + TANH(dn / delta))
          
            ! Set initial conditions
            dens(i, j) = 0.4d0 * 0.5d0*(1.0d0 + TANH(dn / delta)) + 0.1d0  ! n(x, y)
            !dens(i,j) = 10*exp(-(x(i))**2/(5**2))+0.5
            psi(i, j) = log(3.0d0 - 2.97d0 * 0.5d0 * (1.0d0 + TANH(de / delta)))  ! epsilon(x, y)
          
            ux(i, j) = 0.97d0 * 0.5d0 * (1.0d0 + TANH(dvx / delta))           ! Initial ux
            uy(i, j) = 0.97d0 * 0.5d0 * (1.0d0 + TANH(dvy / delta))          ! Initial uy
            !ux(i,j) = 0.0d0
            !uy(i,j) = 0.0d0
            ! Compute the Lorentz factor
            lorentz_factor = sqrt(1.0d0/(1.0d0-(ux(i,j)**2+uy(i,j)**2)))
      
            ! Compute vx and vy
            ux(i, j) = lorentz_factor * ux(i, j)
            uy(i, j) = lorentz_factor * uy(i, j)
      
             ! -- Relativistic corrections: compute Lorentz factor --
             ut(i,j) = sqrt(1.0d0 + ux(i,j)**2 + uy(i,j)**2)
      
             ! -- Zero out derivatives, as in your rotor code --
             dtpsi(i,j) = 0.0d0; dxpsi(i,j) = 0.0d0; dypsi(i,j) = 0.0d0
             dtux(i,j)  = 0.0d0; dxux(i,j)  = 0.0d0; dyux(i,j)  = 0.0d0
             dtuy(i,j)  = 0.0d0; dxuy(i,j)  = 0.0d0; dyuy(i,j)  = 0.0d0
      
             ! -- Define the q-tensors (your viscous variables) similarly --
             q11(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ut(i,j) - 1.0d0/3.0d0)
             q12(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
             q13(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
      
            !  q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
            !  q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
            !  q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
      
             ! -- Charge/current 4-vector, e.g. --
             jt(i,j) = dens(i,j)*ut(i,j)
      
          end do
        end do
    

      elseif (test==4) then

        do j = 1, ny
          do i = 1, nx

            xval = x(i)
            yval = y(j)


            delta = 1.0d0/10.0d0

            dens(i,j) = 0.0d0
            psi(i,j) = log(exp(-xval**2/25.0d0**2) + delta)
            ux(i,j) = 0.0d0
            uy(i,j) = 0.0d0
            ut(i,j) = 1.0d0

            dtpsi(i,j) = 0.0d0
            dtux(i,j) = 0.0d0
            dtuy(i,j) = 0.0d0

            dxpsi(i,j) = 0.0d0
            dxux(i,j) = 0.0d0
            dxuy(i,j) = 0.0d0

            dypsi(i,j) = 0.0d0
            dyux(i,j) = 0.0d0
            dyuy(i,j) = 0.0d0

             ! -- Define the q-tensors (your viscous variables) similarly --
            q11(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ut(i,j) - 1.0d0/3.0d0)
            q12(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
            q13(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
     
            ! q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
            ! q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
            ! q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
     
            ! -- Charge/current 4-vector, e.g. --
            jt(i,j) = dens(i,j)*ut(i,j)

          end do
        end do
    

        elseif (test==5) then

          do j = 1, ny
            do i = 1, nx
  
              xval = x(i)
              yval = y(j)
  
  
              delta = 1.0d0/10.0d0
  
              dens(i,j) = 0.0d0
              psi(i,j) = log(exp(-yval**2/25.0d0**2) + delta)
              ux(i,j) = 0.0d0
              uy(i,j) = 0.0d0
              ut(i,j) = 1.0d0
  
              dtpsi(i,j) = 0.0d0
              dtux(i,j) = 0.0d0
              dtuy(i,j) = 0.0d0
  
              dxpsi(i,j) = 0.0d0
              dxux(i,j) = 0.0d0
              dxuy(i,j) = 0.0d0
  
              dypsi(i,j) = 0.0d0
              dyux(i,j) = 0.0d0
              dyuy(i,j) = 0.0d0
  
               ! -- Define the q-tensors (your viscous variables) similarly --
              q11(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ut(i,j) - 1.0d0/3.0d0)
              q12(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
              q13(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
       
              ! q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
              ! q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
              ! q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
       
              ! -- Charge/current 4-vector, e.g. --
              jt(i,j) = dens(i,j)*ut(i,j)
  
            end do
          end do

          elseif (test==6) then

            do j = 1, ny
              do i = 1, nx
    
                xval = x(i)
                yval = y(j)
    
               dens(i,j) = 0.0d0

                !w = 5.0d0
                w = 2.0d0

              if (xval <= 0.0d0) then
              psi(i,j) = log(1.0d0)
              else
                psi(i,j) = log(0.1d0)
              end if

              psi(i,j) = log( 0.5d0*(1.0d0 + 0.1d0) + 0.5d0*(1.0d0 - 0.1d0) * tanh(-xval / w) )
               
                ux(i,j) = 0.0d0
                uy(i,j) = 0.0d0
                ut(i,j) = 1.0d0
    
                dtpsi(i,j) = 0.0d0
                dtux(i,j) = 0.0d0
                dtuy(i,j) = 0.0d0
    
                dxpsi(i,j) = 0.0d0
                dxux(i,j) = 0.0d0
                dxuy(i,j) = 0.0d0
    
                dypsi(i,j) = 0.0d0
                dyux(i,j) = 0.0d0
                dyuy(i,j) = 0.0d0

    
                 ! -- Define the q-tensors (your viscous variables) similarly --
                q11(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ut(i,j) - 1.0d0/3.0d0)
                q12(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
                q13(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
         
                ! q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
                ! q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
                ! q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
         
                ! -- Charge/current 4-vector, e.g. --
                jt(i,j) = dens(i,j)*ut(i,j)
    
              end do
            end do

            elseif (test==7) then

              do j = 1, ny
                do i = 1, nx
      
                  xval = x(i)
                  yval = y(j)
      
      
      
                  dens(i,j) = 1.0d0
                  psi(i,j) = log(1.0d0)

                 
                  ux(i,j) = 0.0d0
                  uy(i,j) = 0.0d0
                  ut(i,j) = 1.0d0
      
                  dtpsi(i,j) = 0.0d0
  
      
                  dtux(i,j) = 0.0d0
  
      
                  dtuy(i,j) = 0.0d0
  
      
                   ! -- Define the q-tensors (your viscous variables) similarly --
                  q11(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ut(i,j) - 1.0d0/3.0d0)
                  q12(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
                  q13(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
           
                  ! q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
                  ! q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
                  ! q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
           
                  ! -- Charge/current 4-vector, e.g. --
                  jt(i,j) = dens(i,j)*ut(i,j)
      
                end do
              end do
    
              elseif (test==8) then
                ! Smooth sinusoidal test for WENO convergence on [0,1]×[0,1]
                do j = 1, ny
                  do i = 1, nx
        
                    xval = x(i)
                    yval = y(j)
        
    
                    !w = 5.0d0
                   ! w = 6.0d0
    
                  if (xval <= 0.0d0) then
                  dens(i,j) = 1.0d0
                  psi(i,j) = log(0.1d0)
                  ux(i,j) = 0.9d0
                  uy(i,j) = 0.0d0
                  ut(i,j) = 1.0d0/sqrt(1.0d0 - ux(i,j)**2)
                  ux(i,j) = ux(i,j)*ut(i,j)
                  else
                    dens(i,j) = 1.0d0
                    psi(i,j) = log(0.5d0)
                    ux(i,j) = 0.0d0
                    uy(i,j) = 0.0d0
                    ut(i,j) = 1.0d0
                    ux(i,j) = ux(i,j)*ut(i,j)
                  end if
    
                  !psi(i,j) = log( 0.5d0*(1.0d0 + 0.1d0) + 0.5d0*(1.0d0 - 0.1d0) * tanh(-xval / w) )
                   
        
                    dtpsi(i,j) = 0.0d0
                    dtux(i,j) = 0.0d0
                    dtuy(i,j) = 0.0d0
        
                    dxpsi(i,j) = 0.0d0
                    dxux(i,j) = 0.0d0
                    dxuy(i,j) = 0.0d0
        
                    dypsi(i,j) = 0.0d0
                    dyux(i,j) = 0.0d0
                    dyuy(i,j) = 0.0d0
    
        
                     ! -- Define the q-tensors (your viscous variables) similarly --
                    q11(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ut(i,j) - 1.0d0/3.0d0)
                    q12(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
                    q13(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
             
                    ! q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
                    ! q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
                    ! q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))
             
                    ! -- Charge/current 4-vector, e.g. --
                    jt(i,j) = dens(i,j)*ut(i,j)
                  end do
                end do
    
    
    
    endif
    
 
end subroutine initial_conditions
 