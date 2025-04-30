subroutine initial_conditions(test, nx, ny, x, y, psi, dens, ux, uy, ut, dtpsi, dxpsi, dypsi, dtux, dxux, dyux, dtuy,dxuy, dyuy, q11, q12, q13, jt, q01, q02, q03)
   implicit none
   integer, intent(in) :: nx, ny,  test
   double precision, intent(inout) :: psi(nx, ny), dtpsi(nx, ny), dxpsi(nx, ny), dypsi(nx, ny)
   double precision, intent(inout) :: ut(nx, ny)
   double precision, intent(inout) :: q01(nx, ny), q02(nx, ny), q03(nx, ny)
   double precision, intent(inout) :: jt(nx, ny)
   double precision :: epspf(nx, ny)
   double precision, intent(inout) :: ux(nx, ny), dtux(nx, ny), dxux(nx, ny), dyux(nx, ny)
   double precision, intent(inout) :: uy(nx, ny), dtuy(nx, ny), dxuy(nx, ny), dyuy(nx, ny)
   double precision, intent(inout) :: dens(nx, ny)
   double precision, intent(inout) :: q11(nx, ny), q12(nx, ny), q13(nx, ny)
   integer :: i, j
   real*8 :: delta, L, xval, yval, dn, D, B, theta, lorentz_factor, o, k, sigma, y1, y2, uflow, vx, vy, gamma, de, dvx, dvy
   real*8, intent(in) :: x(nx), y(ny)
   


   if (test==0) then




  
    
    
     ! Example: Assign some initial values or conditions (optional) 
      do j = 4, ny-3
        do i=4 ,nx-3
    
              ! Parameters
           delta = 0.05      ! Smoothness of the transitions    
           L = 3.0 
    
          ! Compute grid positions
          xval = x(i)
          yval = y(j)
        
            if ((sqrt(xval**2+yval**2) <= 0.5) .and. (abs(yval)<0.1)) then
              B=0.1
            else
              B=0
            end if
          
          theta = atan2(yval,xval)
          dn = 0.5 - sqrt(xval**2+yval**2)
          ! Compute D(d, delta)
          D = 0.5 * (1.0 + tanh(dn/delta))
        
          ! Set initial conditions
          dens(i, j) = 0.5 * (1.0 + D) + B  ! n(x, y)
          !dens(i,j) = 10*exp(-(x(i))**2/(5**2))+0.5
          psi(i, j) = 0.0
        
          ux(i, j) = -1.0*sqrt(xval**2+yval**2)*sin(theta)*D           ! Initial vx
          uy(i, j) = 1.0*sqrt(xval**2+yval**2)*cos(theta)*D           ! Initial vy
          !ux(i,j) = 0.0
          !uy(i,j) = 0.0
          ! Compute the Lorentz factor
          lorentz_factor = sqrt(1.0/(1.0-(ux(i,j)**2+uy(i,j)**2)))
    
          ! Compute vx and vy
          
          ux(i, j) = lorentz_factor * ux(i, j)
          uy(i, j) = lorentz_factor * uy(i, j)
    
          ! Recalculate ut
          ut(i, j) = sqrt(1.0 + ux(i, j)**2 + uy(i, j)**2)
    
          dtpsi(i,j)=0.0 !Inicializamos las derivadas igual a 0 porque entran en q11 q12 y q13
          !dxpsi(i,j)=0.0
          !dypsi(i,j)=0.0
    
          dtux(i,j)=0.0
          !dxux(i,j)=0.0
          !dyux(i,j)=0.0
    
          dtuy(i,j)=0.0
          !dxuy(i,j)=0.0
          !dyuy(i,j)=0.0
    
          q11(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ut(i,j) - 1.0/3.0)
          q12(i,j) = (exp(psi(i,j)))*(4.0/3.0*ut(i,j)*ux(i,j))
          q13(i,j) = (exp(psi(i,j)))*(4.0/3.0*ut(i,j)*uy(i,j))
    
          q01(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)**2 - 1.0/3.0)
          q02(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
          q03(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
    
    
          jt(i,j) = dens(i,j)*ut(i,j)
    
    
        end do
      end do
    
    
    elseif (test==1) then
    
    
      do j = 4, ny-3
        do i=4 ,nx-3
    
    
          o=0.01
          k=0.05
          sigma=0.2
          y1=-0.5
          y2=0.5
          uflow=1.0/4.0/sqrt(3.0)
    
           xval = x(i)
           yval = y(j)
    
           ! n(x,y)
           dens(i,j) = 1.0 + 0.5 * (  &
         &    tanh((yval - y1)/k) - tanh((yval - y2)/k) )
    
           ! vx(x,y)
           vx = uflow * ( tanh((yval - y1)/k) - tanh((yval - y2)/k) - 1.0 )
    
           ! vy(x,y)
           vy = o * sin(2.0*acos(-1.0)*xval) * (  &
         &       exp(-((yval - y1)/sigma)**2)  &
         &     + exp(-((yval - y2)/sigma)**2) )
    
         lorentz_factor=sqrt(1.0/(1.0 - (vx**2 + vy**2)))
    
           ux(i,j) = vx!*lorentz_factor
           uy(i,j) = vy!*lorentz_factor
    
           ! You can keep epspf = 1.0 and psi = 0.0, as in your example:
           epspf(i,j) = 1.0
           psi(i,j)   = 0.0
    
           ! -- Relativistic corrections: compute Lorentz factor --
           ut(i,j) = sqrt(1.0 + ux(i,j)**2 + uy(i,j)**2)
    
           ! -- Zero out derivatives, as in your rotor code --
           dtpsi(i,j) = 0.0; dxpsi(i,j) = 0.0; dypsi(i,j) = 0.0
           dtux(i,j)  = 0.0; dxux(i,j)  = 0.0; dyux(i,j)  = 0.0
           dtuy(i,j)  = 0.0; dxuy(i,j)  = 0.0; dyuy(i,j)  = 0.0
    
           ! -- Define the q-tensors (your viscous variables) similarly --
           q11(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ut(i,j) - 1.0/3.0)
           q12(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
           q13(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
    
           q01(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)**2 - 1.0/3.0)
           q02(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
           q03(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
    
           ! -- Charge/current 4-vector, e.g. --
           jt(i,j) = dens(i,j)*ut(i,j)
    
        end do
      end do
    
    
    elseif(test==2) then
    
    
      do j = 4, ny-3
        do i=4, nx-3
    
    
          ! Parameters
         delta = 10.0       ! 10 by default Smoothness of the transitions
         gamma = 10.0       ! Controls the squareness of each quadrant
         L = 200.0            ! Size parameter
    
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
          dens(i, j) = 0.4d0 * (1.0d0 + TANH(dn / delta)) + 0.1d0  ! n(x, y)
          !dens(i,j) = 10*exp(-(x(i))**2/(5**2))+0.5
          psi(i, j) = log(3.0d0 - 2.97d0 * 0.5d0 * (1.0d0 + TANH(de / delta)))  ! epsilon(x, y)
        
          ux(i, j) = 0.95d0 * 0.5d0 * (1.0d0 + TANH(dvx / delta))           ! Initial ux
          uy(i, j) = 0.95d0 * 0.5d0 * (1.0d0 + TANH(dvy / delta))          ! Initial uy
          !ux(i,j) = 0.0
          !uy(i,j) = 0.0
          ! Compute the Lorentz factor
          lorentz_factor = sqrt(1.0d0 + ux(i, j)**2 + uy(i, j)**2)
    
          ! Compute vx and vy
          ux(i, j) = lorentz_factor * ux(i, j)
          uy(i, j) = lorentz_factor * uy(i, j)
    
           ! -- Relativistic corrections: compute Lorentz factor --
           ut(i,j) = sqrt(1.0 + ux(i,j)**2 + uy(i,j)**2)
    
           ! -- Zero out derivatives, as in your rotor code --
           dtpsi(i,j) = 0.0; dxpsi(i,j) = 0.0; dypsi(i,j) = 0.0
           dtux(i,j)  = 0.0; dxux(i,j)  = 0.0; dyux(i,j)  = 0.0
           dtuy(i,j)  = 0.0; dxuy(i,j)  = 0.0; dyuy(i,j)  = 0.0
    
           ! -- Define the q-tensors (your viscous variables) similarly --
           q11(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ut(i,j) - 1.0/3.0)
           q12(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
           q13(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
    
           q01(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)**2 - 1.0/3.0)
           q02(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
           q03(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
    
           ! -- Charge/current 4-vector, e.g. --
           jt(i,j) = dens(i,j)*ut(i,j)
    
        end do
      end do
    
    
      elseif(test==3) then
    
    
        do j = 1, ny
          do i = 1, nx
      
      
            ! Parameters
           delta = 10.0! 10 by default Smoothness of the transitions
           gamma = 10.0       ! Controls the squareness of each quadrant
           L = 200.0            ! Size parameter
      
            ! Compute grid positions
            xval = x(i)
            yval = y(j)
          
            ! Compute dn, de, dvx, and dvy
            dn = L - ((xval + L)**gamma + (yval + L)**gamma)**(1.0 / gamma)
            de = L - ((xval - L)**gamma + (yval - L)**gamma)**(1.0 / gamma)
            dvx = L - ((xval + L)**gamma + (yval - L)**gamma)**(1.0 / gamma)
            dvy = L - ((xval - L)**gamma + (yval + L)**gamma)**(1.0 / gamma)
          
            ! Compute D(d, delta)
             D = 0.5d0 * (1.0d0 + TANH(dn / delta))
          
            ! Set initial conditions
            dens(i, j) = 0.4d0 * 0.5*(1.0d0 + TANH(dn / delta)) + 0.1d0  ! n(x, y)
            !dens(i,j) = 10*exp(-(x(i))**2/(5**2))+0.5
            psi(i, j) = log(3.0d0 - 2.97d0 * 0.5d0 * (1.0d0 + TANH(de / delta)))  ! epsilon(x, y)
          
            ux(i, j) = 0.95d0 * 0.5d0 * (1.0d0 + TANH(dvx / delta))           ! Initial ux
            uy(i, j) = 0.95d0 * 0.5d0 * (1.0d0 + TANH(dvy / delta))          ! Initial uy
            !ux(i,j) = 0.0
            !uy(i,j) = 0.0
            ! Compute the Lorentz factor
            lorentz_factor = sqrt(1.0/(1.0-(ux(i,j)**2+uy(i,j)**2)))
      
            ! Compute vx and vy
            ux(i, j) = lorentz_factor * ux(i, j)
            uy(i, j) = lorentz_factor * uy(i, j)
      
             ! -- Relativistic corrections: compute Lorentz factor --
             ut(i,j) = sqrt(1.0 + ux(i,j)**2 + uy(i,j)**2)
      
             ! -- Zero out derivatives, as in your rotor code --
             dtpsi(i,j) = 0.0; dxpsi(i,j) = 0.0; dypsi(i,j) = 0.0
             dtux(i,j)  = 0.0; dxux(i,j)  = 0.0; dyux(i,j)  = 0.0
             dtuy(i,j)  = 0.0; dxuy(i,j)  = 0.0; dyuy(i,j)  = 0.0
      
             ! -- Define the q-tensors (your viscous variables) similarly --
             q11(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ut(i,j) - 1.0/3.0)
             q12(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
             q13(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
      
             q01(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)**2 - 1.0/3.0)
             q02(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
             q03(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
      
             ! -- Charge/current 4-vector, e.g. --
             jt(i,j) = dens(i,j)*ut(i,j)
      
          end do
        end do
    

      elseif (test==4) then

        do j = 1, ny
          do i = 1, nx

            xval = x(i)
            yval = y(j)


            delta = 1.0/10.0

            dens(i,j) = 0.0
            psi(i,j) = log(exp(-xval**2/25.0**2) + delta)
            ux(i,j) = 0.0
            uy(i,j) = 0.0
            ut(i,j) = 1.0

            dtpsi(i,j) = 0.0
            dtux(i,j) = 0.0
            dtuy(i,j) = 0.0

            dxpsi(i,j) = 0.0
            dxux(i,j) = 0.0
            dxuy(i,j) = 0.0

            dypsi(i,j) = 0.0
            dyux(i,j) = 0.0
            dyuy(i,j) = 0.0

             ! -- Define the q-tensors (your viscous variables) similarly --
            q11(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ut(i,j) - 1.0/3.0)
            q12(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
            q13(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
     
            q01(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)**2 - 1.0/3.0)
            q02(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
            q03(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
     
            ! -- Charge/current 4-vector, e.g. --
            jt(i,j) = dens(i,j)*ut(i,j)

          end do
        end do
    

        elseif (test==5) then

          do j = 1, ny
            do i = 1, nx
  
              xval = x(i)
              yval = y(j)
  
  
              delta = 1.0/10.0
  
              dens(i,j) = 0.0
              psi(i,j) = log(exp(-yval**2/25.0**2) + delta)
              ux(i,j) = 0.0
              uy(i,j) = 0.0
              ut(i,j) = 1.0
  
              dtpsi(i,j) = 0.0
              dtux(i,j) = 0.0
              dtuy(i,j) = 0.0
  
              dxpsi(i,j) = 0.0
              dxux(i,j) = 0.0
              dxuy(i,j) = 0.0
  
              dypsi(i,j) = 0.0
              dyux(i,j) = 0.0
              dyuy(i,j) = 0.0
  
               ! -- Define the q-tensors (your viscous variables) similarly --
              q11(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ut(i,j) - 1.0/3.0)
              q12(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
              q13(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
       
              q01(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)**2 - 1.0/3.0)
              q02(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
              q03(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
       
              ! -- Charge/current 4-vector, e.g. --
              jt(i,j) = dens(i,j)*ut(i,j)
  
            end do
          end do

          elseif (test==6) then

            do j = 1, ny
              do i = 1, nx
    
                xval = x(i)
                yval = y(j)
    
    
    
                dens(i,j) = 0.0
                if (xval<=0) then
                  psi(i,j) = log(1.0)
                else
                  psi(i,j) = log(0.8)
                endif
               
                ux(i,j) = 0.0
                uy(i,j) = 0.0
                ut(i,j) = 1.0
    
                dtpsi(i,j) = 0.0

    
                dtux(i,j) = 0.0

    
                dtuy(i,j) = 0.0

    
                 ! -- Define the q-tensors (your viscous variables) similarly --
                q11(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ut(i,j) - 1.0/3.0)
                q12(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
                q13(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
         
                q01(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)**2 - 1.0/3.0)
                q02(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*ux(i,j))
                q03(i,j) = exp(psi(i,j))*(4.0/3.0*ut(i,j)*uy(i,j))
         
                ! -- Charge/current 4-vector, e.g. --
                jt(i,j) = dens(i,j)*ut(i,j)
    
              end do
            end do
    
    
    
    
    endif
    
 
end subroutine initial_conditions
 