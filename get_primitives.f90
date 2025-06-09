subroutine get_primitives(eta, deltaeta, psi, dtpsi, dxpsi, dypsi, ux, dtux, dxux, dyux, uy, dtuy, dxuy, dyuy, ut, q11, q12, q13, dens, jt, nx, ny )
   implicit none
   integer, intent(in) :: nx, ny
   real(kind=8), intent(inout) :: psi(nx, ny), dtpsi(nx, ny), dxpsi(nx, ny), dypsi(nx, ny)
   real(kind=8), intent(inout) :: ut(nx, ny)
   real(kind=8) :: q01(nx, ny), q02(nx, ny), q03(nx, ny)
   real(kind=8), intent(in) :: jt(nx, ny)
   real(kind=8), intent(inout) :: ux(nx, ny), dtux(nx, ny), dxux(nx, ny), dyux(nx, ny)
   real(kind=8), intent(inout) :: uy(nx, ny), dtuy(nx, ny), dxuy(nx, ny), dyuy(nx, ny)
   real(kind=8), intent(in) :: eta, deltaeta
   real(kind=8), intent(in) :: q11(nx, ny), q12(nx, ny), q13(nx, ny)
   integer :: i, j,t,pf
   real(kind=8), intent(out) :: dens(nx, ny)

   real(kind=8) :: a11(nx, ny), a12(nx, ny), a13(nx, ny), a21(nx, ny), a22(nx, ny), a23(nx, ny), a31(nx, ny), a32(nx, ny), a33(nx, ny), det
   real(kind=8) :: b1(nx,ny), b2(nx,ny), b3(nx,ny)
   real(kind=8) :: inv11(nx, ny), inv12(nx, ny), inv13(nx, ny), inv21(nx, ny), inv22(nx, ny), inv23(nx, ny), inv31(nx, ny), inv32(nx, ny), inv33(nx, ny)
   real(kind=8) :: q11hat(nx, ny), q12hat(nx, ny), q13hat(nx, ny), eps_threshold
   real(kind=8) :: epspf(nx, ny), psipf(nx, ny), absvpf(nx, ny), utpf(nx, ny), uxpf(nx, ny), uypf(nx, ny)
   real(kind=8) :: TttPF(nx, ny), TtxPF(nx, ny), TtyPF(nx, ny), TxxPF(nx, ny), TxyPF(nx, ny), TyyPF(nx, ny)

   
  
 !print *, "Psi problema", psi(270,170)
 
  pf = 1  !0 perfect fluid, 1 non-perfect fluid


   if (pf == 0) then

    do j = 4, ny-3
      do i = 4, nx-3


      q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
      q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
      q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))


      epspf(i,j) = -q01(i,j) + sqrt(6.0d0*q01(i,j)**2 + 3.0d0*(q01(i,j)**2-q02(i,j)**2-q03(i,j)**2))
      psipf(i,j) = log(epspf(i,j))
      absvpf(i,j) = (sqrt(q02(i,j)**2 + q03(i,j)**2))/(q01(i,j) + 3.0d0*epspf(i,j))
      utpf(i,j) = 1.0d0/sqrt(1.0d0 - absvpf(i,j)**2)
      uxpf(i,j) = (3.0d0*utpf(i,j)*q02(i,j))/(3.0d0*q01(i,j) + epspf(i,j))
      uypf(i,j) = (3.0d0*utpf(i,j)*q03(i,j))/(3.0d0*q01(i,j) + epspf(i,j))

      psi(i,j) = psipf(i,j)
      ux(i,j) = uxpf(i,j)
      uy(i,j) = uypf(i,j)
      ut(i,j) = utpf(i,j)

     
      ut(i,j) = sqrt(1.0d0 + ux(i,j)**2 + uy(i,j)**2)
      dens(i,j) = jt(i,j)/ut(i,j)


      enddo
    enddo

  elseif(pf == 1 ) then 


    do j = 4, ny-3
      do i = 4, nx-3

        eps_threshold=0.0
      

      a11(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/16.0d0*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)*3.0d0*ut(i,j) + 2.0d0*(ux(i,j)**2)*25.0d0/28.0d0*ut(i,j) + 2.0d0*(uy(i,j)**2)*25.0d0/28.0d0*ut(i,j)) 
      a12(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/4.0d0*ux(i,j)/ut(i,j)*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0) + 2.0d0*ux(i,j)*25.0d0/7.0d0*ut(i,j))
      a13(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/4.0d0*uy(i,j)/ut(i,j)*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0) + 2.0d0*uy(i,j)*25.0d0/7.0d0*ut(i,j)) 

      a21(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/16.0d0*4.0d0/3.0d0*ux(i,j)*ut(i,j)*3.0d0*ut(i,j) + (ut(i,j) + (ux(i,j)**2)/ut(i,j))*25.0d0/28.0d0*ux(i,j)*ut(i,j) + 25.0d0/28.0d0*ux(i,j)*uy(i,j)**2 )
      a22(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/16.0d0*4.0d0*4.0d0/3.0d0*ux(i,j)**2 + (ut(i,j)+(ux(i,j)**2)/ut(i,j))*25.0d0/7.0d0*ut(i,j) - 4.0d0/3.0d0/(ut(i,j)**2)*ux(i,j)**4 - (ux(i,j)**2)/(ut(i,j)**2)*2.0d0/3.0d0*(3.0d0*uy(i,j)**2 + 2.0d0) - (uy(i,j)**2)/(ut(i,j)**2)*1.0d0/3.0d0*(ux(i,j)**2 + 3.0d0*uy(i,j)**2 + 3.0d0)) 
      a23(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/16.0d0*4.0d0*4.0d0/3.0d0*ux(i,j)*uy(i,j) + 25.0d0/7.0d0*ux(i,j)*uy(i,j) + (ux(i,j)**3)/(ut(i,j)**2)*2.0d0/3.0d0*uy(i,j) + (ux(i,j))/(ut(i,j)**2)*2.0d0/3.0d0*uy(i,j) - (uy(i,j))/(3.0d0*ut(i,j)**2)*ux(i,j)*(3.0d0*ux(i,j)**2 + uy(i,j)**2 + 3.0d0) )

      a31(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/4.0d0*uy(i,j)*ut(i,j)**2 + (ut(i,j) + (uy(i,j)**2)/ut(i,j))*25.0d0/28.0d0*uy(i,j)*ut(i,j) + 25.0d0/28.0d0*uy(i,j)*ux(i,j)**2)
      a32(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/4.0d0*ux(i,j)*4.0d0/3.0d0*uy(i,j) + 25.0d0/7.0d0*ux(i,j)*uy(i,j) - uy(i,j)*ux(i,j)/(3.0d0*ut(i,j)**2)*(ux(i,j)**2 + 3.0d0*uy(i,j)**2 + 3.0d0) + 2.0d0/3.0d0*uy(i,j)*ux(i,j)/(ut(i,j)**2)*(uy(i,j)**2 + 1.0d0))
      a33(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/4.0d0*4.0d0/3.0d0*uy(i,j)**2 + 25.0d0/7.0d0*(ut(i,j) + (uy(i,j)**2)/ut(i,j))*ut(i,j) - (ux(i,j)**2)/(3.0d0*ut(i,j)**2)*(3.0d0*ux(i,j)**2 + uy(i,j)**2 + 3.0d0) - 2.0d0/(ut(i,j)**2)*(ux(i,j)**2)*(uy(i,j)**2) - 4.0d0/3.0d0*uy(i,j)/(ut(i,j)**2)*(uy(i,j)**3 + uy(i,j))) 



      ! Definicion de vector b
      b1(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/16.0d0*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0) * (3.0d0*ux(i,j)*dxpsi(i,j) + 3.0d0*uy(i,j)*dypsi(i,j) + 4.0d0*dxux(i,j) + 4.0d0*dyuy(i,j) ) &
          + 2.0d0*ux(i,j)*25.0d0/28.0d0 * ( (ux(i,j)**2 +1.0d0)*dxpsi(i,j) + 4.0d0*ux(i,j)*dxux(i,j) + uy(i,j)*(ux(i,j)*dypsi(i,j) + 4.0d0*dyux(i,j)) ) & 
          + 2.0d0*uy(i,j)*25.0d0/28.0d0 * ( uy(i,j)*ux(i,j)*dxpsi(i,j) + (uy(i,j)**2)*dypsi(i,j) + 4.0d0*uy(i,j)*dyuy(i,j) + 4.0d0*ux(i,j)*dxuy(i,j) + dypsi(i,j) ))

      b2(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/4.0d0/3.0d0*ux(i,j)*ut(i,j) * (3.0d0*ux(i,j)*dxpsi(i,j) + 3.0d0*uy(i,j)*dypsi(i,j) + 4.0d0*dyuy(i,j) + 4.0d0*dxux(i,j)) &
          + ((ut(i,j) + (ux(i,j)**2)/ut(i,j))*25.0d0/28.0d0) * (((ux(i,j)**2)+1.0d0)*dxpsi(i,j) + 4.0d0*ux(i,j)*dxux(i,j) + uy(i,j)*(ux(i,j)*dypsi(i,j) + 4.0d0*dyux(i,j))) &
          + 25.0d0/28.0d0*ux(i,j)*uy(i,j)/ut(i,j) * (uy(i,j)*ux(i,j)*dxpsi(i,j) + (uy(i,j)**2)*dypsi(i,j) + 4.0d0*uy(i,j)*dyuy(i,j) + 4.0d0*ux(i,j)*dxuy(i,j) +dypsi(i,j)) &
          + 2.0d0/3.0d0*ux(i,j)/ut(i,j)**2 * (ut(i,j)*(-2.0d0*(ux(i,j)**2 + 1.0d0)*dxux(i,j) + (ux(i,j)**2)*dyuy(i,j) - 3.0d0*ux(i,j)*uy(i,j)*dyux(i,j) + dyuy(i,j))) & 
          - uy(i,j)/(3.0d0*ut(i,j)**2) * (ut(i,j)*(3.0d0*(ux(i,j)**2)*dxuy(i,j) + ux(i,j)*uy(i,j)*(dxux(i,j)+dyuy(i,j)) + 3.0d0*dyux(i,j)*(uy(i,j)**2 + 1.0d0) + 3.0d0*dxuy(i,j))) )

      b3(i,j) = exp(3.0d0/4.0d0*psi(i,j))*(25.0d0/4.0d0/3.0d0*uy(i,j)*ut(i,j) * (3.0d0*ux(i,j)*dxpsi(i,j) + 3.0d0*uy(i,j)*dypsi(i,j) + 4.0d0*dyuy(i,j) + 4.0d0*dxux(i,j) ) &
          + ((ut(i,j) + (uy(i,j)**2)/ut(i,j))*25.0d0/28.0d0) * (uy(i,j)*ux(i,j)*dxpsi(i,j) + (uy(i,j)**2)*dypsi(i,j) + 4.0d0*uy(i,j)*dyuy(i,j) + 4.0d0*ux(i,j)*dxuy(i,j) +dypsi(i,j)) &
          + 25.0d0/28.0d0*uy(i,j)*ux(i,j)/ut(i,j) * (((ux(i,j)**2)+1.0d0)*dxpsi(i,j) + 4.0d0*ux(i,j)*dxux(i,j) + uy(i,j)*(ux(i,j)*dypsi(i,j) + 4.0d0*dyux(i,j))) &
          + 2.0d0/3.0d0*uy(i,j)/ut(i,j)**2 * (-3.0d0*ut(i,j)*ux(i,j)*uy(i,j)*dxuy(i,j) + ut(i,j)*dxux(i,j)*(uy(i,j)**2 + 1.0d0) - 2.0d0*ut(i,j)*(uy(i,j)**2 + 1.0d0)*dyuy(i,j)) & 
          - ux(i,j)/(3.0d0*ut(i,j)**2) * (ut(i,j)*(3.0d0*(ux(i,j)**2)*dxuy(i,j) + ux(i,j)*uy(i,j)*(dxux(i,j)+dyuy(i,j)) + 3.0d0*dyux(i,j)*(uy(i,j)**2 + 1.0d0) + 3.0d0*dxuy(i,j)))) 
     enddo
   enddo



   ! Inversión de la matriz A
   do j = 4, ny-3
    do i = 4, nx-3
! Calculate the determinant for each (i,j) matrix
       det = a11(i,j) * (a22(i,j) * a33(i,j) - a23(i,j) * a32(i,j)) - a12(i,j) * (a21(i,j) * a33(i,j) - a23(i,j) * a31(i,j)) + a13(i,j) * (a21(i,j) * a32(i,j) - a22(i,j) * a31(i,j))

       
! Check for singularity
       if (abs(det) < 1.0e-35) then 
         print *, char(7)
         print *, 'Matrix is singular or near singular at index (', i, ',', j, ')', det
         stop 

        
       end if



! Compute the inverse matrix components for (i,j)
       inv11(i,j)=(a22(i,j)*a33(i,j)-a23(i,j)*a32(i,j))/det
       inv12(i,j)=-(a21(i,j)*a33(i,j)-a23(i,j)*a31(i,j))/det
       inv13(i,j)=(a21(i,j)*a32(i,j)-a22(i,j)*a31(i,j))/det
       inv21(i,j)=-(a12(i,j)*a33(i,j)-a13(i,j)*a32(i,j))/det
       inv22(i,j)=(a11(i,j)*a33(i,j)-a13(i,j)*a31(i,j))/det
       inv23(i,j)=-(a11(i,j)*a32(i,j)-a12(i,j)*a31(i,j))/det
       inv31(i,j)=(a12(i,j)*a23(i,j)-a13(i,j)*a22(i,j))/det
       inv32(i,j)=-(a11(i,j)*a23(i,j)-a13(i,j)*a21(i,j))/det
       inv33(i,j)=(a11(i,j)*a22(i,j)-a12(i,j)*a21(i,j))/det
    end do
 end do


   t=1 !t=0 ecuacion 33
        !t=1 ecuacion 38-39 (algoritmo adaptativo)

   if (t==0) then
    do j = 4, ny-3
      do i = 4, nx-3


      q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
      q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
      q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))

      dtpsi(i,j) = (inv11(i,j)*(q11(i,j)-q01(i,j))/eta + inv12(i,j)*(q12(i,j)-q02(i,j))/eta + inv13(i,j)*(q13(i,j)-q03(i,j))/eta) - inv11(i,j)*b1(i,j) - inv12(i,j)*b2(i,j) - inv13(i,j)*b3(i,j)
      dtux(i,j) = (inv21(i,j)*(q11(i,j)-q01(i,j))/eta + inv22(i,j)*(q12(i,j)-q02(i,j))/eta + inv23(i,j)*(q13(i,j)-q03(i,j))/eta) - inv21(i,j)*b1(i,j) - inv22(i,j)*b2(i,j) - inv23(i,j)*b3(i,j)
      dtuy(i,j) = (inv31(i,j)*(q11(i,j)-q01(i,j))/eta + inv32(i,j)*(q12(i,j)-q02(i,j))/eta + inv33(i,j)*(q13(i,j)-q03(i,j))/eta) - inv31(i,j)*b1(i,j) - inv32(i,j)*b2(i,j) - inv33(i,j)*b3(i,j)
      
      !psi(i,j) = psi(i,j) + dt*dtpsi(i,j)
      !ux(i,j) = ux(i,j) + dt*dtux(i,j)
      !uy(i,j) = uy(i,j) + dt*dtuy(i,j)


      !ut(i,j) = sqrt(1.0d0 + ux(i,j)**2 + uy(i,j)**2)
      dens(i,j) = jt(i,j)/ut(i,j) !en cada timestep recuperamos densidad a partir de las variables primitivas ux y uy y de la variable conservada Jt (habra que evolucionarla)









    enddo
    enddo



   elseif(t==1) then

    do j = 4, ny-3
      do i = 4, nx-3

      q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2 - 1.0d0/3.0d0)
      q02(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*ux(i,j))
      q03(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)*uy(i,j))


      dtpsi(i,j) = (inv11(i,j)*(q11(i,j)-q01(i,j))/eta + inv12(i,j)*(q12(i,j)-q02(i,j))/eta + inv13(i,j)*(q13(i,j)-q03(i,j))/eta) - inv11(i,j)*b1(i,j) - inv12(i,j)*b2(i,j) - inv13(i,j)*b3(i,j)
      dtux(i,j) = (inv21(i,j)*(q11(i,j)-q01(i,j))/eta + inv22(i,j)*(q12(i,j)-q02(i,j))/eta + inv23(i,j)*(q13(i,j)-q03(i,j))/eta) - inv21(i,j)*b1(i,j) - inv22(i,j)*b2(i,j) - inv23(i,j)*b3(i,j)
      dtuy(i,j) = (inv31(i,j)*(q11(i,j)-q01(i,j))/eta + inv32(i,j)*(q12(i,j)-q02(i,j))/eta + inv33(i,j)*(q13(i,j)-q03(i,j))/eta) - inv31(i,j)*b1(i,j) - inv32(i,j)*b2(i,j) - inv33(i,j)*b3(i,j)
      
   
      
      psipf(i,j) = -(4.0d0/3.0d0*ut(i,j)*(ux(i,j)*exp(psi(i,j))*dxpsi(i,j) + uy(i,j)*exp(psi(i,j))*dypsi(i,j) + exp(psi(i,j))*(dxux(i,j)+dyuy(i,j))))/(exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2-1.0d0/3.0d0))
      uxpf(i,j) = -(ut(i,j)*ux(i,j)*exp(psi(i,j))*psipf(i,j) + exp(psi(i,j))*dxpsi(i,j)*ux(i,j)**2 + exp(psi(i,j))*2.0d0*ux(i,j)*dxux(i,j) + exp(psi(i,j))*dxpsi(i,j)/12.0d0 + exp(psi(i,j))*dypsi(i,j)*ux(i,j)*uy(i,j) + dyux(i,j)*exp(psi(i,j))*uy(i,j) + exp(psi(i,j))*ux(i,j)*dyuy(i,j))/(ut(i,j)*exp(psi(i,j)))
      uypf(i,j) = -(ut(i,j)*uy(i,j)*exp(psi(i,j))*psipf(i,j) + exp(psi(i,j))*dypsi(i,j)*uy(i,j)**2 + exp(psi(i,j))*2.0d0*uy(i,j)*dyuy(i,j) + exp(psi(i,j))*dypsi(i,j)/12.0d0 + exp(psi(i,j))*dxpsi(i,j)*ux(i,j)*uy(i,j) + dxux(i,j)*exp(psi(i,j))*uy(i,j) + exp(psi(i,j))*ux(i,j)*dxuy(i,j))/(ut(i,j)*exp(psi(i,j)))


      
      q11hat(i,j) = eta * (a11(i,j)*(dtpsi(i,j)-psipf(i,j)) + a12(i,j)*(dtux(i,j)-uxpf(i,j)) + a13(i,j)*(dtuy(i,j) - uypf(i,j)) )
      q12hat(i,j) = eta * (a21(i,j)*(dtpsi(i,j)-psipf(i,j)) + a22(i,j)*(dtux(i,j)-uxpf(i,j)) + a23(i,j)*(dtuy(i,j) - uypf(i,j)) )
      q13hat(i,j) = eta * (a31(i,j)*(dtpsi(i,j)-psipf(i,j)) + a32(i,j)*(dtux(i,j)-uxpf(i,j)) + a33(i,j)*(dtuy(i,j) - uypf(i,j)) )








      if (sqrt(q11hat(i,j)**2+q12hat(i,j)**2+q13hat(i,j)**2)>= deltaeta) then 

      dtpsi(i,j) = 1.0d0/eta*(inv11(i,j)*q11hat(i,j) + inv12(i,j)*q12hat(i,j) + inv13(i,j)*q13hat(i,j)) + psipf(i,j)
      dtux(i,j) = 1.0d0/eta*(inv21(i,j)*q11hat(i,j) + inv22(i,j)*q12hat(i,j) + inv23(i,j)*q13hat(i,j)) + uxpf(i,j)
      dtuy(i,j) = 1.0d0/eta*(inv31(i,j)*q11hat(i,j) + inv32(i,j)*q12hat(i,j) + inv33(i,j)*q13hat(i,j)) + uypf(i,j)


      !psi(i,j) = psi(i,j) + dt*dtpsi(i,j)
      !ux(i,j) = ux(i,j) + dt*dtux(i,j)
      !uy(i,j) = uy(i,j) + dt*dtuy(i,j)


      
      !ut(i,j) = sqrt(1.0d0 + ux(i,j)**2 + uy(i,j)**2)
      dens(i,j) = jt(i,j)/ut(i,j)


      elseif (sqrt(q11hat(i,j)**2+q12hat(i,j)**2+q13hat(i,j)**2)< deltaeta) then

        dtpsi(i,j) = 0.0!psipf(i,j)
        dtux(i,j) = 0.0!uxpf(i,j)
        dtuy(i,j) = 0.0!uypf(i,j)
        



  
      epspf(i,j) = -q11(i,j) + sqrt(6.0d0*q11(i,j)**2 + 3.0d0*(q11(i,j)**2-q12(i,j)**2-q13(i,j)**2))
      psipf(i,j) = log(epspf(i,j))
      absvpf(i,j) = (sqrt(q12(i,j)**2 + q13(i,j)**2))/(q11(i,j) + 3.0d0*epspf(i,j))
      utpf(i,j) = 1.0d0/sqrt(1.0d0 - absvpf(i,j)**2)
      uxpf(i,j) = (3.0d0*utpf(i,j)*q12(i,j))/(3.0d0*q11(i,j) + epspf(i,j))
      uypf(i,j) = (3.0d0*utpf(i,j)*q13(i,j))/(3.0d0*q11(i,j) + epspf(i,j))

      psi(i,j) = psipf(i,j)
      ux(i,j) = uxpf(i,j)
      uy(i,j) = uypf(i,j)
      ut(i,j) = utpf(i,j)

      ut(i,j) = sqrt(1.0d0 + ux(i,j)**2 + uy(i,j)**2)
      dens(i,j) = jt(i,j)/ut(i,j)

      endif
    enddo
   enddo
   endif








endif


end subroutine get_primitives
 