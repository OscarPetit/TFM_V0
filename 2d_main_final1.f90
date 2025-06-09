program main
    implicit none
  
    ! Grid and time parameters
    integer, parameter :: nx = 2**9 +1  , ny = 2*nx
    real(kind=8), parameter :: x_min = -0.5d0, x_max = -x_min
    real(kind=8), parameter :: y_min = -1.5d0, y_max = -y_min
    integer, parameter :: it_min = 1, it_max = 600000
  
    real(kind=8) :: dx, dy, dt, time, total_mass, c, deltaeta
    integer :: i, j, it, test

    real(kind=8) :: x(nx), y(ny)

 !-------------------------------------------------------------------
  ! Variables físicas y derivadas
  !-------------------------------------------------------------------
    real(kind=8) :: psi(nx,ny), dens(nx,ny), ux(nx,ny), uy(nx,ny), ut(nx,ny), eta(nx,ny)
    real(kind=8) :: dtpsi(nx,ny), dtux(nx,ny), dtuy(nx,ny)
    real(kind=8) :: dxpsi(nx,ny), dxux(nx,ny), dxuy(nx,ny)
    real(kind=8) :: dypsi(nx,ny), dyux(nx,ny), dyuy(nx,ny)
    real(kind=8) :: psiaux(nx,ny), uxaux(nx,ny), uyaux(nx,ny)
    real(kind=8) :: dtpsiaux(nx,ny), dtuxaux(nx,ny), dtuyaux(nx,ny)
 
  
    real(kind=8) :: q11(nx,ny), q12(nx,ny), q13(nx,ny)
    ! real(kind=8) :: q01(nx,ny), q02(nx,ny), q03(nx,ny)
    real(kind=8) :: jt(nx,ny)

    real(kind=8) :: psi_lx(nx,ny), psi_rx(nx,ny), psi_ly(nx,ny), psi_ry(nx,ny)
    real(kind=8) :: ux_lx(nx,ny), ux_rx(nx,ny), ux_ly(nx,ny), ux_ry(nx,ny)
    real(kind=8) :: uy_lx(nx,ny), uy_rx(nx,ny), uy_ly(nx,ny), uy_ry(nx,ny)
    real(kind=8), parameter :: alpha_visc = 0.0d0 !##################################################numerical visc parameter################
    real(kind=8) :: lap_jt, lap_q11, lap_q12, lap_q13
  
  
    !-------------------------------------------------------------------
    ! Variables de reconstrucción
    !-------------------------------------------------------------------
    real(kind=8) :: dens_lx(nx,ny), dens_rx(nx,ny), dens_ly(nx,ny), dens_ry(nx,ny)
    real(kind=8) :: ut_lx(nx,ny), ut_rx(nx,ny), ut_ly(nx,ny), ut_ry(nx,ny)
    real(kind=8) :: dtpsi_lx(nx,ny), dtpsi_rx(nx,ny), dtpsi_ly(nx,ny), dtpsi_ry(nx,ny)
    real(kind=8) :: dtux_lx(nx,ny), dtux_rx(nx,ny), dtux_ly(nx,ny), dtux_ry(nx,ny)
    real(kind=8) :: dtuy_lx(nx,ny), dtuy_rx(nx,ny), dtuy_ly(nx,ny), dtuy_ry(nx,ny)
    real(kind=8) :: dxpsi_lx(nx,ny), dxpsi_rx(nx,ny), dxpsi_ly(nx,ny), dxpsi_ry(nx,ny)
    real(kind=8) :: dxux_lx(nx,ny), dxux_rx(nx,ny), dxux_ly(nx,ny), dxux_ry(nx,ny)
    real(kind=8) :: dxuy_lx(nx,ny), dxuy_rx(nx,ny), dxuy_ly(nx,ny), dxuy_ry(nx,ny)
    real(kind=8) :: dypsi_lx(nx,ny), dypsi_rx(nx,ny), dypsi_ly(nx,ny), dypsi_ry(nx,ny)
    real(kind=8) :: dyux_lx(nx,ny), dyux_rx(nx,ny), dyux_ly(nx,ny), dyux_ry(nx,ny)
    real(kind=8) :: dyuy_lx(nx,ny), dyuy_rx(nx,ny), dyuy_ly(nx,ny), dyuy_ry(nx,ny)
  
    !-------------------------------------------------------------------
    ! Cantidades conservadas q reconstrucción
    !-------------------------------------------------------------------
    real(kind=8) :: q11_lxR(nx,ny), q11_rxR(nx,ny), q11_lyR(nx,ny), q11_ryR(nx,ny)
    real(kind=8) :: q12_lxR(nx,ny), q12_rxR(nx,ny), q12_lyR(nx,ny), q12_ryR(nx,ny)
    real(kind=8) :: q13_lxR(nx,ny), q13_rxR(nx,ny), q13_lyR(nx,ny), q13_ryR(nx,ny)
    ! real(kind=8) :: q01_lxR(nx,ny), q01_rxR(nx,ny), q01_lyR(nx,ny), q01_ryR(nx,ny)
    ! real(kind=8) :: q02_lxR(nx,ny), q02_rxR(nx,ny), q02_lyR(nx,ny), q02_ryR(nx,ny)
    ! real(kind=8) :: q03_lxR(nx,ny), q03_rxR(nx,ny), q03_lyR(nx,ny), q03_ryR(nx,ny)
  
    ! real(kind=8) :: q11_lxL(nx,ny), q11_rxL(nx,ny), q11_lyL(nx,ny), q11_ryL(nx,ny)
    ! real(kind=8) :: q12_lxL(nx,ny), q12_rxL(nx,ny), q12_lyL(nx,ny), q12_ryL(nx,ny)
    ! real(kind=8) :: q13_lxL(nx,ny), q13_rxL(nx,ny), q13_lyL(nx,ny), q13_ryL(nx,ny)
     real(kind=8) :: q01(nx,ny)
    ! real(kind=8) :: q02_lxL(nx,ny), q02_rxL(nx,ny), q02_lyL(nx,ny), q02_ryL(nx,ny)
    ! real(kind=8) :: q03_lxL(nx,ny), q03_rxL(nx,ny), q03_lyL(nx,ny), q03_ryL(nx,ny)
  
    !-------------------------------------------------------------------
    ! Variables de flujos
    !-------------------------------------------------------------------
    real(kind=8) :: A_lx(nx,ny), A_rx(nx,ny), A_ly(nx,ny), A_ry(nx,ny)
    real(kind=8) :: Qx_lx(nx,ny), Qx_rx(nx,ny), Qx_ly(nx,ny), Qx_ry(nx,ny)
    real(kind=8) :: Qy_lx(nx,ny), Qy_rx(nx,ny), Qy_ly(nx,ny), Qy_ry(nx,ny)
    real(kind=8) :: Qt_lx(nx,ny), Qt_rx(nx,ny), Qt_ly(nx,ny), Qt_ry(nx,ny)
    real(kind=8) :: sigmaxx_lx(nx,ny), sigmaxx_rx(nx,ny), sigmaxx_ly(nx,ny), sigmaxx_ry(nx,ny)
    real(kind=8) :: sigmaxy_lx(nx,ny), sigmaxy_rx(nx,ny), sigmaxy_ly(nx,ny), sigmaxy_ry(nx,ny)
    real(kind=8) :: sigmayy_lx(nx,ny), sigmayy_rx(nx,ny), sigmayy_ly(nx,ny), sigmayy_ry(nx,ny)
    real(kind=8) :: sigmatx_lx(nx,ny), sigmatx_rx(nx,ny), sigmatx_ly(nx,ny), sigmatx_ry(nx,ny)
    real(kind=8) :: sigmaty_lx(nx,ny), sigmaty_rx(nx,ny), sigmaty_ly(nx,ny), sigmaty_ry(nx,ny)
    real(kind=8) :: fxt_lxR(nx,ny), fxt_rxR(nx,ny), fxx_lxR(nx,ny), fxx_rxR(nx,ny), &
         fxy_lxR(nx,ny), fxy_rxR(nx,ny)
    real(kind=8) :: fdensx_lxR(nx,ny), fdensx_rxR(nx,ny)
    ! real(kind=8) :: fxt_lxL(nx,ny), fxt_rxL(nx,ny), fxx_lxL(nx,ny), fxx_rxL(nx,ny), &
    !      fxy_lxL(nx,ny), fxy_rxL(nx,ny)
    ! real(kind=8) :: f0xt_lxL(nx,ny), f0xt_rxL(nx,ny), f0xx_lxL(nx,ny), f0xx_rxL(nx,ny), &
    !      f0xy_lxL(nx,ny), f0xy_rxL(nx,ny)
    ! real(kind=8) :: fdensx_lxL(nx,ny), fdensx_rxL(nx,ny)
    real(kind=8) :: fyt_lyR(nx,ny), fyt_ryR(nx,ny), fyx_lyR(nx,ny), fyx_ryR(nx,ny), &
         fyy_lyR(nx,ny), fyy_ryR(nx,ny)
    real(kind=8) :: fdensy_lyR(nx,ny), fdensy_ryR(nx,ny)
    ! real(kind=8) :: fyt_lyL(nx,ny), fyt_ryL(nx,ny), fyx_lyL(nx,ny), fyx_ryL(nx,ny), &
    !      fyy_lyL(nx,ny), fyy_ryL(nx,ny)
    ! real(kind=8) :: f0yt_lyL(nx,ny), f0yt_ryL(nx,ny), f0yx_lyL(nx,ny), f0yx_ryL(nx,ny), &
    !      f0yy_lyL(nx,ny), f0yy_ryL(nx,ny)
    ! real(kind=8) :: fdensy_lyL(nx,ny), fdensy_ryL(nx,ny)
    real(kind=8) :: jt_lxR(nx,ny), jt_rxR(nx,ny), jt_lyR(nx,ny), jt_ryR(nx,ny)
    ! real(kind=8) :: jt_lxL(nx,ny), jt_rxL(nx,ny), jt_lyL(nx,ny), jt_ryL(nx,ny)
  
    !-------------------------------------------------------------------
    ! Funciones de flujo kt
    !-------------------------------------------------------------------
    real(kind=8) :: FxtR(nx,ny), FxxR(nx,ny), FxyR(nx,ny), FdensxR(nx,ny)
    real(kind=8) :: FytR(nx,ny), FyxR(nx,ny), FyyR(nx,ny), FdensyR(nx,ny)
    ! real(kind=8) :: FxtL(nx,ny), FxxL(nx,ny), FxyL(nx,ny), F0xtL(nx,ny), F0xxL(nx,ny), &
    !      F0xyL(nx,ny), FdensxL(nx,ny)
    ! real(kind=8) :: FytL(nx,ny), FyxL(nx,ny), FyyL(nx,ny), F0ytL(nx,ny), F0yxL(nx,ny), &
    !      F0yyL(nx,ny), FdensyL(nx,ny)
  
     real(kind=8) :: FxtauxR(nx,ny), FxxauxR(nx,ny), FxyauxR(nx,ny), FdensxauxR(nx,ny)
     real(kind=8) :: FytauxR(nx,ny), FyxauxR(nx,ny), FyyauxR(nx,ny), FdensyauxR(nx,ny)
    ! real(kind=8) :: FxtauxL(nx,ny), FxxauxL(nx,ny), FxyauxL(nx,ny), F0xtauxL(nx,ny), &
    !      F0xxauxL(nx,ny), F0xyauxL(nx,ny), FdensxauxL(nx,ny)
    ! real(kind=8) :: FytauxL(nx,ny), FyxauxL(nx,ny), FyyauxL(nx,ny), F0ytauxL(nx,ny), &
    !      F0yxauxL(nx,ny), F0yyauxL(nx,ny), FdensyauxL(nx,ny)
  
    !-------------------------------------------------------------------
    ! Funciones q auxiliares
    !-------------------------------------------------------------------
    real(kind=8) :: q11aux(nx,ny), &
         q12aux(nx,ny), q13aux(nx,ny), jtaux(nx,ny)!, q11aux2(nx,ny), &
         !q12aux2(nx,ny), q13aux2(nx,ny), jtaux2(nx,ny)
  


  
!---------------------------------------
!--------Generacion de la malla---------
!---------------------------------------


  ! Calculate grid spacing
  dx = (x_max - x_min) / (nx - 1)
  dy = (y_max - y_min) / (ny - 1)

  ! Generate grid points in x and y directions
  do i = 1, nx
    x(i) = x_min + (i - 1) * dx
  end do

  do j = 1, ny
    y(j) = y_min + (j - 1) * dy
  end do



!---------------------------------------
!--------Condiciones iniciales----------
!---------------------------------------

  test=1

  !test = 0 -> 2D Rotor (nx=ny=2.0d0**9, L=3)
  !test = 1 -> KH-Instability (ny=2.0d0*nx=2.0d0**9, Lx=0.5d0, Ly=1.5d0)
  !test = 2.0d0,3 -> Oblique-shockwave (nx=ny=2.0d0**9, L=200)
  !test = 4 -> gausiana
do i = 1, nx
  do j = 1, ny
    eta(i,j)=0.0d0
  end do
enddo
  deltaeta=0.001d0 !tolerancia
  c=1.0d0



  call initial_conditions(test, nx, ny, x, y, psi, dens, ux, uy, ut, dtpsi, dxpsi, dypsi, dtux, dxux, dyux, dtuy,dxuy, dyuy, q11, q12, q13, jt)

!------------------------------------------------------------
!----------Reconstrucion de variables interfases 1-----------
!------------------------------------------------------------

  call boundxy(psi, nx, ny)
  call boundxy(ux, nx, ny)
  call boundxy(uy, nx, ny)
  call boundxy(ut, nx, ny)
  call boundxy(dtpsi, nx, ny)
  call boundxy(dtux, nx, ny)
  call boundxy(dtuy, nx, ny)
  call boundxy(dens, nx, ny)

  call boundxy(jt, nx, ny)
  call boundxy(q11, nx, ny)
  call boundxy(q12, nx, ny)
  call boundxy(q13, nx, ny)




  !  call derivadas(nx, ny, dx, dy, psi, dxpsi, dypsi)
  !  call derivadas(nx, ny, dx, dy, ux, dxux, dyux)
  !  call derivadas(nx, ny, dx, dy, uy, dxuy, dyuy)




  
  ! call weno5_x(nx,ny,psi,psi_lx,psi_rx); call weno5_y(nx,ny,psi,psi_ly,psi_ry)
  ! call weno5_x(nx,ny,ux,ux_lx,ux_rx); call weno5_y(nx,ny,ux,ux_ly,ux_ry)
  ! call weno5_x(nx,ny,uy,uy_lx,uy_rx); call weno5_y(nx,ny,uy,uy_ly,uy_ry)
  ! call weno5_x(nx,ny,dens,dens_lx,dens_rx); call weno5_y(nx,ny,dens,dens_ly,dens_ry)
  
  ! call weno5_x(nx,ny,dtpsi,dtpsi_lx,dtpsi_rx); call weno5_y(nx,ny,dtpsi,dtpsi_ly,dtpsi_ry)
  ! call weno5_x(nx,ny,dtux,dtux_lx,dtux_rx); call weno5_y(nx,ny,dtux,dtux_ly,dtux_ry)
  ! call weno5_x(nx,ny,dtuy,dtuy_lx,dtuy_rx); call weno5_y(nx,ny,dtuy,dtuy_ly,dtuy_ry)

  ! call weno5_x(nx,ny,dxpsi,dxpsi_lx,dxpsi_rx); call weno5_y(nx,ny,dxpsi,dxpsi_ly,dxpsi_ry)
  ! call weno5_x(nx,ny,dxux,dxux_lx,dxux_rx); call weno5_y(nx,ny,dxux,dxux_ly,dxux_ry)
  ! call weno5_x(nx,ny,dxuy,dxuy_lx,dxuy_rx); call weno5_y(nx,ny,dxuy,dxuy_ly,dxuy_ry)

  ! call weno5_x(nx,ny,dypsi,dypsi_lx,dypsi_rx); call weno5_y(nx,ny,dypsi,dypsi_ly,dypsi_ry)
  ! call weno5_x(nx,ny,dyux,dyux_lx,dyux_rx); call weno5_y(nx,ny,dyux,dyux_ly,dyux_ry)
  ! call weno5_x(nx,ny,dyuy,dyuy_lx,dyuy_rx); call weno5_y(nx,ny,dyuy,dyuy_ly,dyuy_ry)

 call weno_2d_xL(nx, ny, psi, psi_lx); call weno_2d_xR(nx, ny, psi, psi_rx); call weno_2d_yL(nx, ny, psi, psi_ly); call weno_2d_yR(nx, ny, psi, psi_ry)
 call weno_2d_xL(nx, ny, ux, ux_lx); call weno_2d_xR(nx, ny, ux, ux_rx); call weno_2d_yL(nx, ny, ux, ux_ly); call weno_2d_yR(nx, ny, ux, ux_ry)
 call weno_2d_xL(nx, ny, uy, uy_lx); call weno_2d_xR(nx, ny, uy, uy_rx); call weno_2d_yL(nx, ny, uy, uy_ly); call weno_2d_yR(nx, ny, uy, uy_ry)
 call weno_2d_xL(nx, ny, dens, dens_lx); call weno_2d_xR(nx, ny, dens, dens_rx); call weno_2d_yL(nx, ny, dens, dens_ly); call weno_2d_yR(nx, ny, dens, dens_ry)

 call weno_2d_xL(nx, ny, dtpsi, dtpsi_lx); call weno_2d_xR(nx, ny, dtpsi, dtpsi_rx); call weno_2d_yL(nx, ny, dtpsi, dtpsi_ly); call weno_2d_yR(nx, ny, dtpsi, dtpsi_ry)
 call weno_2d_xL(nx, ny, dxpsi, dxpsi_lx); call weno_2d_xR(nx, ny, dxpsi, dxpsi_rx); call weno_2d_yL(nx, ny, dxpsi, dxpsi_ly); call weno_2d_yR(nx, ny, dxpsi, dxpsi_ry)
 call weno_2d_xL(nx, ny, dypsi, dypsi_lx); call weno_2d_xR(nx, ny, dypsi, dypsi_rx); call weno_2d_yL(nx, ny, dypsi, dypsi_ly); call weno_2d_yR(nx, ny, dypsi, dypsi_ry)

 call weno_2d_xL(nx, ny, dtux, dtux_lx); call weno_2d_xR(nx, ny, dtux, dtux_rx); call weno_2d_yL(nx, ny, dtux, dtux_ly); call weno_2d_yR(nx, ny, dtux, dtux_ry)
 call weno_2d_xL(nx, ny, dxux, dxux_lx); call weno_2d_xR(nx, ny, dxux, dxux_rx); call weno_2d_yL(nx, ny, dxux, dxux_ly); call weno_2d_yR(nx, ny, dxux, dxux_ry)
 call weno_2d_xL(nx, ny, dyux, dyux_lx); call weno_2d_xR(nx, ny, dyux, dyux_rx); call weno_2d_yL(nx, ny, dyux, dyux_ly); call weno_2d_yR(nx, ny, dyux, dyux_ry)
 
 
 call weno_2d_xL(nx, ny, dtuy, dtuy_lx); call weno_2d_xR(nx, ny, dtuy, dtuy_rx); call weno_2d_yL(nx, ny, dtuy, dtuy_ly); call weno_2d_yR(nx, ny, dtuy, dtuy_ry)
 call weno_2d_xL(nx, ny, dxuy, dxuy_lx); call weno_2d_xR(nx, ny, dxuy, dxuy_rx); call weno_2d_yL(nx, ny, dxuy, dxuy_ly); call weno_2d_yR(nx, ny, dxuy, dxuy_ry)
 call weno_2d_xL(nx, ny, dyuy, dyuy_lx); call weno_2d_xR(nx, ny, dyuy, dyuy_rx); call weno_2d_yL(nx, ny, dyuy, dyuy_ly); call weno_2d_yR(nx, ny, dyuy, dyuy_ry)



! do i=1, nx
!   do j=1, ny
!     ut(i,j) = sqrt(1.0d0 + ux(i,j)**2.0d0 + uy(i,j)**2.0d0)
!   end do
! enddo


  
!   ! call weno5_x(nx,ny,ut,ut_lx,ut_rx); call weno5_y(nx,ny,ut,ut_ly,ut_ry)
!   call weno_2d_xL(nx, ny, ut, ut_lx); call weno_2d_xR(nx, ny, ut, ut_rx); call weno_2d_yL(nx, ny, ut, ut_ly); call weno_2d_yR(nx, ny, ut, ut_ry)
!   !call weno_2d( nx, ny, ut, ut_lx, ut_rx, ut_ly, ut_ry)   !meditar si esto se hace WENO o se calcula de nuevo con el método de bajo

  
  
  
  
  
  
     do i = 3, nx-3 
      do j = 3, ny-3
        !ut(i,j)    = sqrt(1.0d0 + ux(i,j)**2.0d0 + uy(i,j)**2.0d0)
        ut_rx(i,j) = sqrt(1.0d0 + ux_rx(i,j)**2.0d0 + uy_rx(i,j)**2.0d0)
        ut_lx(i,j) = sqrt(1.0d0 + ux_lx(i,j)**2.0d0 + uy_lx(i,j)**2.0d0)
        ut_ly(i,j) = sqrt(1.0d0 + ux_ly(i,j)**2.0d0 + uy_ly(i,j)**2.0d0)
        ut_ry(i,j) = sqrt(1.0d0 + ux_ry(i,j)**2.0d0 + uy_ry(i,j)**2.0d0)
      end do
    end do







!---------------------------------------
!--------Abrir ficheros .dat -----------
!---------------------------------------

open(unit=10, file='ux_output.dat', status='replace', action='write')
open(unit=11, file='uy_output.dat', status="replace", action="write")
open(unit=12, file='psi_output.dat', status="replace", action="write")
open(unit=13, file='dens_output.dat', status="replace", action="write")
!open(unit=14, file='jt_output.dat', status="replace", action="write")
open(unit=15, file='ut_output.dat', status="replace", action="write")

!---------------------------------------
!-------Bucle de iteracion temporal-----
!---------------------------------------

  time = 0.0d0
  it = 0
  dt=0.1d0*dx

!-------------------------------------------------------------------------------------------------------
!-----------------------------------------Inicio de programa--------------------------------------------
!-------------------------------------------------------------------------------------------------------


  do it = it_min, it_max
    print *, "Iteration: ", it, "Time: ", time

! if (time<10)then
!   dt=0.001*dx
!   else
!     dt=0.1*dx
!   endif



    
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------
    





!------------------------------------------------------------
!----------------- Calculo de flujos 1 ----------------------
!------------------------------------------------------------

    call calculo_funciones_flujos(nx, ny, psi_lx, psi_rx, psi_ly, psi_ry, ux_lx, ux_rx, ux_ly, ux_ry, uy_lx, uy_rx, uy_ly, uy_ry, ut_lx, ut_rx, ut_ly, ut_ry, dtpsi_lx, dtpsi_rx, dtpsi_ly, dtpsi_ry, dtux_lx, dtux_rx, dtux_ly, dtux_ry, dtuy_lx, dtuy_rx, dtuy_ly, dtuy_ry,dxpsi_lx, dxpsi_rx, dxpsi_ly, dxpsi_ry, dxux_lx, dxux_rx, dxux_ly, dxux_ry, dxuy_lx, dxuy_rx, dxuy_ly, dxuy_ry, dypsi_lx, dypsi_rx, dypsi_ly, dypsi_ry, dyux_lx, dyux_rx, dyux_ly, dyux_ry, dyuy_lx, dyuy_rx, dyuy_ly, dyuy_ry, eta, A_lx, A_rx, A_ly, A_ry, Qx_lx, Qx_rx, Qx_ly, Qx_ry, Qy_lx, Qy_rx, Qy_ly, Qy_ry, Qt_lx, Qt_rx, Qt_ly, Qt_ry, sigmaxx_lx, sigmaxx_rx, sigmaxx_ly, sigmaxx_ry, sigmaxy_lx, sigmaxy_rx, sigmaxy_ly, sigmaxy_ry, sigmayy_lx, sigmayy_rx, sigmayy_ly, sigmayy_ry, sigmatx_lx, sigmatx_rx, sigmatx_ly, sigmatx_ry, sigmaty_lx, sigmaty_rx, sigmaty_ly, sigmaty_ry)
    

    call flujos_xR(eta,nx, ny, psi_lx, psi_rx, ux_lx, ux_rx, uy_lx, uy_rx, ut_lx, ut_rx, A_lx, A_rx, Qx_lx, Qx_rx, Qy_lx, Qy_rx, Qt_lx, Qt_rx, sigmatx_lx, sigmatx_rx, sigmaxx_lx, sigmaxx_rx, sigmaxy_lx, sigmaxy_rx, dens_lx, dens_rx, fxt_lxR, fxt_rxR, fxx_lxR, fxx_rxR, fxy_lxR, fxy_rxR, fdensx_lxR, fdensx_rxR)
    call flujos_yR(eta,nx, ny, psi_ly, psi_ry, ux_ly, ux_ry, uy_ly, uy_ry, ut_ly, ut_ry, A_ly, A_ry, Qx_ly, Qx_ry, Qy_ly, Qy_ry, Qt_ly, Qt_ry, sigmaty_ly, sigmaty_ry, sigmaxy_ly, sigmaxy_ry, sigmayy_ly, sigmayy_ry, dens_ly, dens_ry, fyt_lyR, fyt_ryR, fyx_lyR, fyx_ryR, fyy_lyR, fyy_ryR, fdensy_lyR, fdensy_ryR)



!------------------------------------------------------------
!------------ Calculo conservadas interfase 1 ---------------
!------------------------------------------------------------
    do i = 3, nx-3 
      do j = 3, ny-3

            ! PARA INTERFASE DERECHA R

            q11_rxR(i,j) = (exp(psi_rx(i,j))+A_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)*ut_rx(i,j) - 1.0d0/3.0d0) + 2.0d0*ut_rx(i,j)*Qt_rx(i,j)
            q12_rxR(i,j) = (exp(psi_rx(i,j))+A_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)*ux_rx(i,j)) + ux_rx(i,j)*Qt_rx(i,j) + ut_rx(i,j)*Qx_rx(i,j) + sigmatx_rx(i,j)
            q13_rxR(i,j) = (exp(psi_rx(i,j))+A_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)*uy_rx(i,j)) + uy_rx(i,j)*Qt_rx(i,j) + ut_rx(i,j)*Qy_rx(i,j) + sigmaty_rx(i,j)
        
            q11_lxR(i,j) = (exp(psi_lx(i,j))+A_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)*ut_lx(i,j) - 1.0d0/3.0d0) + 2.0d0*ut_lx(i,j)*Qt_lx(i,j)
            q12_lxR(i,j) = (exp(psi_lx(i,j))+A_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)*ux_lx(i,j)) + ux_lx(i,j)*Qt_lx(i,j) + ut_lx(i,j)*Qx_lx(i,j) + sigmatx_lx(i,j)
            q13_lxR(i,j) = (exp(psi_lx(i,j))+A_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)*uy_lx(i,j)) + uy_lx(i,j)*Qt_lx(i,j) + ut_lx(i,j)*Qy_lx(i,j) + sigmaty_lx(i,j)
        
        
            q11_ryR(i,j) = (exp(psi_ry(i,j))+A_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)*ut_ry(i,j) - 1.0d0/3.0d0) + 2.0d0*ut_ry(i,j)*Qt_ry(i,j)
            q12_ryR(i,j) = (exp(psi_ry(i,j))+A_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)*ux_ry(i,j)) + ux_ry(i,j)*Qt_ry(i,j) + ut_ry(i,j)*Qx_ry(i,j) + sigmatx_ry(i,j)
            q13_ryR(i,j) = (exp(psi_ry(i,j))+A_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)*uy_ry(i,j)) + uy_ry(i,j)*Qt_ry(i,j) + ut_ry(i,j)*Qy_ry(i,j) + sigmaty_ry(i,j)
        
            q11_lyR(i,j) = (exp(psi_ly(i,j))+A_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)*ut_ly(i,j) - 1.0d0/3.0d0) + 2.0d0*ut_ly(i,j)*Qt_ly(i,j)
            q12_lyR(i,j) = (exp(psi_ly(i,j))+A_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)*ux_ly(i,j)) + ux_ly(i,j)*Qt_ly(i,j) + ut_ly(i,j)*Qx_ly(i,j) + sigmatx_ly(i,j)
            q13_lyR(i,j) = (exp(psi_ly(i,j))+A_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)*uy_ly(i,j)) + uy_ly(i,j)*Qt_ly(i,j) + ut_ly(i,j)*Qy_ly(i,j) + sigmaty_ly(i,j)
        
        
            jt_rxR(i,j) = dens_rx(i,j)*ut_rx(i,j)
            jt_lxR(i,j) = dens_lx(i,j)*ut_lx(i,j)
        
            jt_ryR(i,j) = dens_ry(i,j)*ut_ry(i,j)
            jt_lyR(i,j) = dens_ly(i,j)*ut_ly(i,j)
        
        
            ! q01_rxR(i,j) = exp(psi_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)**2.0d0 - 1.0d0/3.0d0)
            ! q02_rxR(i,j) = exp(psi_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)*ux_rx(i,j))
            ! q03_rxR(i,j) = exp(psi_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)*uy_rx(i,j))
        
            ! q01_lxR(i,j) = exp(psi_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)**2.0d0 - 1.0d0/3.0d0)
            ! q02_lxR(i,j) = exp(psi_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)*ux_lx(i,j))
            ! q03_lxR(i,j) = exp(psi_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)*uy_lx(i,j))
        
            ! q01_ryR(i,j) = exp(psi_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)**2.0d0 - 1.0d0/3.0d0)
            ! q02_ryR(i,j) = exp(psi_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)*ux_ry(i,j))
            ! q03_ryR(i,j) = exp(psi_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)*uy_ry(i,j))
        
            ! q01_lyR(i,j) = exp(psi_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)**2.0d0 - 1.0d0/3.0d0)
            ! q02_lyR(i,j) = exp(psi_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)*ux_ly(i,j))
            ! q03_lyR(i,j) = exp(psi_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)*uy_ly(i,j))


        enddo
    enddo



!------------------------------------------------------------
!------------ Calculo funciones de flujo 1 ------------------
!------------------------------------------------------------
    do i = 3, nx-3 
      do j = 3, ny-3

      ! INTERFASE DERECHA R

      FxtR(i,j) = 0.5d0*(fxt_rxR(i,j) + fxt_lxR(i,j) - c*(q11_rxR(i,j) - q11_lxR(i,j)))
      FxxR(i,j) = 0.5d0*(fxx_rxR(i,j) + fxx_lxR(i,j) - c*(q12_rxR(i,j) - q12_lxR(i,j)))
      FxyR(i,j) = 0.5d0*(fxy_rxR(i,j) + fxy_lxR(i,j) - c*(q13_rxR(i,j) - q13_lxR(i,j)))
  
      ! F0xtR(i,j) = 0.5d0*(f0xt_rxR(i,j) + f0xt_lxR(i,j) - c*(q01_rxR(i,j) - q01_lxR(i,j)))
      ! F0xxR(i,j) = 0.5d0*(f0xx_rxR(i,j) + f0xx_lxR(i,j) - c*(q02_rxR(i,j) - q02_lxR(i,j)))
      ! F0xyR(i,j) = 0.5d0*(f0xy_rxR(i,j) + f0xy_lxR(i,j) - c*(q03_rxR(i,j) - q03_lxR(i,j)))
  
      FdensxR(i,j) = 0.5d0*(fdensx_rxR(i,j) + fdensx_lxR(i,j) - c*(jt_rxR(i,j) - jt_lxR(i,j)))
  
      !En y
  
      FytR(i,j) = 0.5d0*(fyt_ryR(i,j) + fyt_lyR(i,j) - c*(q11_ryR(i,j) - q11_lyR(i,j)))
      FyxR(i,j) = 0.5d0*(fyx_ryR(i,j) + fyx_lyR(i,j) - c*(q12_ryR(i,j) - q12_lyR(i,j)))
      FyyR(i,j) = 0.5d0*(fyy_ryR(i,j) + fyy_lyR(i,j) - c*(q13_ryR(i,j) - q13_lyR(i,j)))
  
      ! F0ytR(i,j) = 0.5d0*(f0yt_ryR(i,j) + f0yt_lyR(i,j) - c*(q01_ryR(i,j) - q01_lyR(i,j)))
      ! F0yxR(i,j) = 0.5d0*(f0yx_ryR(i,j) + f0yx_lyR(i,j) - c*(q02_ryR(i,j) - q02_lyR(i,j)))
      ! F0yyR(i,j) = 0.5d0*(f0yy_ryR(i,j) + f0yy_lyR(i,j) - c*(q03_ryR(i,j) - q03_lyR(i,j)))
  
  
      FdensyR(i,j) = 0.5d0*(fdensy_ryR(i,j) + fdensy_lyR(i,j) - c*(jt_ryR(i,j) - jt_lyR(i,j)))



  
       enddo
    enddo

   
!-------------------------------------------------------------
!-----------Evolucion variables conservadas 1-----------------
!-------------------------------------------------------------
    do i = 4, nx-3 
      do j = 4, ny-3
  
        lap_jt = ( jt(i+1,j) - 2.0d0*jt(i,j) + jt(i-1,j) ) / dx**2  &
        + ( jt(i,j+1) - 2.0d0*jt(i,j) + jt(i,j-1) ) / dy**2

        lap_q11 = ( q11(i+1,j) - 2.0d0*q11(i,j) + q11(i-1,j) ) / dx**2  &
         + ( q11(i,j+1) - 2.0d0*q11(i,j) + q11(i,j-1) ) / dy**2

        lap_q12 = ( q12(i+1,j) - 2.0d0*q12(i,j) + q12(i-1,j) ) / dx**2  &
         + ( q12(i,j+1) - 2.0d0*q12(i,j) + q12(i,j-1) ) / dy**2

        lap_q13 = ( q13(i+1,j) - 2.0d0*q13(i,j) + q13(i-1,j) ) / dx**2  &
         + ( q13(i,j+1) - 2.0d0*q13(i,j) + q13(i,j-1) ) / dy**2
  
  
       jtaux(i,j) = jt(i,j)-dt/dx*(FdensxR(i,j)-FdensxR(i-1,j))-dt/dy*(FdensyR(i,j)-FdensyR(i,j-1)) + alpha_visc*dx*dy * lap_jt
  

       q11aux(i,j) = q11(i,j) - dt/dx*(FxtR(i,j)-FxtR(i-1,j)) -dt/dy*(FytR(i,j)-FytR(i,j-1))+ alpha_visc*dx*dy*  lap_q11
       q12aux(i,j) = q12(i,j) - dt/dx*(FxxR(i,j)-FxxR(i-1,j)) -dt/dy*(FyxR(i,j)-FyxR(i,j-1))+ alpha_visc*dx*dy*  lap_q12
       q13aux(i,j) = q13(i,j) - dt/dx*(FxyR(i,j)-FxyR(i-1,j)) -dt/dy*(FyyR(i,j)-FyyR(i,j-1))+ alpha_visc*dx*dy*  lap_q13
  
  

  
  

        enddo
    enddo

!-------------------------------------------------------------
!-----------Recuperacion de variables primitivas 1 ----------- !usamos auxiliares / correctores
!-------------------------------------------------------------

    do i = 4, nx-3 
      do j = 4, ny-3
      psiaux(i,j) = psi(i,j)
      uxaux(i,j) = ux(i,j)
      uyaux(i,j) = uy(i,j)

      dtpsiaux(i,j) = dtpsi(i,j)
      dtuxaux(i,j) = dtux(i,j)
      dtuyaux(i,j) = dtuy(i,j)

      
      psi(i,j) = psi(i,j) + dt*dtpsi(i,j)
      ux(i,j) = ux(i,j) + dt*dtux(i,j)
      uy(i,j) = uy(i,j) + dt*dtuy(i,j)


      ut(i,j) = sqrt(1.0d0 + ux(i,j)**2.0d0 + uy(i,j)**2.0d0)

      dens(i,j) = jtaux(i,j)/ut(i,j) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    


      enddo
    enddo
    

    !call get_primitives(eta, deltaeta, psi, dtpsi, dxpsi, dypsi, ux, dtux, dxux, dyux, uy, dtuy, dxuy, dyuy, ut, q11, q12, q13, dens, jt, nx, ny )



    

    !-----------------------------------------------------------------------REPETICIÓN 1------------------------------------------------------

!----------------------------------------
!--------Condiciones de contorno 2.0d0-------
!----------------------------------------

    call boundxy(psi, nx, ny)
    call boundxy(ux, nx, ny)
    call boundxy(uy, nx, ny)
    call boundxy(ut, nx, ny)
    call boundxy(dtpsi, nx, ny)
    call boundxy(dtux, nx, ny)
    call boundxy(dtuy, nx, ny)
    call boundxy(dens, nx, ny)
    call boundxy(jtaux, nx, ny)
    call boundxy(q11aux, nx, ny)
    call boundxy(q12aux, nx, ny)
    call boundxy(q13aux, nx, ny)




    call derivadas_x(nx, ny, dx, psi, eta, dxpsi); call derivadas_y(nx, ny, dy, psi, eta, dypsi)
    call derivadas_x(nx, ny, dx, ux, eta, dxux); call derivadas_y(nx, ny, dy, ux, eta, dyux)
    call derivadas_x(nx, ny, dx, uy, eta, dxuy); call derivadas_y(nx, ny, dy, uy, eta, dyuy)


     call compute_dtpsi(nx, ny, psi, ux, uy, dxpsi, dypsi, dxux, dyux, dxuy, dyuy, dtpsi,eta,q11aux,q12aux,q13aux,deltaeta)
     call compute_dtux(nx, ny, psi, ux, uy, dxpsi, dypsi, dxux, dyux, dxuy, dyuy, dtux,eta,q11aux,q12aux,q13aux,deltaeta)
     call compute_dtuy(nx, ny, psi, ux, uy, dxpsi, dypsi, dxux, dyux, dxuy, dyuy, dtuy,eta,q11aux,q12aux,q13aux,deltaeta)
  
     call boundxy(psi, nx, ny)
     call boundxy(ux, nx, ny)
     call boundxy(uy, nx, ny)
     call boundxy(ut, nx, ny)
     call boundxy(dtpsi, nx, ny)
     call boundxy(dtux, nx, ny)
     call boundxy(dtuy, nx, ny)
     call boundxy(dens, nx, ny)
     call boundxy(jtaux, nx, ny)
     call boundxy(q11aux, nx, ny)
     call boundxy(q12aux, nx, ny)
     call boundxy(q13aux, nx, ny)
 
  
  
    !  call weno5_x(nx,ny,psi,psi_lx,psi_rx); call weno5_y(nx,ny,psi,psi_ly,psi_ry)
    !  call weno5_x(nx,ny,ux,ux_lx,ux_rx); call weno5_y(nx,ny,ux,ux_ly,ux_ry)
    !  call weno5_x(nx,ny,uy,uy_lx,uy_rx); call weno5_y(nx,ny,uy,uy_ly,uy_ry)
    !  call weno5_x(nx,ny,dens,dens_lx,dens_rx); call weno5_y(nx,ny,dens,dens_ly,dens_ry)
     
    !  call weno5_x(nx,ny,dtpsi,dtpsi_lx,dtpsi_rx); call weno5_y(nx,ny,dtpsi,dtpsi_ly,dtpsi_ry)
    !  call weno5_x(nx,ny,dtux,dtux_lx,dtux_rx); call weno5_y(nx,ny,dtux,dtux_ly,dtux_ry)
    !  call weno5_x(nx,ny,dtuy,dtuy_lx,dtuy_rx); call weno5_y(nx,ny,dtuy,dtuy_ly,dtuy_ry)

    !  call weno5_x(nx,ny,dxpsi,dxpsi_lx,dxpsi_rx); call weno5_y(nx,ny,dxpsi,dxpsi_ly,dxpsi_ry)
    !  call weno5_x(nx,ny,dxux,dxux_lx,dxux_rx); call weno5_y(nx,ny,dxux,dxux_ly,dxux_ry)
    !  call weno5_x(nx,ny,dxuy,dxuy_lx,dxuy_rx); call weno5_y(nx,ny,dxuy,dxuy_ly,dxuy_ry)

    !  call weno5_x(nx,ny,dypsi,dypsi_lx,dypsi_rx); call weno5_y(nx,ny,dypsi,dypsi_ly,dypsi_ry)
    !  call weno5_x(nx,ny,dyux,dyux_lx,dyux_rx); call weno5_y(nx,ny,dyux,dyux_ly,dyux_ry)
    !  call weno5_x(nx,ny,dyuy,dyuy_lx,dyuy_rx); call weno5_y(nx,ny,dyuy,dyuy_ly,dyuy_ry)

    call weno_2d_xL(nx, ny, psi, psi_lx); call weno_2d_xR(nx, ny, psi, psi_rx); call weno_2d_yL(nx, ny, psi, psi_ly); call weno_2d_yR(nx, ny, psi, psi_ry)
    call weno_2d_xL(nx, ny, ux, ux_lx); call weno_2d_xR(nx, ny, ux, ux_rx); call weno_2d_yL(nx, ny, ux, ux_ly); call weno_2d_yR(nx, ny, ux, ux_ry)
    call weno_2d_xL(nx, ny, uy, uy_lx); call weno_2d_xR(nx, ny, uy, uy_rx); call weno_2d_yL(nx, ny, uy, uy_ly); call weno_2d_yR(nx, ny, uy, uy_ry)
    call weno_2d_xL(nx, ny, dens, dens_lx); call weno_2d_xR(nx, ny, dens, dens_rx); call weno_2d_yL(nx, ny, dens, dens_ly); call weno_2d_yR(nx, ny, dens, dens_ry)
  
    call weno_2d_xL(nx, ny, dtpsi, dtpsi_lx); call weno_2d_xR(nx, ny, dtpsi, dtpsi_rx); call weno_2d_yL(nx, ny, dtpsi, dtpsi_ly); call weno_2d_yR(nx, ny, dtpsi, dtpsi_ry)
    call weno_2d_xL(nx, ny, dxpsi, dxpsi_lx); call weno_2d_xR(nx, ny, dxpsi, dxpsi_rx); call weno_2d_yL(nx, ny, dxpsi, dxpsi_ly); call weno_2d_yR(nx, ny, dxpsi, dxpsi_ry)
    call weno_2d_xL(nx, ny, dypsi, dypsi_lx); call weno_2d_xR(nx, ny, dypsi, dypsi_rx); call weno_2d_yL(nx, ny, dypsi, dypsi_ly); call weno_2d_yR(nx, ny, dypsi, dypsi_ry)
  
    call weno_2d_xL(nx, ny, dtux, dtux_lx); call weno_2d_xR(nx, ny, dtux, dtux_rx); call weno_2d_yL(nx, ny, dtux, dtux_ly); call weno_2d_yR(nx, ny, dtux, dtux_ry)
    call weno_2d_xL(nx, ny, dxux, dxux_lx); call weno_2d_xR(nx, ny, dxux, dxux_rx); call weno_2d_yL(nx, ny, dxux, dxux_ly); call weno_2d_yR(nx, ny, dxux, dxux_ry)
    call weno_2d_xL(nx, ny, dyux, dyux_lx); call weno_2d_xR(nx, ny, dyux, dyux_rx); call weno_2d_yL(nx, ny, dyux, dyux_ly); call weno_2d_yR(nx, ny, dyux, dyux_ry)
    
    
    call weno_2d_xL(nx, ny, dtuy, dtuy_lx); call weno_2d_xR(nx, ny, dtuy, dtuy_rx); call weno_2d_yL(nx, ny, dtuy, dtuy_ly); call weno_2d_yR(nx, ny, dtuy, dtuy_ry)
    call weno_2d_xL(nx, ny, dxuy, dxuy_lx); call weno_2d_xR(nx, ny, dxuy, dxuy_rx); call weno_2d_yL(nx, ny, dxuy, dxuy_ly); call weno_2d_yR(nx, ny, dxuy, dxuy_ry)
    call weno_2d_xL(nx, ny, dyuy, dyuy_lx); call weno_2d_xR(nx, ny, dyuy, dyuy_rx); call weno_2d_yL(nx, ny, dyuy, dyuy_ly); call weno_2d_yR(nx, ny, dyuy, dyuy_ry)
  


    !  do i=1, nx
    !   do j=1, ny
    !     ut(i,j) = sqrt(1.0d0 + ux(i,j)**2.0d0 + uy(i,j)**2.0d0)
    !   end do
    ! enddo
    
    
      
    
    ! call weno5_x(nx,ny,ut,ut_lx,ut_rx); call weno5_y(nx,ny,ut,ut_ly,ut_ry)
   ! call weno_2d_xL(nx, ny, ut, ut_lx); call weno_2d_xR(nx, ny, ut, ut_rx); call weno_2d_yL(nx, ny, ut, ut_ly); call weno_2d_yR(nx, ny, ut, ut_ry)
    !call weno_2d( nx, ny, ut, ut_lx, ut_rx, ut_ly, ut_ry)   !meditar si esto se hace WENO o se calcula de nuevo con el método de bajo
  
    
    
    
  
  
  
  
  
    do i = 3, nx-3 
      do j = 3, ny-3
        !ut(i,j)    = sqrt(1.0d0 + ux(i,j)**2.0d0 + uy(i,j)**2.0d0)
        ut_rx(i,j) = sqrt(1.0d0 + ux_rx(i,j)**2.0d0 + uy_rx(i,j)**2.0d0)
        ut_lx(i,j) = sqrt(1.0d0 + ux_lx(i,j)**2.0d0 + uy_lx(i,j)**2.0d0)
        ut_ly(i,j) = sqrt(1.0d0 + ux_ly(i,j)**2.0d0 + uy_ly(i,j)**2.0d0)
        ut_ry(i,j) = sqrt(1.0d0 + ux_ry(i,j)**2.0d0 + uy_ry(i,j)**2.0d0)
      end do
    end do


!------------------------------------------------------------
!----------------- Calculo de flujos 2.0d0 ----------------------
!------------------------------------------------------------

    call calculo_funciones_flujos(nx, ny, psi_lx, psi_rx, psi_ly, psi_ry, ux_lx, ux_rx, ux_ly, ux_ry, uy_lx, uy_rx, uy_ly, uy_ry, ut_lx, ut_rx, ut_ly, ut_ry, dtpsi_lx, dtpsi_rx, dtpsi_ly, dtpsi_ry, dtux_lx, dtux_rx, dtux_ly, dtux_ry, dtuy_lx, dtuy_rx, dtuy_ly, dtuy_ry,dxpsi_lx, dxpsi_rx, dxpsi_ly, dxpsi_ry, dxux_lx, dxux_rx, dxux_ly, dxux_ry, dxuy_lx, dxuy_rx, dxuy_ly, dxuy_ry, dypsi_lx, dypsi_rx, dypsi_ly, dypsi_ry, dyux_lx, dyux_rx, dyux_ly, dyux_ry, dyuy_lx, dyuy_rx, dyuy_ly, dyuy_ry, eta, A_lx, A_rx, A_ly, A_ry, Qx_lx, Qx_rx, Qx_ly, Qx_ry, Qy_lx, Qy_rx, Qy_ly, Qy_ry, Qt_lx, Qt_rx, Qt_ly, Qt_ry, sigmaxx_lx, sigmaxx_rx, sigmaxx_ly, sigmaxx_ry, sigmaxy_lx, sigmaxy_rx, sigmaxy_ly, sigmaxy_ry, sigmayy_lx, sigmayy_rx, sigmayy_ly, sigmayy_ry, sigmatx_lx, sigmatx_rx, sigmatx_ly, sigmatx_ry, sigmaty_lx, sigmaty_rx, sigmaty_ly, sigmaty_ry)
  
    call flujos_xR(eta,nx, ny, psi_lx, psi_rx, ux_lx, ux_rx, uy_lx, uy_rx, ut_lx, ut_rx, A_lx, A_rx, Qx_lx, Qx_rx, Qy_lx, Qy_rx, Qt_lx, Qt_rx, sigmatx_lx, sigmatx_rx, sigmaxx_lx, sigmaxx_rx, sigmaxy_lx, sigmaxy_rx, dens_lx, dens_rx, fxt_lxR, fxt_rxR, fxx_lxR, fxx_rxR, fxy_lxR, fxy_rxR, fdensx_lxR, fdensx_rxR)
    call flujos_yR(eta,nx, ny, psi_ly, psi_ry, ux_ly, ux_ry, uy_ly, uy_ry, ut_ly, ut_ry, A_ly, A_ry, Qx_ly, Qx_ry, Qy_ly, Qy_ry, Qt_ly, Qt_ry, sigmaty_ly, sigmaty_ry, sigmaxy_ly, sigmaxy_ry, sigmayy_ly, sigmayy_ry, dens_ly, dens_ry, fyt_lyR, fyt_ryR, fyx_lyR, fyx_ryR, fyy_lyR, fyy_ryR, fdensy_lyR, fdensy_ryR)


!------------------------------------------------------------
!------------ Calculo primitivas interfase 2.0d0 ----------------
!------------------------------------------------------------
    do i = 3, nx-3 
      do j = 3, ny-3

            ! PARA INTERFASE DERECHA R

          
        q11_rxR(i,j) = (exp(psi_rx(i,j))+A_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)*ut_rx(i,j) - 1.0d0/3.0d0) + 2.0d0*ut_rx(i,j)*Qt_rx(i,j)
        q12_rxR(i,j) = (exp(psi_rx(i,j))+A_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)*ux_rx(i,j)) + ux_rx(i,j)*Qt_rx(i,j) + ut_rx(i,j)*Qx_rx(i,j) + sigmatx_rx(i,j)
        q13_rxR(i,j) = (exp(psi_rx(i,j))+A_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)*uy_rx(i,j)) + uy_rx(i,j)*Qt_rx(i,j) + ut_rx(i,j)*Qy_rx(i,j) + sigmaty_rx(i,j)
    
        q11_lxR(i,j) = (exp(psi_lx(i,j))+A_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)*ut_lx(i,j) - 1.0d0/3.0d0) + 2.0d0*ut_lx(i,j)*Qt_lx(i,j)
        q12_lxR(i,j) = (exp(psi_lx(i,j))+A_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)*ux_lx(i,j)) + ux_lx(i,j)*Qt_lx(i,j) + ut_lx(i,j)*Qx_lx(i,j) + sigmatx_lx(i,j)
        q13_lxR(i,j) = (exp(psi_lx(i,j))+A_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)*uy_lx(i,j)) + uy_lx(i,j)*Qt_lx(i,j) + ut_lx(i,j)*Qy_lx(i,j) + sigmaty_lx(i,j)
    
    
        q11_ryR(i,j) = (exp(psi_ry(i,j))+A_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)*ut_ry(i,j) - 1.0d0/3.0d0) + 2.0d0*ut_ry(i,j)*Qt_ry(i,j)
        q12_ryR(i,j) = (exp(psi_ry(i,j))+A_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)*ux_ry(i,j)) + ux_ry(i,j)*Qt_ry(i,j) + ut_ry(i,j)*Qx_ry(i,j) + sigmatx_ry(i,j)
        q13_ryR(i,j) = (exp(psi_ry(i,j))+A_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)*uy_ry(i,j)) + uy_ry(i,j)*Qt_ry(i,j) + ut_ry(i,j)*Qy_ry(i,j) + sigmaty_ry(i,j)
    
        q11_lyR(i,j) = (exp(psi_ly(i,j))+A_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)*ut_ly(i,j) - 1.0d0/3.0d0) + 2.0d0*ut_ly(i,j)*Qt_ly(i,j)
        q12_lyR(i,j) = (exp(psi_ly(i,j))+A_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)*ux_ly(i,j)) + ux_ly(i,j)*Qt_ly(i,j) + ut_ly(i,j)*Qx_ly(i,j) + sigmatx_ly(i,j)
        q13_lyR(i,j) = (exp(psi_ly(i,j))+A_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)*uy_ly(i,j)) + uy_ly(i,j)*Qt_ly(i,j) + ut_ly(i,j)*Qy_ly(i,j) + sigmaty_ly(i,j)
    
    
        jt_rxR(i,j) = dens_rx(i,j)*ut_rx(i,j)
        jt_lxR(i,j) = dens_lx(i,j)*ut_lx(i,j)
    
        jt_ryR(i,j) = dens_ry(i,j)*ut_ry(i,j)
        jt_lyR(i,j) = dens_ly(i,j)*ut_ly(i,j)
    
    
        ! q01_rxR(i,j) = exp(psi_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)**2.0d0 - 1.0d0/3.0d0)
        ! q02_rxR(i,j) = exp(psi_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)*ux_rx(i,j))
        ! q03_rxR(i,j) = exp(psi_rx(i,j))*(4.0d0/3.0d0*ut_rx(i,j)*uy_rx(i,j))
    
        ! q01_lxR(i,j) = exp(psi_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)**2.0d0 - 1.0d0/3.0d0)
        ! q02_lxR(i,j) = exp(psi_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)*ux_lx(i,j))
        ! q03_lxR(i,j) = exp(psi_lx(i,j))*(4.0d0/3.0d0*ut_lx(i,j)*uy_lx(i,j))
    
        ! q01_ryR(i,j) = exp(psi_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)**2.0d0 - 1.0d0/3.0d0)
        ! q02_ryR(i,j) = exp(psi_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)*ux_ry(i,j))
        ! q03_ryR(i,j) = exp(psi_ry(i,j))*(4.0d0/3.0d0*ut_ry(i,j)*uy_ry(i,j))
    
        ! q01_lyR(i,j) = exp(psi_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)**2.0d0 - 1.0d0/3.0d0)
        ! q02_lyR(i,j) = exp(psi_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)*ux_ly(i,j))
        ! q03_lyR(i,j) = exp(psi_ly(i,j))*(4.0d0/3.0d0*ut_ly(i,j)*uy_ly(i,j))


        enddo
    enddo



!------------------------------------------------------------
!------------ Calculo funciones de flujo 2.0d0 ------------------
!------------------------------------------------------------

    do i = 3, nx-3 
      do j = 3, ny-3

      ! INTERFASE DERECHA R

      FxtauxR(i,j) = 0.5d0*(fxt_rxR(i,j) + fxt_lxR(i,j) - c*(q11_rxR(i,j) - q11_lxR(i,j)))
      FxxauxR(i,j) = 0.5d0*(fxx_rxR(i,j) + fxx_lxR(i,j) - c*(q12_rxR(i,j) - q12_lxR(i,j)))
      FxyauxR(i,j) = 0.5d0*(fxy_rxR(i,j) + fxy_lxR(i,j) - c*(q13_rxR(i,j) - q13_lxR(i,j)))
  
      ! F0xtauxR(i,j) = 0.5d0*(f0xt_rxR(i,j) + f0xt_lxR(i,j) - c*(q01_rxR(i,j) - q01_lxR(i,j)))
      ! F0xxauxR(i,j) = 0.5d0*(f0xx_rxR(i,j) + f0xx_lxR(i,j) - c*(q02_rxR(i,j) - q02_lxR(i,j)))
      ! F0xyauxR(i,j) = 0.5d0*(f0xy_rxR(i,j) + f0xy_lxR(i,j) - c*(q03_rxR(i,j) - q03_lxR(i,j)))
  
      FdensxauxR(i,j) = 0.5d0*(fdensx_rxR(i,j) + fdensx_lxR(i,j) - c*(jt_rxR(i,j) - jt_lxR(i,j)))


  
      !En y
  
      FytauxR(i,j) = 0.5d0*(fyt_ryR(i,j) + fyt_lyR(i,j) - c*(q11_ryR(i,j) - q11_lyR(i,j)))
      FyxauxR(i,j) = 0.5d0*(fyx_ryR(i,j) + fyx_lyR(i,j) - c*(q12_ryR(i,j) - q12_lyR(i,j)))
      FyyauxR(i,j) = 0.5d0*(fyy_ryR(i,j) + fyy_lyR(i,j) - c*(q13_ryR(i,j) - q13_lyR(i,j)))
  
      ! F0ytauxR(i,j) = 0.5d0*(f0yt_ryR(i,j) + f0yt_lyR(i,j) - c*(q01_ryR(i,j) - q01_lyR(i,j)))
      ! F0yxauxR(i,j) = 0.5d0*(f0yx_ryR(i,j) + f0yx_lyR(i,j) - c*(q02_ryR(i,j) - q02_lyR(i,j)))
      ! F0yyauxR(i,j) = 0.5d0*(f0yy_ryR(i,j) + f0yy_lyR(i,j) - c*(q03_ryR(i,j) - q03_lyR(i,j)))
  
  
      FdensyauxR(i,j) = 0.5d0*(fdensy_ryR(i,j) + fdensy_lyR(i,j) - c*(jt_ryR(i,j) - jt_lyR(i,j)))




       enddo
    enddo





!print *, "AAAA", maxval(q11), maxval(q12), maxval(q13)
!------------------------------------------------------------
!------------ Calculo variables conservadas 2.0d0 ---------------
!------------------------------------------------------------
    do i = 4, nx-3 
      do j = 4, ny-3
  
      
        lap_jt = ( jtaux(i+1,j) - 2.0d0*jtaux(i,j) + jtaux(i-1,j) ) / dx**2  &
        + ( jtaux(i,j+1) - 2.0d0*jtaux(i,j) + jtaux(i,j-1) ) / dy**2

        lap_q11 = ( q11aux(i+1,j) - 2.0d0*q11aux(i,j) + q11aux(i-1,j) ) / dx**2  &
         + ( q11aux(i,j+1) - 2.0d0*q11aux(i,j) + q11aux(i,j-1) ) / dy**2

        lap_q12 = ( q12aux(i+1,j) - 2.0d0*q12aux(i,j) + q12aux(i-1,j) ) / dx**2  &
         + ( q12aux(i,j+1) - 2.0d0*q12aux(i,j) + q12aux(i,j-1) ) / dy**2

        lap_q13 = ( q13aux(i+1,j) - 2.0d0*q13aux(i,j) + q13aux(i-1,j) ) / dx**2  &
         + ( q13aux(i,j+1) - 2.0d0*q13aux(i,j) + q13aux(i,j-1) ) / dy**2

      jt(i,j) = jt(i,j)- dt/dx/2.0d0*(FdensxR(i,j)-FdensxR(i-1,j) + FdensxauxR(i,j)-FdensxauxR(i-1,j)) -dt/dy/2.0d0*(FdensyR(i,j)-FdensyR(i,j-1) + FdensyauxR(i,j)-FdensyauxR(i,j-1))+ alpha_visc*dx*dy* lap_jt
  

      q11(i,j) = q11(i,j) - dt/dx/2.0d0*(FxtR(i,j)-FxtR(i-1,j) + FxtauxR(i,j)-FxtauxR(i-1,j)) -dt/dy/2.0d0*(FytR(i,j)-FytR(i,j-1) + FytauxR(i,j)-FytauxR(i,j-1))+ alpha_visc*dx*dy*  lap_q11
      q12(i,j) = q12(i,j) - dt/dx/2.0d0*(FxxR(i,j)-FxxR(i-1,j) + FxxauxR(i,j)-FxxauxR(i-1,j)) -dt/dy/2.0d0*(FyxR(i,j)-FyxR(i,j-1) + FyxauxR(i,j)-FyxauxR(i,j-1))+ alpha_visc*dx*dy*  lap_q12
      q13(i,j) = q13(i,j) - dt/dx/2.0d0*(FxyR(i,j)-FxyR(i-1,j) + FxyauxR(i,j)-FxyauxR(i-1,j)) -dt/dy/2.0d0*(FyyR(i,j)-FyyR(i,j-1) + FyyauxR(i,j)-FyyauxR(i,j-1))+ alpha_visc*dx*dy*  lap_q13
  

        enddo
      enddo


!-------------------------------------------------------------
!-----------Recuperacion de variables primitivas 2.0d0 -----------
!-------------------------------------------------------------


      !   print *, "BBBB", maxval(q11), maxval(q12), maxval(q13)
      do i = 4, nx-3 
        do j = 4, ny-3
    psi(i,j) = psiaux(i,j) + dt/2.0d0*(dtpsi(i,j)+dtpsiaux(i,j))
    ux(i,j) = uxaux(i,j) + dt/2.0d0*(dtux(i,j)+dtuxaux(i,j))
    uy(i,j) = uyaux(i,j) + dt/2.0d0*(dtuy(i,j)+dtuyaux(i,j))
    

    ut(i,j) = sqrt(1.0d0 + ux(i,j)**2.0d0 + uy(i,j)**2.0d0)
    dens(i,j) = jt(i,j)/ut(i,j)

    !test

    enddo
  enddo





  
    !call get_primitives(eta, deltaeta, psi, dtpsi, dxpsi, dypsi, ux, dtux, dxux, dyux, uy, dtuy, dxuy, dyuy, ut, q11, q12, q13, dens, jt, nx, ny )



  
    
    call boundxy(psi, nx, ny)
    call boundxy(ux, nx, ny)
    call boundxy(uy, nx, ny)
    call boundxy(ut, nx, ny)
    call boundxy(dtpsi, nx, ny)
    call boundxy(dtux, nx, ny)
    call boundxy(dtuy, nx, ny)
    call boundxy(dens, nx, ny)
    call boundxy(jt, nx, ny)
    call boundxy(q11, nx, ny)
    call boundxy(q12, nx, ny)
    call boundxy(q13, nx, ny)




  
    call derivadas_x(nx, ny, dx, psi, eta, dxpsi); call derivadas_y(nx, ny, dy, psi, eta, dypsi)
    call derivadas_x(nx, ny, dx, ux, eta, dxux); call derivadas_y(nx, ny, dy, ux, eta, dyux)
    call derivadas_x(nx, ny, dx, uy, eta, dxuy); call derivadas_y(nx, ny, dy, uy, eta, dyuy)

 !aqui esta el error de la subida
     call compute_dtpsi(nx, ny, psi, ux, uy, dxpsi, dypsi, dxux, dyux, dxuy, dyuy, dtpsi,eta,q11,q12,q13,deltaeta)
     call compute_dtux(nx, ny, psi, ux, uy, dxpsi, dypsi, dxux, dyux, dxuy, dyuy, dtux,eta,q11,q12,q13,deltaeta)
     call compute_dtuy(nx, ny, psi, ux, uy, dxpsi, dypsi, dxux, dyux, dxuy, dyuy, dtuy,eta,q11,q12,q13,deltaeta)



     call boundxy(psi, nx, ny)
     call boundxy(ux, nx, ny)
     call boundxy(uy, nx, ny)
     call boundxy(ut, nx, ny)
     call boundxy(dtpsi, nx, ny)
     call boundxy(dtux, nx, ny)
     call boundxy(dtuy, nx, ny)
     call boundxy(dens, nx, ny)
     call boundxy(jtaux, nx, ny)
     call boundxy(jt, nx, ny)
     call boundxy(q11, nx, ny)
     call boundxy(q12, nx, ny)
     call boundxy(q13, nx, ny)
  

  


      ! call weno5_x(nx,ny,psi,psi_lx,psi_rx); call weno5_y(nx,ny,psi,psi_ly,psi_ry)
      ! call weno5_x(nx,ny,ux,ux_lx,ux_rx); call weno5_y(nx,ny,ux,ux_ly,ux_ry)
      ! call weno5_x(nx,ny,uy,uy_lx,uy_rx); call weno5_y(nx,ny,uy,uy_ly,uy_ry)
      ! call weno5_x(nx,ny,dens,dens_lx,dens_rx); call weno5_y(nx,ny,dens,dens_ly,dens_ry)
      
      ! call weno5_x(nx,ny,dtpsi,dtpsi_lx,dtpsi_rx); call weno5_y(nx,ny,dtpsi,dtpsi_ly,dtpsi_ry)
      ! call weno5_x(nx,ny,dtux,dtux_lx,dtux_rx); call weno5_y(nx,ny,dtux,dtux_ly,dtux_ry)
      ! call weno5_x(nx,ny,dtuy,dtuy_lx,dtuy_rx); call weno5_y(nx,ny,dtuy,dtuy_ly,dtuy_ry)

      ! call weno5_x(nx,ny,dxpsi,dxpsi_lx,dxpsi_rx); call weno5_y(nx,ny,dxpsi,dxpsi_ly,dxpsi_ry)
      ! call weno5_x(nx,ny,dxux,dxux_lx,dxux_rx); call weno5_y(nx,ny,dxux,dxux_ly,dxux_ry)
      ! call weno5_x(nx,ny,dxuy,dxuy_lx,dxuy_rx); call weno5_y(nx,ny,dxuy,dxuy_ly,dxuy_ry)

      ! call weno5_x(nx,ny,dypsi,dypsi_lx,dypsi_rx); call weno5_y(nx,ny,dypsi,dypsi_ly,dypsi_ry)
      ! call weno5_x(nx,ny,dyux,dyux_lx,dyux_rx); call weno5_y(nx,ny,dyux,dyux_ly,dyux_ry)
      ! call weno5_x(nx,ny,dyuy,dyuy_lx,dyuy_rx); call weno5_y(nx,ny,dyuy,dyuy_ly,dyuy_ry)

     call weno_2d_xL(nx, ny, psi, psi_lx); call weno_2d_xR(nx, ny, psi, psi_rx); call weno_2d_yL(nx, ny, psi, psi_ly); call weno_2d_yR(nx, ny, psi, psi_ry)
     call weno_2d_xL(nx, ny, ux, ux_lx); call weno_2d_xR(nx, ny, ux, ux_rx); call weno_2d_yL(nx, ny, ux, ux_ly); call weno_2d_yR(nx, ny, ux, ux_ry)
     call weno_2d_xL(nx, ny, uy, uy_lx); call weno_2d_xR(nx, ny, uy, uy_rx); call weno_2d_yL(nx, ny, uy, uy_ly); call weno_2d_yR(nx, ny, uy, uy_ry)
     call weno_2d_xL(nx, ny, dens, dens_lx); call weno_2d_xR(nx, ny, dens, dens_rx); call weno_2d_yL(nx, ny, dens, dens_ly); call weno_2d_yR(nx, ny, dens, dens_ry)
   
     call weno_2d_xL(nx, ny, dtpsi, dtpsi_lx); call weno_2d_xR(nx, ny, dtpsi, dtpsi_rx); call weno_2d_yL(nx, ny, dtpsi, dtpsi_ly); call weno_2d_yR(nx, ny, dtpsi, dtpsi_ry)
     call weno_2d_xL(nx, ny, dxpsi, dxpsi_lx); call weno_2d_xR(nx, ny, dxpsi, dxpsi_rx); call weno_2d_yL(nx, ny, dxpsi, dxpsi_ly); call weno_2d_yR(nx, ny, dxpsi, dxpsi_ry)
     call weno_2d_xL(nx, ny, dypsi, dypsi_lx); call weno_2d_xR(nx, ny, dypsi, dypsi_rx); call weno_2d_yL(nx, ny, dypsi, dypsi_ly); call weno_2d_yR(nx, ny, dypsi, dypsi_ry)
   
     call weno_2d_xL(nx, ny, dtux, dtux_lx); call weno_2d_xR(nx, ny, dtux, dtux_rx); call weno_2d_yL(nx, ny, dtux, dtux_ly); call weno_2d_yR(nx, ny, dtux, dtux_ry)
     call weno_2d_xL(nx, ny, dxux, dxux_lx); call weno_2d_xR(nx, ny, dxux, dxux_rx); call weno_2d_yL(nx, ny, dxux, dxux_ly); call weno_2d_yR(nx, ny, dxux, dxux_ry)
     call weno_2d_xL(nx, ny, dyux, dyux_lx); call weno_2d_xR(nx, ny, dyux, dyux_rx); call weno_2d_yL(nx, ny, dyux, dyux_ly); call weno_2d_yR(nx, ny, dyux, dyux_ry)
     
     
     call weno_2d_xL(nx, ny, dtuy, dtuy_lx); call weno_2d_xR(nx, ny, dtuy, dtuy_rx); call weno_2d_yL(nx, ny, dtuy, dtuy_ly); call weno_2d_yR(nx, ny, dtuy, dtuy_ry)
     call weno_2d_xL(nx, ny, dxuy, dxuy_lx); call weno_2d_xR(nx, ny, dxuy, dxuy_rx); call weno_2d_yL(nx, ny, dxuy, dxuy_ly); call weno_2d_yR(nx, ny, dxuy, dxuy_ry)
     call weno_2d_xL(nx, ny, dyuy, dyuy_lx); call weno_2d_xR(nx, ny, dyuy, dyuy_rx); call weno_2d_yL(nx, ny, dyuy, dyuy_ly); call weno_2d_yR(nx, ny, dyuy, dyuy_ry)
   

    !  do i=1, nx
    !   do j=1, ny
    !     ut(i,j) = sqrt(1.0d0 + ux(i,j)**2.0d0 + uy(i,j)**2.0d0)
    !   end do
    ! enddo
    
    
      
    ! call weno5_x(nx,ny,ut,ut_lx,ut_rx); call weno5_y(nx,ny,ut,ut_ly,ut_ry)
    !call weno_2d_xL(nx, ny, ut, ut_lx); call weno_2d_xR(nx, ny, ut, ut_rx); call weno_2d_yL(nx, ny, ut, ut_ly); call weno_2d_yR(nx, ny, ut, ut_ry)
    !call weno_2d( nx, ny, ut, ut_lx, ut_rx, ut_ly, ut_ry)   !meditar si esto se hace WENO o se calcula de nuevo con el método de bajo
  
    
    
    
  
  
  
  
  
  
    
    do i = 3, nx-3 
      do j = 3, ny-3
        !ut(i,j)    = sqrt(1.0d0 + ux(i,j)**2.0d0 + uy(i,j)**2.0d0)
        ut_rx(i,j) = sqrt(1.0d0 + ux_rx(i,j)**2.0d0 + uy_rx(i,j)**2.0d0)
        ut_lx(i,j) = sqrt(1.0d0 + ux_lx(i,j)**2.0d0 + uy_lx(i,j)**2.0d0)
        ut_ly(i,j) = sqrt(1.0d0 + ux_ly(i,j)**2.0d0 + uy_ly(i,j)**2.0d0)
        ut_ry(i,j) = sqrt(1.0d0 + ux_ry(i,j)**2.0d0 + uy_ry(i,j)**2.0d0)
      end do
    end do



    do i=4,nx-3
      do j=4,ny-3
        q01(i,j) = exp(psi(i,j))*(4.0d0/3.0d0*ut(i,j)**2.0d0 - 1.0d0/3.0d0)
      enddo
    enddo


!-------------------------------------------------------------
!-----------Final de bucle y escribir en archivos-------------
!-------------------------------------------------------------
    
! Initialize total_mass to zero at the start of the calculation
 !print *, "DEBUG", maxval(dtpsi), maxval(dtux), maxval(dtuy)  


    do i = 4, nx-3 
      do j = 4, ny-3
      total_mass = q11(i,j) * dx * dy
    end do
  end do

  call boundxy(q11, nx, ny)
  call boundxy(q01, nx, ny)



  print *, "MASS",total_mass, "DENS_max", MAXVAL(dens)
  print*, "qdiff", maxval(abs(q11 - q01)), minval(abs(q11 - q01))
  print*, "Fluxes", maxval(A_lx), maxval(Qt_lx), maxval(Qx_lx), maxval(Qy_lx), maxval(sigmatx_lx), maxval(sigmaty_lx), maxval(sigmaxx_lx), maxval(sigmaxy_lx), maxval(sigmayy_lx)
  print*, "Timespsi", maxval(dtpsiaux), maxval(dtpsi)
  print *, "----------------------------------------------------------------------------" 
  print *, "----------------------------------------------------------------------------" 
  
  
  time = time + dt


 ! PARAR EL PROGRAMA SI HAY NAN

  if (time>30) then
    print *, char(7)
    stop 
  endif

  if (any(dens /= dens) .or. (any(psi/=psi)))  then
    print *, "NaN detected at iteration ", it, ". Stopping simulation."
    print *, char(7)
    stop
  end if




  if ((mod(it,2500)==0) .or. (it == 1)) then



    ! Write data to files
    do i = 1, nx 
    do j = 1, ny
           write(10, *) time, x(i), y(j), ux(i,j)!ux(i, j)  ! Write ux values
           write(11, *) time, x(i), y(j), dypsi(i,j)!uy(i, j)  ! Write uy values
           write(12, *) time, x(i), y(j), psi(i, j)  ! Write psi values
           write(13, *) time, x(i), y(j), dens(i,j) !Write dens
           !write(14, *) time, i, j-4, jt(i,j)
           write(15, *) time, x(i), y(j), dxpsi(i,j) 
        end do
    end do
    
    
    endif
    
    
    
    end do
    
    ! Close files
    close(10)
    close(11)
    close(12)
    close(13)
    !close(14)
    close(15)

!--------------------------------------------------------------------------
! Deallocate all arrays at the end of the program
!--------------------------------------------------------------------------
    ! deallocate( psi, dens, ux, uy, ut )
    ! deallocate( dtpsi, dtux, dtuy, dxpsi, dxux, dxuy, dypsi, dyux, dyuy )
    ! deallocate( q11, q12, q13, q01, q02, q03, jt )
    
    ! ! Reconstruction arrays
    ! deallocate( dens_lx, dens_rx, dens_ly, dens_ry, &
    !             psi_lx,  psi_rx,  psi_ly,  psi_ry, &
    !             ux_lx,   ux_rx,   ux_ly,   ux_ry, &
    !             uy_lx,   uy_rx,   uy_ly,   uy_ry, &
    !             ut_lx,   ut_rx,   ut_ly,   ut_ry, &
    !             dtpsi_lx, dtpsi_rx, dtpsi_ly, dtpsi_ry, &
    !             dtux_lx,  dtux_rx,  dtux_ly,  dtux_ry, &
    !             dtuy_lx,  dtuy_rx,  dtuy_ly,  dtuy_ry, &
    !             dxpsi_lx, dxpsi_rx, dxpsi_ly, dxpsi_ry, &
    !             dxux_lx,  dxux_rx,  dxux_ly,  dxux_ry, &
    !             dxuy_lx,  dxuy_rx,  dxuy_ly,  dxuy_ry, &
    !             dypsi_lx, dypsi_rx, dypsi_ly, dypsi_ry, &
    !             dyux_lx,  dyux_rx,  dyux_ly,  dyux_ry, &
    !             dyuy_lx,  dyuy_rx,  dyuy_ly,  dyuy_ry )
    
    ! ! Q reconstruction arrays
    ! deallocate( q11_lxR, q11_rxR, q11_lyR, q11_ryR, &
    !             q12_lxR, q12_rxR, q12_lyR, q12_ryR, &
    !             q13_lxR, q13_rxR, q13_lyR, q13_ryR, &
    !             q01_lxR, q01_rxR, q01_lyR, q01_ryR, &
    !             q02_lxR, q02_rxR, q02_lyR, q02_ryR, &
    !             q03_lxR, q03_rxR, q03_lyR, q03_ryR, &
    !             q11_lxL, q11_rxL, q11_lyL, q11_ryL, &
    !             q12_lxL, q12_rxL, q12_lyL, q12_ryL, &
    !             q13_lxL, q13_rxL, q13_lyL, q13_ryL, &
    !             q01_lxL, q01_rxL, q01_lyL, q01_ryL, &
    !             q02_lxL, q02_rxL, q02_lyL, q02_ryL, &
    !             q03_lxL, q03_rxL, q03_lyL, q03_ryL )
    
    ! ! Flux arrays
    ! deallocate( A_lx, A_rx, A_ly, A_ry, &
    !             Qx_lx, Qx_rx, Qx_ly, Qx_ry, &
    !             Qy_lx, Qy_rx, Qy_ly, Qy_ry, &
    !             Qt_lx, Qt_rx, Qt_ly, Qt_ry, &
    !             sigmaxx_lx, sigmaxx_rx, sigmaxx_ly, sigmaxx_ry, &
    !             sigmaxy_lx, sigmaxy_rx, sigmaxy_ly, sigmaxy_ry, &
    !             sigmayy_lx, sigmayy_rx, sigmayy_ly, sigmayy_ry, &
    !             sigmatx_lx, sigmatx_rx, sigmatx_ly, sigmatx_ry, &
    !             sigmaty_lx, sigmaty_rx, sigmaty_ly, sigmaty_ry )
    
    ! ! More flux arrays
    ! deallocate( fxt_lxR, fxt_rxR, fxx_lxR, fxx_rxR, fxy_lxR, fxy_rxR, &
    !             f0xt_lxR, f0xt_rxR, f0xx_lxR, f0xx_rxR, f0xy_lxR, f0xy_rxR, &
    !             fdensx_lxR, fdensx_rxR, &
    !             fxt_lxL, fxt_rxL, fxx_lxL, fxx_rxL, fxy_lxL, fxy_rxL, &
    !             f0xt_lxL, f0xt_rxL, f0xx_lxL, f0xx_rxL, f0xy_lxL, f0xy_rxL, &
    !             fdensx_lxL, fdensx_rxL, &
    !             fyt_lyR, fyt_ryR, fyx_lyR, fyx_ryR, fyy_lyR, fyy_ryR, &
    !             f0yt_lyR, f0yt_ryR, f0yx_lyR, f0yx_ryR, f0yy_lyR, f0yy_ryR, &
    !             fdensy_lyR, fdensy_ryR, &
    !             fyt_lyL, fyt_ryL, fyx_lyL, fyx_ryL, fyy_lyL, fyy_ryL, &
    !             f0yt_lyL, f0yt_ryL, f0yx_lyL, f0yx_ryL, f0yy_lyL, f0yy_ryL, &
    !             fdensy_lyL, fdensy_ryL, &
    !             jt_lxR, jt_rxR, jt_lyR, jt_ryR, &
    !             jt_lxL, jt_rxL, jt_lyL, jt_ryL )
    
    ! ! kt flux arrays
    ! deallocate( FxtR, FxxR, FxyR, F0xtR, F0xxR, F0xyR, FdensxR, &
    !             FytR, FyxR, FyyR, F0ytR, F0yxR, F0yyR, FdensyR, &
    !             FxtL, FxxL, FxyL, F0xtL, F0xxL, F0xyL, FdensxL, &
    !             FytL, FyxL, FyyL, F0ytL, F0yxL, F0yyL, FdensyL, &
    !             FxtauxR, FxxauxR, FxyauxR, F0xtauxR, F0xxauxR, F0xyauxR, FdensxauxR, &
    !             FytauxR, FyxauxR, FyyauxR, F0ytauxR, F0yxauxR, F0yyauxR, FdensyauxR, &
    !             FxtauxL, FxxauxL, FxyauxL, F0xtauxL, F0xxauxL, F0xyauxL, FdensxauxL, &
    !             FytauxL, FyxauxL, FyyauxL, F0ytauxL, F0yxauxL, F0yyauxL, FdensyauxL )
    
    ! ! Auxiliary Q arrays
    ! deallocate( q01aux, q02aux, q03aux, q11aux, q12aux, q13aux, jtaux )
    


    


    



end program main