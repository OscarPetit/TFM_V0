subroutine boundxy(n, nx, ny)
   implicit none
   integer, intent(in) :: nx, ny
   real*8, intent(inout) :: n(nx, ny)

   integer :: i, j, t

   t=0 !0 copy
       !1 periodic
   
   if (t==0) then
 
   ! Left boundary ghost cells
   do j = 1, ny
     n(1, j) = n(4, j)
     n(2, j) = n(4, j)
     n(3, j) = n(4, j)

     n(nx, j) = n(nx-3, j)
     n(nx-1, j) = n(nx-3, j)
     n(nx-2, j) = n(nx-3, j)

   end do
 
   ! Bottom boundary ghost cells (1 = 3, 2 = 3)
   do i = 1, nx
     n(i, 3) = n(i, 4)
     n(i, 1) = n(i, 4)
     n(i, 2) = n(i, 4)

     n(i, ny) = n(i, ny-3)
     n(i, ny-1) = n(i, ny-3)
     n(i, ny-2) = n(i, ny-3)
    
   end do


  elseif (t==1) then
    
    do j = 1, ny
      n(1,j)    = n(nx-5,j)
      n(2,j)    = n(nx-4,j)
      n(3,j)    = n(nx-3,j)
     
      n(nx-2,j) = n(4,j)
      n(nx-1,j) = n(5,j)
      n(nx,  j) = n(6,j)
      
   end do
 
   do i = 1, nx
      n(i,1)    = n(i,ny-5)
      n(i,2)    = n(i,ny-4)
      n(i,3)    = n(i,ny-3)
     
      n(i,ny-2) = n(i,4)
      n(i,ny-1) = n(i,5)
      n(i,ny  ) = n(i,6)
      
   end do

endif
 
end subroutine boundxy
 