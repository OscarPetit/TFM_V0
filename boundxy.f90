subroutine boundxy(n, nx, ny)
   implicit none
   integer, intent(in) :: nx, ny
   real(kind=8), intent(inout) :: n(nx, ny)

   integer :: i, j, t

   t=1 !0 copy
       !1 periodic
        !2 periodic y absorbing x
        !3 absorbing
   
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


  elseif (t==2) then
    do j = 1, ny
      n(1,j)    = n(4,j)
      n(2,j)    = n(4,j)
      n(3,j)    = n(4,j)
     
      n(nx-2,j) = n(nx-3,j)
      n(nx-1,j) = n(nx-3,j)
      n(nx,  j) = n(nx-3,j)
      
   end do

   do i = 1, nx
    n(i,1)    = n(i,ny-5)
    n(i,2)    = n(i,ny-4)
    n(i,3)    = n(i,ny-3)
   
    n(i,ny-2) = n(i,4)
    n(i,ny-1) = n(i,5)
    n(i,ny  ) = n(i,6)
   end do


  elseif (t==3) then

    !— X‐direction, three ghost layers —!
    do j = 4, ny-3
  ! left ghosts i=3,2,1 from interior i=4,5
  n(3, j) = 2.0d0*n(4, j) - n(5, j)
  n(2, j) = 2.0d0*n(3, j) - n(4, j)
  n(1, j) = 2.0d0*n(2, j) - n(3, j)

  ! right ghosts i=nx-2,nx-1,nx from interior i=nx-3,nx-4
  n(nx-2, j) = 2.0d0*n(nx-3, j) - n(nx-4, j)
  n(nx-1, j) = 2.0d0*n(nx-2, j) - n(nx-3, j)
  n(nx  , j) = 2.0d0*n(nx-1, j) - n(nx-2, j)
end do

!— Y‐direction, three ghost layers —!
do i = 4, nx-3
  ! bottom ghosts j=3,2,1 from interior j=4,5
  n(i, 3 ) = 2.0d0*n(i, 4 ) - n(i, 5 )
  n(i, 2 ) = 2.0d0*n(i, 3 ) - n(i, 4 )
  n(i, 1 ) = 2.0d0*n(i, 2 ) - n(i, 3 )

  ! top ghosts j=ny-2,ny-1,ny from interior j=ny-3,ny-4
  n(i, ny-2) = 2.0d0*n(i, ny-3) - n(i, ny-4)
  n(i, ny-1) = 2.0d0*n(i, ny-2) - n(i, ny-3)
  n(i, ny  ) = 2.0d0*n(i, ny-1) - n(i, ny-2)
end do








endif
 
end subroutine boundxy
 