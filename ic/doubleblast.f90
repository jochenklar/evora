subroutine ic(q)
  use global
  implicit none

  real(8), intent(out) ,dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: x,y,z
  integer :: ix,iy,iz

  real(8) :: d,u,p

  integer :: ierror

  !  tend  = 0.038

  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)

           if (x .lt. 0.1) then
              p = 1000.0
           elseif (x .gt. 0.9) then
              p = 100.0
           else
              p = 0.01
           endif

           q(1,ix,iy,iz) = 1.0
           q(2,ix,iy,iz) = 0
           q(3,ix,iy,iz) = 0
           q(4,ix,iy,iz) = 0
           q(5,ix,iy,iz) = p / (g - 1)
           q(6,ix,iy,iz) = p / (d ** (g - 1))
        enddo
     enddo
  enddo

end subroutine ic

