subroutine ic(q)
  use global
  implicit none

  real(8), intent(out) ,dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8), parameter :: d = 1.0, templ = 400.0, tempr = 800.0

  real(8) :: x,y,z
  integer :: ix,iy,iz

  real(8) :: e,r,el,er,r0

  el = templ / g1 / temp0 * d
  er = tempr / g1 / temp0 * d

  r0 = 0.5 * (xmax - xmin)

  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)
           
           r = x

           if (r .le. r0) then
              e = el
           else
              e = er
           endif

           e = e

           q(1,ix,iy,iz) = d

           q(2,ix,iy,iz) = 0.0
           q(3,ix,iy,iz) = 0.0
           q(4,ix,iy,iz) = 0.0

           q(5,ix,iy,iz) = e
           q(6,ix,iy,iz) = g1 * e / (d**g1)
        enddo
     enddo
  enddo

end subroutine ic

