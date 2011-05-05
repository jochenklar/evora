subroutine ic(q)
  use global
  implicit none

  ! shockcloud interaction from Ziegler (2005)
  real(8), parameter :: dl=3.86859,ul=0.0,pl=167.345 
  real(8), parameter :: dr=1.0,ur=-11.2536,pr=1.0
  real(8), parameter :: ds=10.0, x_diss=0.1, x_sphere=0.3

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: x,y,z
  integer :: ix,iy,iz

  real(8) :: d,u,p,r

  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)

           r = sqrt((x-x_sphere)**2 + y**2 + z**2)

           if (x <= x_diss) then
              d = dl
              u = ul
              p = pl
           else
              d = dr
              if (r <= 0.15) d = ds
              u = ur
              p = pr
           endif

           q(1,ix,iy,iz) = d 
           q(2,ix,iy,iz) = d * u
           q(3,ix,iy,iz) = 0.0
           q(4,ix,iy,iz) = 0.0
           q(5,ix,iy,iz) = .5 * d * u ** 2 + p / (g - 1)
           q(6,ix,iy,iz) = p / (d ** (g - 1))
        enddo
     enddo
  enddo

end subroutine ic
