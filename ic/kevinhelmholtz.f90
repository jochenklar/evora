subroutine ic(q)
  use global
  implicit none

  real(8), intent(out) ,dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: x,y,z,r
  integer :: ix,iy,iz

  real(8) :: d,u,v,p
  integer :: dir

  integer :: ierror

  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)

           r = abs(y - 0.5)

           if (r .le. 0.25) then
              d = 2.0
              u = 0.5
           else
              d = 1.0
              u = -0.5
           endif
           
           v = 0.1 * sin(4 * pi * x) &
                * (exp(- ((y - 0.25)**2) / 0.05) &
                + exp(- ((y - 0.75)**2) / 0.05))

           p = 2.5
           
           q(1,ix,iy,iz) = d
           q(2,ix,iy,iz) = d * u
           q(3,ix,iy,iz) = d * v
           q(4,ix,iy,iz) = 0
           q(5,ix,iy,iz) = .5 * d * ((u**2)+(v**2)) + p / (g - 1)
           q(6,ix,iy,iz) = p / (d ** (g - 1))
        enddo
     enddo
  enddo

end subroutine ic

