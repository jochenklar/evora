subroutine ic(q)
  use global
  implicit none

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: x,y,z
  integer :: ix,iy,iz

  real(8) :: d,p,u,v,w,e,temp,tinit

  real(8) :: kx,kx2,ampx,ampx2
  namelist /ic_para/ kx,kx2,ampx,ampx2
  integer :: ierror

! read the parameter
  if (master) then
     ampx  = 2.0 ; kx = 1.0
     ampx2 = 0.0 ; kx2 = 1.0

     open (unit=10, file=parameterfile, status="old", iostat=ierror)
     read (10,nml=ic_para)
     close(10)
 
     ampx = ampx * cosmo
     ampx2 = ampx2 * cosmo

     kx = 2 * pi * kx 
     kx2 = 2 * pi * kx2

     tinit = 0.01/cosmo**2
  endif

  if (parallel) then
     call bcastdouble(kx)
     call bcastdouble(kx2)
     call bcastdouble(ampx)
     call bcastdouble(ampx2)
     call bcastdouble(tinit)
  endif

  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)

           d = 1 + ampx * cos(kx*x) + ampx2 * cos(kx2*x)

           u = - sqrt(cosmo) * &
                (ampx * sin(kx*x) / kx + ampx2 * sin(kx2*x) / kx2)

           temp = tinit * d**g1

           p = d * temp  / temp0

           e = 0.5 * u**2 * d + p / g1

           q(1,ix,iy,iz) = d
           q(2,ix,iy,iz) = d * u
           q(3,ix,iy,iz) = 0.0
           q(4,ix,iy,iz) = 0.0
           q(5,ix,iy,iz) = e
           q(6,ix,iy,iz) = p / (d ** g1)

        enddo
     enddo
  enddo

end subroutine ic
