subroutine ic(q)
  use global
  implicit none

  logical, parameter :: noise = .false.

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: x,y,z
  integer :: ix,iy,iz

  real(8) :: d,u,v,w,p,r,e,c

  real(8) :: d0,u0,v0,w0,p0
  real(8) :: r1,x1,y1,z1,d1,u1,v1,w1,p1
  real(8) :: r2,x2,y2,z2,d2,u2,v2,w2,p2 
  namelist /ic_para/ d0,u0,v0,w0,p0,r1,x1,y1,z1,d1,u1,v1,w1,p1, &
       r2,x2,y2,z2,d2,u2,v2,w2,p2
  integer :: ierror
  
  ! read the parameter
  if (master) then
     open (unit=10, file=parameterfile, status="old", iostat=ierror)
     read (10,nml=ic_para)
     close(10) 
  endif

  call bcastdouble(d0)
  call bcastdouble(u0)
  call bcastdouble(v0)
  call bcastdouble(w0)
  call bcastdouble(p0)
  call bcastdouble(r1)
  call bcastdouble(x1)
  call bcastdouble(y1)
  call bcastdouble(z1)
  call bcastdouble(d1)
  call bcastdouble(u1)
  call bcastdouble(v1)
  call bcastdouble(w1)
  call bcastdouble(p1)
  call bcastdouble(r2)
  call bcastdouble(x2)
  call bcastdouble(y2)
  call bcastdouble(z2)
  call bcastdouble(d2)
  call bcastdouble(u2)
  call bcastdouble(v2)
  call bcastdouble(w2)
  call bcastdouble(p2)

  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)

           d = d0
           p = p0
           u = u0
           v = v0 
           w = w0

           r = sqrt((x-x1)**2 + (y-y1)**2 + (z-z1)**2)

           if (r <= r1) then
              d = d1
              u = u1
              v = v1
              w = w1
              p = p1
           endif

           r = sqrt((x-x2)**2 + (y-y2)**2 + (z-z2)**2)

           if (r <= r2) then
              d = d2
              u = u2
              v = v2
              w = w2
              p = p2
           endif              
          
           if (noise) then
              c = sqrt(g * p / d)
              call noiseomator(u,v,w,c,x,y,z)   
           endif

           e = 0.5 * (u**2+v**2+w**2) * d + p / g1

           q(1,ix,iy,iz) = d 
           q(2,ix,iy,iz) = d * u
           q(3,ix,iy,iz) = d * v
           q(4,ix,iy,iz) = d * w
           q(5,ix,iy,iz) = e
           q(6,ix,iy,iz) = (p / d ** (g - 1))
        enddo
     enddo
  enddo

end subroutine ic


subroutine noiseomator(u,v,w,c,x,y,z)
  use global
  implicit none
  
  real(8) :: c,d,u,v,w,x,y,z

  real(8), parameter :: funk=1.0d-2

  real(8) :: noise,k

  box = xmax - xmin

  k = 2.0 * pi / 4.0 * box

  noise = funk * c * (sin(k*x) * sin(k*y) * sin(k*z))

  u = u + noise
  v = v + noise
  w = w + noise

end subroutine noiseomator
