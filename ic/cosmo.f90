subroutine ic(q)
  use global
  implicit none

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: x,y,z
  integer :: i,ix,iy,iz,ikx

  real(8) :: d,p,u,v,w,e,temp,tinit,fpeebles,vscale

  real(8) :: kx,kx2,ky,kz,ampx,ampx2,ampy,ampz,maxmode
  namelist /ic_para/ kx,kx2,ky,kz,ampx,ampx2,ampy,ampz,maxmode
  integer :: ierror,io

  character (len=100) :: line
  integer :: pn

  complex(8), dimension(:,:,:), allocatable :: phi
  real(8),    dimension(:,:)  , allocatable :: phibuf2,phibuf4
  real(8),    dimension(:)    , allocatable :: logpower,logk

  real(8), dimension(3) :: f

  real(8), external :: gaussian
  real(8) :: k,power

  if (parallel) then
     print *,'ERROR: cosmo does not work in parallel mode yet' 
     stop
  endif
  if (.not. gravity) then
     print *,'ERROR: cosmo does not work without gravity' 
     stop
  endif

  ! initialize random generator 
  call initrandom(1)

  ! allocate helper arrays
  allocate(phi(lx:mx,ly:my,lz:mz))

  ! read the parameter
  ampx = 0.0 ; kx = 1.0
  ampx2 = 0.0 ; kx2 = 1.0
  ampy = 0.0 ; ky = 1.0
  ampz = 0.0 ; kz = 1.0
  maxmode = 10000.0

  open (unit=10, file=parameterfile, status="old", iostat=ierror)
  read (10,nml=ic_para)
  close(10)
  
  ampx = ampx * cosmo
  ampx2 = ampx2 * cosmo
  ampy = ampy * cosmo
  ampz = ampz * cosmo
  
  kx = 2 * pi * kx 
  kx2 = 2 * pi * kx2
  ky = 2 * pi * ky
  kz = 2 * pi * kz
  
  tinit = 0.01/cosmo**2

  call peebles(fpeebles)
  vscale = fpeebles * cosmodot / cosmo

  ! read in powerspectrum
  open (unit=10, file="ic/powerspectrum", status="old", iostat=ierror)

  ! count lines
  pn = -1 ; io = 1
  do while (io >= 0)
     pn = pn + 1
     read (10, '(A)', iostat = io) line
  end do
  
  ! allocate arrays
  allocate(logk(pn),logpower(pn))

  ! read stuff
  rewind(10)
  do i=1,pn
     read (10,*) k,power
     logk(i)     = log10(k)
     logpower(i) = log10(power)
  enddo
  close(10)

  ! construct white noise density field
  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           phi(ix,iy,iz) = cmplx(gaussian(),0.0,8)
        enddo
     enddo
  enddo

  ! renorm noise
  phi = phi / sqrt(dble(n3))

  ! make fft
  call fft_c2c_forward(phi)
  
  ! multiply with powerspectrum
  call multipower(phi,logk,logpower,pn,maxmode)

  ! assign amplitudes by hand
  if (nx .gt. 1 .and. ampx .ne. 0.0) then
     phi(2,1,1)  = cmplx(-ampx*0.5,0.0,8)
     phi(nx,1,1) = cmplx(-ampx*0.5,0.0,8)
  endif
  if (ny .gt. 1 .and. ampy .ne. 0.0) then
     phi(1,2,1)  = cmplx(-ampy*0.5,0.0,8)
     phi(1,ny,1) = cmplx(-ampy*0.5,0.0,8)
  endif
  if (nz .gt. 1 .and. ampz .ne. 0.0) then
     phi(1,1,1)  = cmplx(-ampz*0.5,0.0,8)
     phi(1,1,nz) = cmplx(-ampz*0.5,0.0,8)
  endif

  ! make fft backwards
  call fft_c2c_backward(phi)

  ! assign to density
  q(1,:,:,:) = dble(phi) + 1.0

  ! get potential
  phi(lx:mx,ly:my,lz:mz) = cmplx(q(1,lx:mx,ly:my,lz:mz) - 1.0,0.0,8)
  call fft_c2c_forward(phi)
  call multigreen(phi)
  call fft_c2c_backward(phi)

  ! construct velocity field
  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call grav_force(phi,phibuf2,phibuf4,ix,iy,iz,f)

           u = vscale * f(1)
           v = vscale * f(2)
           w = vscale * f(3)

           q(2,ix,iy,iz) = q(1,ix,iy,iz) * u
           q(3,ix,iy,iz) = q(1,ix,iy,iz) * v
           q(4,ix,iy,iz) = q(1,ix,iy,iz) * w
        enddo
     enddo
  enddo

  ! construct energy and entropy
  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           temp = tinit * q(1,ix,iy,iz)**g1

           p = q(1,ix,iy,iz) * temp / temp0

           e = 0.5 * ((q(2,ix,iy,iz)**2)+(q(3,ix,iy,iz)**2)+(q(4,ix,iy,iz)**2)) &
                / q(1,ix,iy,iz) + p / g1

           q(5,ix,iy,iz) = e
           q(6,ix,iy,iz) = p / (q(1,ix,iy,iz)**g1)

        enddo
     enddo
  enddo

  deallocate(phi)

end subroutine ic

subroutine initrandom(seed)
  use global
  implicit none

  integer :: i
  integer :: seed
  integer, dimension(42) :: put_seed

  put_seed = seed * (/ (i,i = 1,42) /) * (mpirank + 1)

  call random_seed(PUT = put_seed)

end subroutine initrandom

real(8) function gaussian()
  use global
  implicit none
 
  real,    dimension(2) :: u
 
  ! Box-Muller transform
  call random_number(u)
  gaussian  = sqrt(- 2.0 * log(dble(u(1)))) * cos(2.0 * pi * dble(u(2)))

end function gaussian
