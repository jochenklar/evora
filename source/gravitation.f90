subroutine gravpotential(q,phi)
  use global
  implicit none

  real(8),    intent(in),  dimension(nq,lx:mx,ly:my,lz:mz) :: q
  complex(8), intent(out), dimension(lx:mx,ly:my,lz:mz)    :: phi

  real(8) :: meanrho

  ! compute the source term
  if (supercomoving) then
     phi(lx:mx,ly:my,lz:mz) = cosmo * grav * (q(1,lx:mx,ly:my,lz:mz)/fb - 1.0)
  else
     meanrho = sum(q(1,lx:mx,ly:my,lz:mz))
     if (parallel) call allreducesumdouble(meanrho)
     meanrho = meanrho / n3
     phi(lx:mx,ly:my,lz:mz) = grav * (q(1,lx:mx,ly:my,lz:mz) - meanrho)
  endif
  
  ! forward transformation
  call fft_c2c_forward(phi)

  ! multiply with greens function
  call multigreen(phi)

  ! backward transformation
  call fft_c2c_backward(phi)

end subroutine gravpotential

subroutine multigreen(phi)
  use global
  implicit none

  complex(8), intent(inout), dimension(lx:mx,ly:my,lz:mz) :: phi

  integer :: ix,iy,iz
  real(8) :: green,sin2x,sin2y,sin2z,kx,ky,kz

  do iz=lz,mz
     call kvalue(iz,nz,nyqz,kz)
     sin2z = sin(kfz2 * kz)**2

     do iy=ly,my
        call kvalue(iy,ny,nyqy,ky)
        sin2y  = sin(kfy2 * ky)**2

        do ix=lx,mx
           call kvalue(ix,nx,nyqx,kx)
           sin2x = sin(kfx2 * kx)**2

           if (ix .eq. 1 .and. iy .eq. 1 .and. iz .eq. 1) then
              green = 0.0
           else 
              green = - h2 / (sin2x+sin2y+sin2z)
           endif

           phi(ix,iy,iz) = phi(ix,iy,iz) / dble(n3) * dble(green)

        enddo
     enddo
  enddo
end subroutine multigreen

subroutine multipower(phi,logk,logpower,pn,maxmode)
  use global
  implicit none

  complex(8), intent(inout), dimension(lx:mx,ly:my,lz:mz) :: phi
  integer, intent(in) :: pn
  real(8), intent(in) :: maxmode
  real(8), intent(in), dimension(pn) :: logk,logpower

  integer :: ix,iy,iz,i
  real(8) :: kx,ky,kz,kx2,ky2,kz2,logkk,pk,delta,kf,dpk,logpk,kk,mink

  ! fundamental mode in Mpc/h
  kf = 2.0 * pi / box

  ! spaceing in log k space
  delta = (logk(pn) - logk(1)) / (pn-1)

  ! derivative for large k
  dpk = (logpower(pn) - logpower(pn-3))/(logk(pn) - logk(pn-3))

  ! mink
  mink = 2 * pi / maxmode

  do iz=lz,mz
     call kvalue(iz,nz,nyqz,kz)
     kz2 = kz * kz
     do iy=ly,my
        call kvalue(iy,ny,nyqy,ky)
        ky2 = ky * ky
        do ix=lx,mx
           call kvalue(ix,nx,nyqx,kx)
           kx2 = kx * kx
           if (ix .eq. 1 .and. iy .eq. 1 .and. iz .eq. 1) then
              phi(ix,iy,iz) = 0.0
           else 
              ! kompute k for cell
              kk = sqrt(kx2 + ky2 + kz2) * kf

              if (kk .ge. mink) then
                 ! get k space index
                 logkk = log10(kk)
                 i = floor((logkk - logk(1)) / delta) + 1

                 ! get corresponding powerspectrum
                 if (i .le. pn) then
                    ! interpolate
                    logpk = logpower(i) * &
                         (logk(i+1) - logkk) / (logk(i+1) - logk(i)) &
                         + logpower(i+1) * &
                         (logkk - logk(i)) / (logk(i+1) - logk(i))
                 else
                    ! extrapolate
                    logpk = logpower(pn) + (logkk - logk(pn)) * dpk
                 endif
                 
                 pk = sqrt(10**logpk)
              else
                 ! null components
                 pk = 0.0
              endif

              ! multiply with powerspectrum
              phi(ix,iy,iz) = phi(ix,iy,iz) * pk
           endif
        enddo
     enddo
  enddo

  ! re-norm density 
  if (onedimension) then
     phi = phi / box**0.5
  else if (twodimension) then 
     phi = phi / box
  else
     phi = phi / box**1.5
  endif

end subroutine multipower

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gravintegration(q,phi,phibuf2,phibuf4)
  use global
  implicit none

  real(8),    intent(inout), dimension(nq,lx:mx,ly:my,lz:mz) :: q
  complex(8), intent(in),    dimension(lx:mx,ly:my,lz:mz)    :: phi
  real(8),    intent(in),    dimension(lx:mx,ly:my)          :: phibuf2,phibuf4

  real(8), dimension(3) :: f
  integer :: ix,iy,iz

  do iz=lz,mz
     do iy=ly,my
        do ix=lx,mx
           call grav_force(phi,phibuf2,phibuf4,ix,iy,iz,f)

           q(5,ix,iy,iz) = q(5,ix,iy,iz) &
                + dt * ( q(2,ix,iy,iz) * f(1) &
                + q(3,ix,iy,iz) * f(2) &
                + q(4,ix,iy,iz) * f(3) )

           q(2:4,ix,iy,iz) = q(2:4,ix,iy,iz) + dt * q(1,ix,iy,iz) * f(1:3)          
        enddo
     enddo
  enddo

end subroutine gravintegration

subroutine grav_force(phi,phibuf2,phibuf4,ix,iy,iz,f)
  use global
  implicit none

  complex(8), intent(in),  dimension(lx:mx,ly:my,lz:mz) :: phi
  real(8),    intent(in),  dimension(lx:mx,ly:my)       :: phibuf2,phibuf4
  integer,    intent(in)                                :: ix,iy,iz
  real(8),    intent(out), dimension(3) :: f

  if (nx .eq. 1) then 
     f(1) = 0
  else
     if (ix .eq. lx) then
        f(1) = - (dble(phi(lx0,iy,iz)) - dble(phi(mx,iy,iz))) * dx2inv
     else if (ix .eq. mx) then
        f(1) = - (dble(phi(lx,iy,iz)) - dble(phi(mx0,iy,iz))) * dx2inv
     else
        f(1) = - (dble(phi(ix+1,iy,iz)) - dble(phi(ix-1,iy,iz))) * dx2inv
     endif
  endif

  if (ny .eq. 1) then 
     f(2) = 0
  else
     if (iy .eq. ly) then
        f(2) = - (dble(phi(ix,ly0,iz)) - dble(phi(ix,my,iz))) * dy2inv
     else if (iy .eq. my) then
        f(2) = - (dble(phi(ix,ly,iz)) - dble(phi(ix,my0,iz))) * dy2inv
     else
        f(2) = - (dble(phi(ix,iy+1,iz)) - dble(phi(ix,iy-1,iz))) * dy2inv
     endif
  endif

  if (nz .eq. 1) then
     f(3) = 0
  else
     if (parallel) then
        if (iz .eq. lz) then
           f(3) = - (dble(phi(ix,iy,lz0)) - phibuf2(ix,iy)) * dz2inv
        else if (iz .eq. mz) then
           f(3) = - (phibuf4(ix,iy) - dble(phi(ix,iy,mz0))) * dz2inv
        else
           f(3) = - (dble(phi(ix,iy,iz+1)) - dble(phi(ix,iy,iz-1))) * dz2inv
        endif
     else
        if (iz .eq. lz) then
           f(3) = - (dble(phi(ix,iy,lz0)) - dble(phi(ix,iy,mz))) * dz2inv
        else if (iz .eq. mz) then
           f(3) = - (dble(phi(ix,iy,lz)) - dble(phi(ix,iy,mz0))) * dz2inv
        else
           f(3) = - (dble(phi(ix,iy,iz+1)) - dble(phi(ix,iy,iz-1))) * dz2inv
        endif
     endif
  endif
end subroutine grav_force
   
