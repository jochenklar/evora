subroutine cosmostart()
  use global
  implicit none

  real(8), external :: dcosmodt

  redshift = zstart
  cosmo = 1 / (1 + zstart)  
  cosmodot = dcosmodt(t,cosmo)

end subroutine cosmostart

subroutine compcosmo()
  use global
  implicit none

  real(8), external :: dcosmodt

  oldcosmo = cosmo

  call rungekutta(cosmo,t,dt,dcosmodt)
  cosmodot = dcosmodt(t,cosmo)
  expand = cosmo * cosmodot
  redshift = 1.0 / cosmo - 1.0

end subroutine compcosmo

subroutine expansion(q)
  use global
  implicit none

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: press
  integer :: ix,iy,iz

   do iz=lz,mz
     do iy=ly,my
        do ix=lx,mx
           call pressure(q(:,ix,iy,iz),press)
           q(5,ix,iy,iz) = q(5,ix,iy,iz) - dt * g35 * expand * press / g1
           q(6,ix,iy,iz) = q(6,ix,iy,iz) - dt * g35 * expand * q(6,ix,iy,iz)
        enddo
     enddo
  enddo

end subroutine expansion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine peebles(fpeebles)
  use global
  implicit none

  ! this routine computes the peebles growth factor dlnD/dlna

  real(8), intent(out) :: fpeebles

  real(8) :: f1,f2,integral

  real(8), external :: peeblesint,adot2

  ! the term with dadot/da
  f1 = (lambda0 * (cosmo**2) - 0.5 * omega0 / cosmo) / adot2(cosmo)
  
  ! the term with dI/da
  call trapezoid(integral,1d-8,cosmo,1d-8,peeblesint)
  f2 = cosmo * peeblesint(cosmo) / integral

  fpeebles = f1 + f2 - 1.0

end subroutine peebles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8) function adot2(a)
  use global
  implicit none
  
  real(8), intent(in) :: a

  adot2 = omega0 * (1.0 / a - 1.0) + lambda0 * ((a**2.0) - 1.0) + 1.0

endfunction adot2

real(8) function dcosmodt(t,a)
  implicit none
  
  real(8), intent(in) :: t,a
  
  real(8), external :: adot2

  dcosmodt = (a**2) * sqrt(adot2(a))

end function dcosmodt

real(8) function peeblesint(a)
  use global
  implicit none
  
  real(8), intent(in) :: a
  
  real(8), external :: adot2

  peeblesint = adot2(a)**(-1.5)

end function peeblesint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rungekutta(f,x,dx,ode)
  implicit none

  real(8), intent(inout) :: f,x
  real(8), intent(in)    :: dx
  real(8), external      :: ode

  real(8) :: xx,dx2,ff,dfdx1,dfdx2,dfdx3

  dfdx1 = ode(x,f)

  dx2 = 0.5 * dx
 
  xx = x + dx2

  ff = f + dx2 * dfdx1
  dfdx2 = ode(xx,ff)

  ff = f + dx2 * dfdx2
  dfdx3 = ode(xx,ff)
  
  ff = f + dx * dfdx3
  dfdx3 = dfdx2 + dfdx3
  dfdx2 = ode(x+dx,f)

  f = f + (dfdx1 + dfdx2 + 2.0 * dfdx3) * dx / 6 

end subroutine rungekutta

subroutine trapezoid(int,a,b,delta,f)
  implicit none

  real(8), intent(out) :: int
  real(8), intent(in) :: a,b,delta
  real(8), external :: f

  real(8) :: x,dx,fl,fr

  int = 0.0
  x = a
  dx = delta
  fr = f(x)

  do while (x .lt. b)
     fl = fr
     if (x + dx .gt. b) then
        dx = b - x
        x = b
     else
        x = x + dx
     endif
     fr = f(x)
     int = int + 0.5 * dx * (fl + fr)
  enddo

end subroutine trapezoid
