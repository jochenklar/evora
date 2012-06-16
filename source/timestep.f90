!########################################################################!
!#                                                                      #!
!# timestep.f90                                                         #!
!# Copyright (C) 2007-2012 Jochen Klar                                  #!
!#                                                                      #!
!# This file is part of evora - a hydrodynamic code for the computation #!
!# of the cosmological evolution of a polytropic fluid including the    #!
!# influence of gravity, primordial chemical processes, radiative       #!
!# cooling, heating by a UV background, and thermal conduction.         #!
!#                                                                      #!
!# evora is free software: you can redistribute it and/or modify        #!
!# it under the terms of the GNU General Public License as published by #!
!# the Free Software Foundation, either version 3 of the License, or    #!
!# (at your option) any later version.                                  #!
!#                                                                      #!
!# evora is distributed in the hope that it will be useful,             #!
!# but WITHOUT ANY WARRANTY; without even the implied warranty of       #!
!# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #!
!# GNU General Public License for more details.                         #!
!#                                                                      #!
!# You should have received a copy of the GNU General Public License    #!
!# along with evora.  If not, see <http://www.gnu.org/licenses/>.       #!  
!#                                                                      #!
!########################################################################!

subroutine hydrostep(w)
  use global
  implicit none

  real(8), intent(in), dimension(nw,lx2:mx2,ly2:my2,lz2:mz2) :: w
     
  real(8) :: v,vx,vy,vz,a
  integer :: ix,iy,iz

  v = 1e6

  do iy=ly,my
     do ix=lx,mx
        do iz=lz,mz
           a = sqrt(g * w(5,ix,iy,iz)/w(1,ix,iy,iz))
           vx = dx / (abs(w(2,ix,iy,iz)) + a)
           vy = dy / (abs(w(3,ix,iy,iz)) + a)
           vz = dz / (abs(w(4,ix,iy,iz)) + a)
           v = min(vx,vy,vz,v)
        enddo
     enddo
  enddo

  dt = min(cfl * v,dt)

  if (dt .le. 0.0) THEN
     print *,'ERROR: hydrostep is negative'
     print *,'infinite pressure ?'
     call cleanstop
  endif

end subroutine hydrostep

subroutine cosmostep()
  use global
  implicit none

  dt = min(cosmo / cosmodot * c_cosmo,dt)

  if (dt .le. 0.0) THEN
     print *,'ERROR: cosmostep is negative'
     call cleanstop
  endif

end subroutine cosmostep

subroutine thermalstep(q,w)
  use global
  implicit none

  real(8), intent(in), dimension(nq,lx:mx,ly:my,lz:mz) :: q
  real(8), intent(in), dimension(nw,lx2:mx2,ly2:my2,lz2:mz2) :: w

  real(8) :: v,press,temp,kappa
  integer :: ix,iy,iz

  v = 1d99

  do iz=lz,mz
     do iy=ly,my 
        do ix=lx,mx
           call pressure(q(:,ix,iy,iz),press)
           call press2temp(press,q(1,ix,iy,iz),q(7:nq,ix,iy,iz),temp)
           kappa = kappa0 * (temp**2.5)
           v = min(w(1,ix,iy,iz)/kappa,v)
        enddo
     enddo
  enddo

  v = c_cond * (dx**2.0) * v / temp0

  dt = min(v,dt)

  if (dt .le. 0.0) THEN
     print *,'ERROR: thermalstep is negative'
     call cleanstop
  endif

end subroutine thermalstep

subroutine coolstep(q,w)
  use global
  implicit none

  real(8), intent(in), dimension(nq,lx:mx,ly:my,lz:mz) :: q
  real(8), intent(in), dimension(nw,lx2:mx2,ly2:my2,lz2:mz2) :: w

  real(8), dimension(nchem) :: chemrate
  real(8), dimension(ncool) :: coolrate
  
  real(8) :: v,temp,lambda,gamma,coolio,cooltime

  integer :: ix,iy,iz

  v = 1e6

  do iz=lz,mz
     do iy=ly,my
        do ix=lx,mx
           call press2temp(w(5,ix,iy,iz),q(1,ix,iy,iz),q(7:nq,ix,iy,iz),temp)
           call comprates(temp,chemrate,coolrate)
           call cieabundance(q(7:nq,ix,iy,iz),chemrate,q(1,ix,iy,iz))
           call cool(q(7:nq,ix,iy,iz),coolrate,lambda,gamma)
           coolio = gamma - lambda
           cooltime = abs(w(5,ix,iy,iz) / coolio)
           v = min(v,cooltime)
        enddo
     enddo
  enddo

  dt = min(dt,c_cooldt * v)

  if (dt .le. 0.0) THEN
     print *,'ERROR: coolstep is negative'
     call cleanstop
  endif

endsubroutine coolstep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine updatetime(go,dump,k)
  use global
  implicit none

  logical :: go,dump
  integer :: k

  if (supercomoving) then
     if (redshift .le. zend) then
        call timestepiteration(zend)
        go = .false.
     else if (k .le. ndump .and. redshift .lt. dump_redshift(k)) then
        call timestepiteration(dump_redshift(k))
        dump = .true.
     endif
  else
     if (tend - t .le. dt) then
        dt = tend - t
        go = .false.
     endif
  endif

  t = t + dt

endsubroutine updatetime

subroutine timestepiteration(targetredshift)
  use global
  implicit none

  real(8), intent(in) :: targetredshift

  integer, parameter :: cosmoend_nit=10
  real(8), parameter :: cosmoend_tol=1d-20

  real(8) :: targetcosmo

  real(8) ::l,r,fl,fr,m,fm

  integer :: i,s

  targetcosmo = 1.0 / (1.0 + targetredshift)

  l = 0.0
  r = dt
  fl = targetcosmo - oldcosmo
  fr = targetcosmo - cosmo

  s = 0
  do i=1,cosmoend_nit
     m = (fr * l - fl * r)/(fr - fl)

     if (abs(r-l) .lt. cosmoend_tol*abs(r+l)) exit

     dt = m
     cosmo = oldcosmo
     call compcosmo

     fm = targetcosmo - cosmo

     if (abs(fm) .lt. cosmoend_tol) exit

     if (fm * fr .gt. 0.0) then
        r = m ; fr = fm
        if (s .eq. -1) fl = fl * 0.5
        s = -1
     else if (fl * fm .gt. 0.0) then
        l = m ; fl= fm
        if (s .eq. 1) fr = fr * 0.5
        s = 1
     else 
        exit
     endif
     if (i .eq. cosmoend_nit) then
        print *,'ERROR timestep iteration does not converge!'
        call cleanstop
     endif
  enddo

end subroutine timestepiteration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
