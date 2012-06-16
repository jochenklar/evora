!########################################################################!
!#                                                                      #!
!# riemann.f90                                                          #!
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

subroutine riemannsolver(wl,wr,f,dir)
  use global
  implicit none

  real(8), dimension(nw):: wl,wr 
  real(8), dimension(nf):: f
  integer :: dir

  if (solver .eq. 1) then
     call hll(wl,wr,f,dir)
  elseif (solver .eq. 2) then
     call hllc(wl,wr,f,dir)
  else
     print *,"ERROR choose proper solver!"
     call cleanstop
  endif
  
end subroutine riemannsolver
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hll(wl,wr,f,dir)
  use global
  implicit none

  ! left and right state as seen by the interface 
  real(8), dimension(nw):: wl,wr 

  real(8), dimension(nf):: f
  integer :: dir

  real(8) :: ar, al, sl, sr
  real(8), dimension(nq) :: ql, qr
  real(8), dimension(nf) :: fl, fr
  
  if (wl(5) <= 0 .or. wr(5) <= 0) then
     print *,'ERROR: hll solver unphysical pressure'
     print *,'wl',wl
     print *,'wr',wr
     stop
  endif

  ! computation of wavespeeds
  al = sqrt(g * wl(5) / wl(1))
  ar = sqrt(g * wr(5) / wr(1))
    
  sl = min(wl(dir+1) - al , wr(dir+1) - ar)
  sr = max(wl(dir+1) + al , wr(dir+1) + ar)

  ! computation of flux
  if (0 .lt. sl) then
     call flux(wl,f,dir)
  else if (0 .gt. sr) then
     call flux(wr,f,dir)
  else if (sl .le. 0 .and. 0 .le. sr) then
     call flux(wl,fl,dir)
     call flux(wr,fr,dir)
     call state(wl,ql)
     call state(wr,qr)
     f(1:nf) = (sr * fl(1:nf) - sl * fr(1:nf) & 
          + sl * sr * (qr(1:nf) - ql(1:nf))) / (sr - sl)
  else
     print *,'ERROR: hll solver'
     stop
  endif

end subroutine hll

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hllc(wl,wr,f,dir)
  use global
  implicit none

  ! left and right state as seen by the interface 
  real(8), dimension(nw):: wl,wr 
  real(8), dimension(nf):: f
  integer :: dir

  real(8) :: ar, al, sl, sr, sm
  real(8), dimension(nq):: ql, qr, qx
  real(8), dimension(nf):: fl, fr
  
  real(8) :: facl,facr

  if (wl(5) <= 0 .or. wr(5) <= 0) then
     print *,'ERROR: hllc solver unphysical pressure'
     print *,'wl',wl
     print *,'wr',wr
     stop
  endif

  ! computation of wavespeeds
  al = sqrt(g * wl(5) / wl(1))
  ar = sqrt(g * wr(5) / wr(1))
  
  sl = min(wl(dir+1) - al , wr(dir+1) - ar)
  sr = max(wl(dir+1) + al , wr(dir+1) + ar)

  facl = wl(1) * (sl - wl(dir+1))
  facr = wr(1) * (sr - wr(dir+1))
  sm = (wr(5) - wl(5) + facl * wl(dir+1) - facr * wr(dir+1))/(facl - facr)

  ! computation of flux
  if (0 .lt. sl) then
     ! left state
     call flux(wl,f,dir)
  else if (0 .gt. sr) then
     ! right state
     call flux(wr,f,dir)
  else if (sl .lt. 0 .and. 0 .lt. sm) then
     ! left of the contact
     call flux(wl,fl,dir)
     call state(wl,ql)
     call statex(wl,qx,sl,sm,dir)
     f(1:nf) = fl(1:nf) + sl * (qx(1:nf) - ql(1:nf))
  else if (sm .le. 0 .and. 0 .lt. sr) then
     ! right of the contact
     call flux(wr,fr,dir)
     call state(wr,qr)
     call statex(wr,qx,sr,sm,dir)
     f(1:nf) = fr(1:nf) + sr * (qx(1:nf) - qr(1:nf))
  else
     print *,'ERROR: hllc solver',sm,mpirank
     stop
  endif

end subroutine hllc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine flux(w,f,dir)
  use global
  implicit none
  
  real(8), dimension(nw):: w
  real(8), dimension(nf):: f
  integer :: dir

  real(8) :: u2
  integer :: i


  u2 = w(2)**2 + w(3)**2 + w(4)**2

  if (dir == 1) then 
     ! x-direction
     f(1) = w(1) * w(2)
     f(2) = w(1) * w(2)**2 + w(5)
     f(3) = w(1) * w(2) * w(3)
     f(4) = w(1) * w(2) * w(4)     
     f(5) = w(2) * (0.5 * w(1) * u2 + gg1 * w(5))
     f(6) = w(2) * w(5) / (w(1)**g1)
  else if (dir == 2) then 
     ! y-direction
     f(1) = w(1) * w(3)
     f(2) = w(1) * w(3) * w(2)
     f(3) = w(1) * w(3)**2 + w(5)
     f(4) = w(1) * w(3) * w(4)
     f(5) = w(3) * (0.5 *  w(1) * u2 + gg1 * w(5))
     f(6) = w(3) * w(5) / (w(1)**g1)
  else if (dir == 3) then
     ! z-direction
     f(1) = w(1) * w(4)
     f(2) = w(1) * w(4) * w(2)
     f(3) = w(1) * w(4) * w(3)
     f(4) = w(1) * w(4)**2 + w(5)
     f(5) = w(4) * (0.5 * w(1) * u2 + gg1 * w(5))
     f(6) = w(4) * w(5) / (w(1)**g1)
  else 
     print *,'ERROR flux computation'
     stop
  endif
  
  if (noncie) then
     do i=7,nf
        f(i) = f(1) * w(i-1)
     enddo
  endif
end subroutine flux

subroutine state(w,q)
  use global
  implicit none
  
  real(8), dimension(nw):: w
  real(8), dimension(nq):: q

  real(8) :: u2
  integer :: i

  ! density
  q(1) = w(1)

  ! momenta
  q(2) = w(1) * w(2)
  q(3) = w(1) * w(3)
  q(4) = w(1) * w(4)

  ! energy
  u2 = (w(2)**2) + (w(3)**2) + (w(4)**2)
  q(5) = 0.5 * w(1) * u2 + w(5) / g1

  ! entropy
  q(6) = w(5) / (w(1)**g1)

  ! number densities
  if (noncie) then
     do i=7,nf
        q(i) = w(1) * w(i-1)
     enddo
  endif

end subroutine state

subroutine statex(w,qx,s,sm,dir)
  ! computes the intermediate state for the hllc riemann solver
  use global
  implicit none
  
  real(8), dimension(nw):: w
  real(8), dimension(nq):: qx
  real(8) :: s,sm
  integer :: dir
  
  real(8) :: energy
  integer :: i

  ! density
  qx(1) = w(1) * (s - w(dir+1)) / (s - sm)

  ! momenta
  if (dir .eq. 1) then
     qx(2) = sm   * qx(1)
     qx(3) = w(3) * qx(1)
     qx(4) = w(4) * qx(1)
  else if (dir .eq. 2) then
     qx(2) = w(2) * qx(1)
     qx(3) = sm   * qx(1)
     qx(4) = w(4) * qx(1)
  else if (dir .eq. 3) then
     qx(2) = w(2) * qx(1)
     qx(3) = w(3) * qx(1)
     qx(4) = sm   * qx(1)
  endif

  ! energy
  energy = 0.5 * w(1) * ((w(2)**2) + (w(3)**2) + (w(4)**2)) &
       + w(5) / g1
  qx(5) = qx(1) * (energy / w(1) + (sm - w(dir+1)) &
       * (sm + (w(5) / (w(1) * (s - w(dir+1))))))

  ! entropy
  qx(6) = qx(1) * w(5) / (w(1)**g)

  ! number densities
  if (noncie) then
     do i=7,nf
        qx(i) = qx(1) * w(i-1)
     enddo
  endif

end subroutine statex
