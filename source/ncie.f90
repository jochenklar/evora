!########################################################################!
!#                                                                      #!
!# ncie.f90                                                             #!
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

subroutine ncie(rho,press,n)
  use global
  implicit none   

  real(8), intent(in)                   :: rho
  real(8), intent(inout)                :: press
  real(8), intent(inout), dimension(ns) :: n

  real(8), dimension(nchem) :: chemrate
  real(8), dimension(ncool) :: coolrate

  real(8), dimension(ns1) :: delta
  real(8), dimension(3) :: ion,rec

  real(8) :: temp,lambda,gamma,coolio,cooltime
  real(8) :: ct,cdt,cdtrho

  ct = 0.0

  do while (ct .lt. dt)

     ! compute rates
     if (cooling .or. ct .eq. 0.0) then
        call press2temp(press,rho,n,temp)
        call comprates(temp,chemrate,coolrate)
     endif
          
     if (cooling) then
        ! coolingtime computed via cooling function
        call cool(n,coolrate,lambda,gamma)
        coolio = gamma - lambda
        cooltime = abs(press / coolio)
     else
        ! chemical time computed via chemical network
        ! computing explicit differnces for timestep
        ion(1) = (chemrate(1) * n(6) + chemrate(7)) * n(1)
        ion(2) = (chemrate(2) * n(6) + chemrate(8)) * n(3)
        ion(3) = (chemrate(3) * n(6) + chemrate(9)) * n(4)
        rec(1) = chemrate(4) * n(2) * n(6)
        rec(2) = chemrate(5) * n(4) * n(6)
        rec(3) = chemrate(6) * n(5) * n(6)
        delta(1) = rec(1) - ion(1)
        delta(2) = ion(1) - rec(1)
        delta(3) = rec(2) - ion(2)
        delta(4) = ion(2) - rec(2) - ion(3) + rec(3)
        delta(5) = ion(3) - rec(3)
        cooltime = maxval(n(1:ns1) / abs(delta))
     endif

     ! compute timestep
     cdt = min(dt-ct,c_chem * cooltime,dt*substep)

     ! update chemical timestep
     ct = ct + cdt

     ! compute cdt times (?) rho
     cdtrho = cdt/rho

     ! calling modified patankar routines
     call patankarH(n(1:2),n(6),cdt,chemrate)
     call patankarHe(n(3:5),n(6),cdt,chemrate)
     n(6) = n(2) + n(4) + 2.0 * n(5) 
     call ncheck(n)
     n = max(n,nfloor)

     ! update temperature
     if (cooling) &
          press = (press + cdt * gamma)/(1.0 + cdt * lambda / press)
  enddo

end subroutine ncie

subroutine patankarH(n,ne,cdt,chemrate)
  use global
  implicit none

  real(8), intent(inout), dimension(1:2)    :: n
  real(8), intent(in)                       :: ne,cdt
  real(8), intent(in),    dimension(nchem)  :: chemrate

  real(8) :: d1inv,d2,p12,p21
  real(8), dimension(1:2) :: x

  p12 = chemrate(4) * ne * cdt
  p21 = (chemrate(1) * ne + chemrate(7)) * cdt

  d1inv = 1.0 / (1.0 + p21)
  d2 = 1.0 + p12
  
  x(2) = (n(2) + p21 * n(1) * d1inv) / (d2 - p21 * p12 * d1inv)
  x(1) = (n(1) + p12 * x(2)) * d1inv

  n(1:2) = x(1:2)

end subroutine patankarH

subroutine patankarHe(n,ne,cdt,chemrate)
  use global
  implicit none

  real(8), intent(inout), dimension(3:5)    :: n
  real(8), intent(in)                       :: ne,cdt
  real(8), intent(in),    dimension(nchem)  :: chemrate

  real(8) :: d3inv,d4,d5inv,p34,p43,p45,p54
  real(8), dimension(3:5) :: x

  p34 = chemrate(5) * ne * cdt
  p43 = (chemrate(2) * ne + chemrate(8)) * cdt
  p45 = chemrate(6) * ne * cdt
  p54 = (chemrate(3) * ne + chemrate(9)) * cdt

  d3inv = 1.0 / (1.0 + p43)
  d4 = 1.0 + p34 + p54
  d5inv = 1.0 / (1.0 + p45)

  x(4) = (n(4) + p43 * n(3) * d3inv + p45 * n(5) * d5inv ) &
       / (d4 - p43 * p34 * d3inv - p45 * p54 * d5inv)
  x(3) = (n(3) + p34 * x(4)) * d3inv
  x(5) = (n(5) + p54 * x(4)) * d5inv

  n(3:5) = x(3:5)

end subroutine patankarHe
