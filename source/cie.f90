!########################################################################!
!#                                                                      #!
!# cie.f90                                                              #!
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

subroutine cie(rho,press,n)
  use global
  implicit none

  real(8), intent(in)                   :: rho
  real(8), intent(inout)                :: press
  real(8), intent(inout), dimension(ns) :: n

  real(8), dimension(nchem) :: chemrate
  real(8), dimension(ncool) :: coolrate

  real(8) :: temp,lambda,gamma,coolio,cooltime
  real(8) :: ct,cdt

  if (cooling) then
     ct = 0.0
     do while (ct .lt. dt)
        call press2temp(press,rho,n,temp)
        call comprates(temp,chemrate,coolrate)
        call cieabundance(n,chemrate)
        call cool(n,coolrate,lambda,gamma)

        coolio = gamma - lambda
        cooltime = abs(press / coolio)
        cdt = min(c_chem * cooltime,dt-ct,dt*substep) 
        ct = ct + cdt
        press = (press + cdt * gamma)/(1.0 + cdt * lambda / press)
     enddo
  else
     call press2temp(press,rho,n,temp)
     call comprates(temp,chemrate,coolrate)
     call cieabundance(n,chemrate)
  endif

end subroutine cie

subroutine cieabundance(n,chemrate)
  use global
  implicit none
 
  real(8), intent(inout), dimension(ns) :: n
  real(8), intent(in), dimension(nchem) :: chemrate

  real(8) :: fracH,fracHe1,fracHe2

  if (background .and. uv .gt. 0.0) then
     call implicitcie(n,chemrate)
  else  
     fracH   = chemrate(1)/chemrate(4)
     fracHe1 = chemrate(2)/chemrate(5)
     fracHe2 = chemrate(3)/chemrate(6)

     n(1) = nH / (1.0 + fracH)
     n(2) = n(1) * fracH
     n(3) = nHe / (1.0 + fracHe1 * (1.0 + fracHe2))  
     n(4) = n(3) * fracHe1
     n(5) = n(4) * fracHe2
     n(6) = n(2) + n(4) + 2 * n(5)
  endif
     
  call ncheck(n)
  n = max(n,nfloor)

end subroutine cieabundance

subroutine implicitcie(n,chemrate)
  use global
  implicit none 

  real(8), intent(inout), dimension(ns) :: n
  real(8), intent(in), dimension(nchem) :: chemrate

  real(8), parameter :: cietol=1d-10
  integer, parameter :: nit =100

  real(8) :: ne,nHII,nHeII,nHeIII

  call rf_cie(ne,nfloor,nH + 2 * nHe,chemrate)

  call cieions(nHII,nHeII,nHeIII,ne,chemrate)

  n(1) = nH - nHII
  n(2) = nHII
  n(3) = nHe - nHeII - nHeIII
  n(4) = nHeII
  n(5) = nHeIII
  n(6) = ne

  n = max(n,nfloor)

end subroutine implicitcie

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine rf_cie(m,in1,in2,chemrate)
  use global
  implicit none

  real(8), intent(inout)                 :: m
  real(8), intent(in)                    :: in1,in2 
  real(8), intent(in) , dimension(nchem) :: chemrate

  real(8) :: l,r,fr,fl,fm
  integer :: i,s

  real(8), external :: fne ! external function !

  l = in1 ; r = in2
  fl = fne(l,chemrate) ; fr = fne(r,chemrate)
 
  if (fl * fr .gt. 0.0) then
     print *,"ERROR no root between in1 and in2",l,r,fl,fr
     call cleanstop
  endif

  s = 0
  do i=1,cie_nit
     m = (fr * l - fl * r)/(fr - fl)
     if (abs(r-l) .lt. cie_tol*abs(r+l)) exit
     fm = fne(m,chemrate)
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
     if (i .eq. cie_nit) then
        print *,'ERROR cie w background does not converge!'
        call cleanstop
     endif
  enddo

end subroutine rf_cie

real(8) function fne(ne,chemrate)
  use global
  implicit none

  real(8), intent(in)                   :: ne
  real(8), intent(in), dimension(nchem) :: chemrate

  real(8) :: nHII,nHeII,nHeIII

  call cieions(nHII,nHeII,nHeIII,ne,chemrate)

  fne = nHII + nHeII + 2 * nHeIII - ne

end function fne

subroutine cieions(nHII,nHeII,nHeIII,ne,chemrate)
  use global
  implicit none

  real(8), intent(out)                  :: nHII,nHeII,nHeIII
  real(8), intent(in)                   :: ne
  real(8), intent(in), dimension(nchem) :: chemrate

  real(8) :: ion1,ion2,ion3,xi1,xi2,xi3

  ! ionisation
  ion1 = chemrate(1) * ne + chemrate(7)
  ion2 = chemrate(2) * ne + chemrate(8)
  ion3 = chemrate(3) * ne + chemrate(9)

  if (ion1 .lt. ratefloor) then
     nHII = nfloor
  else
     xi1 = chemrate(4) * ne / ion1
     nHII = nH / (1.0 + xi1)
  endif

  if (ion2 .lt. ratefloor) then
     nHeII = nfloor
     nHeIII = nfloor
  else
     xi2 = chemrate(5) * ne / ion2
     if (ion3 .lt. ratefloor) then
        nHeIII = nfloor
        nHeII = nHe / (1.0 + xi2)
     else
        xi3 = chemrate(6) * ne / ion3
        nHeIII = nHe / (1.0 + (xi2 + 1.0) * xi3)
        nHeII = nHe / (1.0 + xi2 + 1.0 / xi3)
     endif
 endif

end subroutine cieions
