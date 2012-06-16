!########################################################################!
!#                                                                      #!
!# themaltest.f90                                                       #!
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

subroutine ic(q)
  use global
  implicit none

  real(8), intent(out) ,dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8), parameter :: d = 1.0, templ = 400.0, tempr = 800.0

  real(8) :: x,y,z
  integer :: ix,iy,iz

  real(8) :: e,r,el,er,r0

  el = templ / g1 / temp0 * d
  er = tempr / g1 / temp0 * d

  r0 = 0.5 * (xmax - xmin)

  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)
           
           r = x

           if (r .le. r0) then
              e = el
           else
              e = er
           endif

           e = e

           q(1,ix,iy,iz) = d

           q(2,ix,iy,iz) = 0.0
           q(3,ix,iy,iz) = 0.0
           q(4,ix,iy,iz) = 0.0

           q(5,ix,iy,iz) = e
           q(6,ix,iy,iz) = g1 * e / (d**g1)
        enddo
     enddo
  enddo

end subroutine ic

