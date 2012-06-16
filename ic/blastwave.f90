!########################################################################!
!#                                                                      #!
!# blastwave.f90                                                        #!
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

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: x,y,z 
  integer :: ix,iy,iz

  real(8) :: e,r,dv
  integer :: counter

  real(8) :: r0,d0,p0,e0
  namelist /ic_para/ r0,d0,p0,e0
  integer :: ierror

  ! read the parameter
  if (master) then
     open (unit=10, file=parameterfile, status="old", iostat=ierror)
     read (10,nml=ic_para)
     close(10) 
  endif

  call bcastdouble(r0)
  call bcastdouble(d0)
  call bcastdouble(p0)
  call bcastdouble(e0)

  counter = 0

  do iz=lz,mz
     do iy=ly,my
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)

           r = sqrt(x**2 + y**2 + z**2)

           if (r < r0) then
              e = e0
              counter = counter + 1
           else
              e = 0.0
           endif

           q(1,ix,iy,iz) = d0
           q(2,ix,iy,iz) = 0
           q(3,ix,iy,iz) = 0
           q(4,ix,iy,iz) = 0
           q(5,ix,iy,iz) = e
           q(6,ix,iy,iz) = e * g1 / (d0 ** g1)
        enddo
     enddo
  enddo

  call allreducesumint(counter)
  dv = dx*dy*dz*counter

  do iz=lz,mz
     do iy=ly,my
        do ix=lx,mx
           if (q(5,ix,iy,iz) > 0.0) then
              q(5,ix,iy,iz) = q(5,ix,iy,iz) / dv
              q(6,ix,iy,iz) = q(6,ix,iy,iz) / dv
           else
              q(5,ix,iy,iz) = p0 / g1
              q(6,ix,iy,iz) = p0 / (d0 ** g1)
           endif
        enddo
     enddo
  enddo

end subroutine ic

