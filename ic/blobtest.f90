!########################################################################!
!#                                                                      #!
!# blobtest.f90                                                         #!
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

  real(8) :: d,u,v,w,p,r

  real(8) :: d0,u0,v0,w0,p0,r1,x1,y1,z1,d1,u1,v1,w1,p1
  namelist /ic_para/ d0,u0,v0,w0,p0,r1,x1,y1,z1,d1,u1,v1,w1,p1 
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

  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)

           r = sqrt((x-x1)**2 + (y-y1)**2 + (z-z1)**2)

           if (r <= r1) then
              d = d1
              u = u1
              v = v1
              w = w1
              p = p1
           else
              d = d0
              u = u0
              v = v0
              w = w0
              p = p0
           endif
          
           q(1,ix,iy,iz) = d 
           q(2,ix,iy,iz) = d * u
           q(3,ix,iy,iz) = d * v
           q(4,ix,iy,iz) = d * w
           q(5,ix,iy,iz) = .5 * d * (u**2+v**2+w**2) + p / (g - 1)
           q(6,ix,iy,iz) = p / (d ** (g - 1))
        enddo
     enddo
  enddo

end subroutine ic
