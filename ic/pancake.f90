!########################################################################!
!#                                                                      #!
!# pancake.f90                                                          #!
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

  real(8), intent(inout), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: x,y,z
  integer :: ix,iy,iz

  real(8) :: d,p,u,v,w,e,temp,tinit,fpeebles,vscale

  real(8) :: kx,kx2,ky,kz,ampx,ampx2,ampy,ampz,maxmode
  namelist /ic_para/ kx,kx2,ky,kz,ampx,ampx2,ampy,ampz,maxmode
  integer :: ierror

  ! read the parameter
  if (master) then
     ampx  = 0.0 ; kx = 1.0
     ampx2 = 0.0 ; kx2 = 1.0
     ampy  = 0.0 ; ky = 1.0
     ampz  = 0.0 ; kz = 1.0
     
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
  endif

  if (parallel) then
     call bcastdouble(kx)
     call bcastdouble(kx2)
     call bcastdouble(ky)
     call bcastdouble(kz)
     call bcastdouble(ampx)
     call bcastdouble(ampx2)
     call bcastdouble(ampy)
     call bcastdouble(ampz)
     call bcastdouble(tinit)
     call bcastdouble(vscale)
  endif

  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)

           d = 1.0 + ampx * cos(kx*x) + ampx2 * cos(kx2*x) &
                + ampy * cos(ky*y) + ampz * cos(kz*z)

           d = d * fb

           u = - vscale * (ampx * sin(kx*x) / kx + ampx2 * sin(kx2*x) / kx2)
           v = - vscale * ampy * sin(ky*y) / ky
           w = - vscale * ampz * sin(kz*z) / kz

           if (jeansfloor) then
              p = jeansfloor0 * d*d
           else
              temp = tinit * d**g1
              p = d * temp / temp0
           endif

           e = 0.5 * ((u**2)+(v**2)+(w**2)) * d + p / g1

           q(1,ix,iy,iz) = d
           q(2,ix,iy,iz) = d * u
           q(3,ix,iy,iz) = d * v
           q(4,ix,iy,iz) = d * w
           q(5,ix,iy,iz) = e
           q(6,ix,iy,iz) = p / (d ** g1)

        enddo
     enddo
  enddo

end subroutine ic
