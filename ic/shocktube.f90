!########################################################################!
!#                                                                      #!
!# shocktube.f90                                                        #!
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

  real(8) :: d,u,p
  integer :: dir

  real(8) :: dl,ul,pl,dr,ur,pr,x_diss
  namelist /ic_para/ dl,ul,pl,dr,ur,pr,x_diss,dir
  integer :: ierror

  ! read the parameter
  if (master) then
     open (unit=10, file=parameterfile, status="old", iostat=ierror)
     read (10,nml=ic_para)
     close(10) 
  endif

  if (parallel) then
     call bcastdouble(dl)
     call bcastdouble(ul)
     call bcastdouble(pl)
     call bcastdouble(dr)
     call bcastdouble(ur)
     call bcastdouble(pr)
     call bcastdouble(x_diss)
  endif

  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           call coordinates(ix,iy,iz,x,y,z)

           if (dir .eq. 1 .and. x < x_diss &
              .or. dir .eq. 2 .and. y < x_diss &
              .or. dir .eq. 3 .and. z < x_diss ) then
              d = dl
              u = ul
              p = pl
           else
              d = dr
              u = ur
              p = pr
           endif

           q(1,ix,iy,iz) = d

           if (dir .eq. 1) then
              q(2,ix,iy,iz) = d * u
              q(3,ix,iy,iz) = 0.0
              q(4,ix,iy,iz) = 0.0
           else if (dir .eq. 2) then
              q(2,ix,iy,iz) = 0.0
              q(3,ix,iy,iz) = d * u
              q(4,ix,iy,iz) = 0.0
           else if (dir .eq. 3) then
              q(2,ix,iy,iz) = 0.0
              q(3,ix,iy,iz) = 0.0
              q(4,ix,iy,iz) = d * u
           else
              print *,"Error choose proper direction",dir
              call cleanstop
           endif

           q(5,ix,iy,iz) = 0.5 * d * (u**2) + p / g1
           q(6,ix,iy,iz) = p / (d ** g1)

        enddo
     enddo
  enddo

end subroutine ic

