!########################################################################!
!#                                                                      #!
!# chemistry.f90                                                        #!
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

subroutine initchem(q)
  use global
  implicit none

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: press
  integer :: i,ix,iy,iz

  if (master) then
     chiH = rhoH / mH
     chiHe = rhoHe / mHe
     call compcoef()
  endif

  call bcastdouble(chiH)
  call bcastdouble(chiHe)
  call bcastcoef()

  do iz=lz,mz
     do iy=ly,my
        do ix=lx,mx

           ! compute number densities of H and He in cell
           nH  = chiH  * q(1,ix,iy,iz)
           nHe = chiHe * q(1,ix,iy,iz)

           ! compute pressure
           call pressure(q(:,ix,iy,iz),press)

           ! fully neutral
           do i=7,nq
              if (i .eq. 7) then
                 q(i,ix,iy,iz) = nH
              else if (i .eq. 9) then
                 q(i,ix,iy,iz) = nHe
              else               
                 q(i,ix,iy,iz) = nfloor
              endif
           enddo

        enddo
     enddo
  enddo
end subroutine initchem

subroutine chem(q)
  use global
  implicit none

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: press,oldpress,dpress,pfloor
  integer :: ix,iy,iz

  do iz=lz,mz
     do iy=ly,my
        do ix=lx,mx
           ! compute pressure
           call pressure(q(:,ix,iy,iz),press)
           oldpress = press

           if (chemistry) then
              ! compute number densities of H and He in cell
              nH  = chiH  * q(1,ix,iy,iz)
              nHe = chiHe * q(1,ix,iy,iz)
              
              ! do it
              if (noncie) then
                 call ncie(q(1,ix,iy,iz),press,q(7:nq,ix,iy,iz))
              else
                 call cie(q(1,ix,iy,iz),press,q(7:nq,ix,iy,iz))
              endif
           endif

           ! update conservative variables
           if (cooling .or. tempfloor .or. jeansfloor) then
              if (tempfloor) then
                 call temp2press(temp_fl,q(1,ix,iy,iz),q(7:nq,ix,iy,iz),pfloor)
                 press = max(press,pfloor)
              endif

              if (jeansfloor) then
                 pfloor = jeansfloor0 * q(1,ix,iy,iz) * q(1,ix,iy,iz)
                 press = max(press,pfloor)
              endif

              dpress = (press - oldpress)

              q(5,ix,iy,iz) = q(5,ix,iy,iz) + dpress * g1inv
              q(6,ix,iy,iz) = q(6,ix,iy,iz) + dpress / (q(1,ix,iy,iz)**g1)
           endif
        enddo
     enddo
  enddo
end subroutine chem

subroutine compbackground()
  use global
  implicit none

  real(8) :: reddex

  if (background) then
     if (bg .eq. 1) then
        ! constant background
        uv = j0
     else if (bg .eq. 2) then
        ! simple model of black (1981)
        if (redshift .lt. 6 .and. redshift .gt. 3) then
           uv = j0 * 4.0 * cosmo
        else if (redshift .le. 3 .and. redshift .gt. 1) then
           uv = j0
        else
           uv = 0.0
        endif
     else if (bg .eq. 3) then
        ! model inspired by Bianchi et al (nonzero all the time)
        if (redshift .gt. 8.0) then
           uv = 1d-5 * j0
        else if (redshift .le. 8.0 .and. redshift .gt. 6.0) then
           reddex = (- redshift + 6.0) * 0.35
           uv = 0.5 * j0 * 10**reddex
        else if (redshift .le. 6.0 .and. redshift .gt. 3.0) then
           reddex = (- redshift + 3.0) * 0.1
           uv = j0 * 10**reddex
        else if (redshift .le. 3.0 .and. redshift .gt. 1.0) then
           uv = j0
        else if (redshift .le. 1.0) then
           uv = j0*0.1 * (10**redshift)
        endif
     endif
  endif

end subroutine compbackground

subroutine ncheck(n)
  use global
  implicit none

  real(8), dimension(ns) :: n

  integer :: i

  do i=1,ns1
     if ((n(i) .ge. 0.0 .and. n(i) .le. 1d20)  .eqv. .false.) then
        print *,"Ok, we're producing wrong abundances's!",i
        print '(6E10.2)',n
        call cleanstop
     endif
  enddo

end subroutine ncheck
