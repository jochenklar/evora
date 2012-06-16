!########################################################################!
!#                                                                      #!
!# muscl.f90                                                            #!
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

subroutine primitive(w,q)
  use global
  implicit none

  real(8), dimension(nq,lx:mx,ly:my,lz:mz):: q
  real(8), dimension(nw,lx2:mx2,ly2:my2,lz2:mz2):: w

  integer :: ix,iy,iz

  do iz=lz,mz
     do iy=ly,my
        do ix=lx,mx
           call q2w(q(:,ix,iy,iz),w(:,ix,iy,iz))
        enddo
     enddo
  enddo

end subroutine primitive

subroutine compevo(w,evo,s)
  use global
  implicit none

  real(8), dimension(nw,lx2:mx2,ly2:my2,lz2:mz2)   :: w
  real(8), dimension(nw,lx1:mx1,ly1:my1,lz1:mz1)   :: evo
  real(8), dimension(nw,3,lx1:mx1,ly1:my1,lz1:mz1) :: s
  real(8), dimension(nw,3) :: a

  real(8) :: dt2

  integer :: ix,iy,iz,i

  dt2 = - 0.5 * dt

  if (onedimension) then ! speedy version for 1D
     do ix=lx1,mx1
        call comp_slope(w(:,ix-1,1,1), &
             w(:,ix,1,1), w(:,ix+1,1,1),s(:,1,ix,1,1))
        a(:,1) = w(2,ix,1,1) * s(:,1,ix,1,1)
        a(5,1) = a(5,1) + g * w(5,ix,1,1) * s(2,1,ix,1,1)
        a(2,1) = a(2,1) + s(5,1,ix,1,1) / w(1,ix,1,1)

        evo(:,ix,1,1) = dt2 * a(:,1)*dxinv
     enddo
  else if (twodimension) then ! speedy version for 2D
     do iy=ly1,my1
        do ix=lx1,mx1
           call comp_slope(w(:,ix-1,iy,1), &
                w(:,ix,iy,1), w(:,ix+1,iy,1),s(:,1,ix,iy,1))
           call comp_slope(w(:,ix,iy-1,1), &
                w(:,ix,iy,1), w(:,ix,iy+1,1),s(:,2,ix,iy,1))

           a(:,1) = w(2,ix,iy,1) * s(:,1,ix,iy,1)
           a(5,1) = a(5,1) + g * w(5,ix,iy,1) * s(2,1,ix,iy,1)
           a(2,1) = a(2,1) + s(5,1,ix,iy,1) / w(1,ix,iy,1)

           a(:,2) = w(3,ix,iy,1) * s(:,2,ix,iy,1)
           a(5,2) = a(5,2) + g * w(5,ix,iy,1) * s(3,2,ix,iy,1)
           a(3,2) = a(3,2) + s(5,2,ix,iy,1) / w(1,ix,iy,1)

           evo(:,ix,iy,1) = dt2 * (a(:,1)*dxinv + a(:,2)*dyinv)
        enddo
     enddo
  else ! regular version for 3D
     do iz=lz1,mz1
        do iy=ly1,my1
           do ix=lx1,mx1
              call comp_slope(w(:,ix-1,iy,iz), &
                   w(:,ix,iy,iz), w(:,ix+1,iy,iz),s(:,1,ix,iy,iz))
              call comp_slope(w(:,ix,iy-1,iz), &
                   w(:,ix,iy,iz), w(:,ix,iy+1,iz),s(:,2,ix,iy,iz))
              call comp_slope(w(:,ix,iy,iz-1), &
                   w(:,ix,iy,iz), w(:,ix,iy,iz+1),s(:,3,ix,iy,iz))

              do i=1,3
                  a(:,i) = w(i+1,ix,iy,iz) * s(:,i,ix,iy,iz)
                  a(5,i) = a(5,i) + g * w(5,ix,iy,iz) * s(i+1,i,ix,iy,iz)
                  a(i+1,i) = a(i+1,i) + s(5,i,ix,iy,iz) / w(1,ix,iy,iz)
              enddo

              evo(:,ix,iy,iz) = dt2 * (a(:,1)*dxinv &
                   + a(:,2)*dyinv + a(:,3)*dzinv)
           enddo
        enddo
     enddo
  endif

end subroutine compevo

subroutine comp_slope(wl,w,wr,s)
  use global
  implicit none

  real(8), dimension(nw) :: wl,w,wr,s

  real(8), dimension(nw) :: dwl,dwr,dwc,pro
  integer :: i
  logical :: dwlpos,dwrpos

  dwl = w - wl
  dwr = wr - w

  if (limiter .eq. 1) then ! minmod
     do i=1,nw
        dwlpos = dwl(i) .ge. 0.0
        dwrpos = dwr(i) .ge. 0.0

        if (dwlpos .neqv. dwrpos) then
           s(i) = 0.0
        else if (dwrpos) then
           s(i) = min(dwl(i),dwr(i))
        else
           s(i) = max(dwl(i),dwr(i))
        endif
     enddo
  else if (limiter .eq. 2) then ! moncen
     dwc = (wr - wl) * 0.5
     do i=1,nw
        dwlpos = dwl(i) .ge. 0.0
        dwrpos = dwr(i) .ge. 0.0

        if (dwlpos .neqv. dwrpos) then
           s(i) = 0.0
        else if (dwrpos) then
           s(i) = min(2.0*dwl(i),dwc(i),2.0*dwr(i))
        else
           s(i) = max(2.0*dwl(i),dwc(i),2.0*dwr(i))
        endif
     enddo
  else if (limiter .eq. 3) then ! vanleer
     pro = dwl * dwr
     do i=1,nw
        if (pro(i) .le. 0.0) then
           s(i) = 0.0
        else
           s(i) = 2.0 * pro(i) / (dwr(i) + dwl(i))
        endif
     enddo
  else
     print *,'ERROR: wrong slope limiter!'
     stop
  endif

end subroutine comp_slope

subroutine sweep(w,q,evo,s)
  use global
  implicit none

  real(8), dimension(nq,lx:mx,ly:my,lz:mz)         :: q
  real(8), dimension(nw,lx2:mx2,ly2:my2,lz2:mz2)   :: w
  real(8), dimension(nw,lx1:mx1,ly1:my1,lz1:mz1)   :: evo
  real(8), dimension(nw,3,lx1:mx1,ly1:my1,lz1:mz1) :: s

  ! left and right state as seen by the interface
  real(8), dimension(nw):: wl,wr

  ! conservative fluxes as seen by the cell center
  real(8), dimension(nf):: fl,fr

  real(8) :: delta

  integer :: ix,iy,iz,dir

  ! xsweep
  if (nx .ne. 1) then
     dir = 1
     delta = dt * dxinv
     do iz=lz,mz
        do iy=ly,my
           wl = w(:,lx1,iy,iz)
           wr = w(:,lx,iy,iz)
           if (secondorder) then
              wl = wl + 0.5 * s(:,1,lx1,iy,iz) + evo(:,lx1,iy,iz)
              wr = wr - 0.5 * s(:,1,lx ,iy,iz) + evo(:,lx ,iy,iz)
           endif
           call riemannsolver(wl,wr,fr,dir)
           do ix=lx,mx
              fl = fr
              wl = w(:,ix,iy,iz)
              wr = w(:,ix+1,iy,iz)
              if (secondorder)  then
                 wl = wl + 0.5 * s(:,1,ix  ,iy,iz) + evo(:,ix  ,iy,iz)
                 wr = wr - 0.5 * s(:,1,ix+1,iy,iz) + evo(:,ix+1,iy,iz)
              endif
              call riemannsolver(wl,wr,fr,dir)
              call integration(fl,fr,q(:,ix,iy,iz),delta)
           enddo
        enddo
     enddo
  endif

  ! ysweep
  if (ny .ne. 1) then
     dir = 2
     delta = dt * dyinv
     do iz=lz,mz
        do ix=lx,mx
           wl = w(:,ix,ly1,iz)
           wr = w(:,ix,ly,iz)
           if (secondorder) then
              wl = wl + 0.5 * s(:,2,ix,ly1,iz) + evo(:,ix,ly1,iz)
              wr = wr - 0.5 * s(:,2,ix,ly ,iz) + evo(:,ix,ly ,iz)
           endif
           call riemannsolver(wl,wr,fr,dir)
           do iy=ly,my
              fl = fr
              wl = w(:,ix,iy,iz)
              wr = w(:,ix,iy+1,iz)
              if (secondorder) then
                 wl = wl + 0.5 * s(:,2,ix,iy  ,iz) + evo(:,ix,iy  ,iz)
                 wr = wr - 0.5 * s(:,2,ix,iy+1,iz) + evo(:,ix,iy+1,iz)
              endif
              call riemannsolver(wl,wr,fr,dir)
              call integration(fl,fr,q(:,ix,iy,iz),delta)
           enddo
        enddo
     enddo
  endif

  ! zsweep
  if (nz .ne. 1) then
     dir = 3
     delta = dt * dzinv
     do iy=ly,my
        do ix=lx,mx
           wl = w(:,ix,iy,lz1)
           wr = w(:,ix,iy,lz)
           if (secondorder) then
              wl = wl + 0.5 * s(:,3,ix,iy,lz1) + evo(:,ix,iy,lz1)
              wr = wr - 0.5 * s(:,3,ix,iy,lz ) + evo(:,ix,iy,lz )
           endif
           call riemannsolver(wl,wr,fr,dir)
           do iz=lz,mz
              fl = fr
              wl = w(:,ix,iy,iz)
              wr = w(:,ix,iy,iz+1)
              if (secondorder) then
                 wl = wl + 0.5 * s(:,3,ix,iy,iz  ) + evo(:,ix,iy,iz  )
                 wr = wr - 0.5 * s(:,3,ix,iy,iz+1) + evo(:,ix,iy,iz+1)
              endif
              call riemannsolver(wl,wr,fr,dir)
              call integration(fl,fr,q(:,ix,iy,iz),delta)
           enddo
        enddo
     enddo
  endif

end subroutine sweep

subroutine integration(fl,fr,q,delta)
  use global
  implicit none

  real(8), intent(in), dimension(nf):: fl,fr
  real(8), intent(inout), dimension(nq):: q
  real(8), intent(in) :: delta

  if (noncie) then
     q(1:nf) = q(1:nf) + delta * (fl(1:nf) - fr(1:nf))
  else
     q(1:6) = q(1:6) + delta * (fl(1:6) - fr(1:6))
  endif

end subroutine integration
