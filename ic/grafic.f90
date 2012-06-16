!########################################################################!
!#                                                                      #!
!# grafic.f90                                                           #!
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

  real(8), parameter :: tol = 1d-6

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8), parameter :: dmean=1.0, tinit=100.0 

  real, dimension(:,:,:), allocatable :: buffer
  integer :: np1,np2,np3
  real :: g_dx,tmp,g_cosmo,g_omega0,g_lambda0,g_hubble0

  real(8) :: x,y,z
  integer :: ix,iy,iz
  
  real(8) :: p,e,kin,temp,vscale

  integer :: i,ierror

  ! some checks first 
  if (master) then
     open(11,file="input/ic/ic_deltab",form='unformatted', iostat=ierror)
     if (ierror .ne. 0) then
        print *,"ERROR file not found input/ic_deltab"
        call cleanstop
     endif
     read(11) np1,np2,np3,g_dx,tmp,tmp,tmp,g_cosmo,g_omega0,g_lambda0,g_hubble0
     close(11)
     if (np1 .ne. nx .or. np2 .ne. ny .or. np3 .ne. nz) then
        print *,"ERROR dimensions do not match in grafics file"
        print *,np1,np2,np3
        print *,nx,ny,nz
        call cleanstop
     endif
     if (abs(nx*g_dx*h - box) .gt. tol) then
        print *,"ERROR boxsize does not match in grafics file"
        print *,nx*g_dx*h,box
        call cleanstop
     endif
     if (abs(g_omega0-omega0) .gt. tol .or. abs(g_lambda0-lambda0) .gt. tol &
          .or. abs(g_hubble0-hubble0) .gt. tol)  then
        print *,"WARNING cosmology does not match in grafics file"
        print *,dble(g_omega0),omega0
        print *,dble(g_lambda0),lambda0
        print *,dble(g_hubble0),hubble0
        call cleanstop
     endif
     zstart = 1.0 / g_cosmo - 1.0
     call cosmostart()
  endif

  call bcasttime()
  if (master) call compunits()
  call bcastunits()

  allocate(buffer(lx:mx,ly:my,lz:mz))

  vscale = cosmo / (box * 100.0)

  do i=1,mpisize
     if (mpirank .eq. (i-1)) then

        open(11,file="data/ic/ic_deltab",form='unformatted')
        open(12,file="data/ic/ic_velbx",form='unformatted')
        open(13,file="data/ic/ic_velby",form='unformatted')
        open(14,file="data/ic/ic_velbz",form='unformatted')

        read(11) ! skip first line (parameters)
        read(12) ! skip first line (parameters)
        read(13) ! skip first line (parameters)
        read(14) ! skip first line (parameters)

        do iz=1,lz-1
           read(11) ! skip some lines
           read(12) ! skip some lines
           read(13) ! skip some lines
           read(14) ! skip some lines
        end do

        ! density
        do iz=lz,mz
           read(11) ((buffer(ix,iy,iz),ix=1,nx),iy=1,ny)
        end do
        q(1,lx:mx,ly:my,lz:mz) = dble(buffer(lx:mx,ly:my,lz:mz)) + 1.0

        ! velocity x
        do iz=lz,mz
           read(12) ((buffer(ix,iy,iz),ix=1,nx),iy=1,ny)
        end do
        q(2,lx:mx,ly:my,lz:mz) = &
             dble(buffer(lx:mx,ly:my,lz:mz)) * q(1,lx:mx,ly:my,lz:mz) * vscale

        ! velocity y
        do iz=lz,mz
           read(13) ((buffer(ix,iy,iz),ix=1,nx),iy=1,ny)
        end do
        q(3,lx:mx,ly:my,lz:mz) = &
             dble(buffer(lx:mx,ly:my,lz:mz)) * q(1,lx:mx,ly:my,lz:mz) * vscale

       ! velocity z
        do iz=lz,mz
           read(14) ((buffer(ix,iy,iz),ix=1,nx),iy=1,ny)
        end do
        q(4,lx:mx,ly:my,lz:mz) = &
             dble(buffer(lx:mx,ly:my,lz:mz)) * q(1,lx:mx,ly:my,lz:mz) * vscale

        close(11)
        close(12)
        close(13)
        close(14)
     endif
     call barrier
  enddo

  ! compute energy
  do iz=lz,mz
     do iy=ly,my       
        do ix=lx,mx
           temp = tinit * q(1,ix,iy,iz)**g1
           p = q(1,ix,iy,iz) * temp  / temp0
           kin = ((q(2,ix,iy,iz)**2) + (q(3,ix,iy,iz)**2) &
           + (q(4,ix,iy,iz)**2)) * 0.5 / q(1,ix,iy,iz)
           e = kin + p / g1

           q(5,ix,iy,iz) = e
           q(6,ix,iy,iz) = p / (q(1,ix,iy,iz)**g1)
        enddo
     enddo
  enddo

end subroutine ic
