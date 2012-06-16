!########################################################################!
!#                                                                      #!
!# fft.f90                                                         #!
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

subroutine gravstart()
  use mpi
  use global
  implicit none

  ! creating fftw plans
  call fftw3d_f77_mpi_create_plan(there,mpi_comm_world,nx,ny,nz, &
       FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE);
  call fftw3d_f77_mpi_create_plan(back,mpi_comm_world,nx,ny,nz, &
       FFTW_BACKWARD,FFTW_ESTIMATE+FFTW_IN_PLACE);
  call fftwnd_f77_mpi_local_sizes(there,local_nx,local_x_start, &
       local_ny_after_transpose,local_y_start_after_transpose,total_local_size)

end subroutine gravstart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gravstop()
  use mpi
  use global
  implicit none

  ! destroying fftw plans
  call fftwnd_f77_mpi_destroy_plan(there)
  call fftwnd_f77_mpi_destroy_plan(back)

end subroutine gravstop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fft_c2c_forward(a)
  use mpi
  use global
  implicit none

  complex(8), dimension(lx:mx,ly:my,lz:mz) :: a

  call fftwnd_f77_mpi(there,1,a,0,0,FFTW_NORMAL_ORDER)

end subroutine fft_c2c_forward

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fft_c2c_backward(a)
  use mpi
  use global
  implicit none

  complex(8), dimension(lx:mx,ly:my,lz:mz) :: a
  
  call fftwnd_f77_mpi(back,1,a,0,0,FFTW_NORMAL_ORDER)

end subroutine fft_c2c_backward

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fft_c2c_power(a)
  use global
  implicit none

  complex(8), dimension(1:nx,1:ny,1:nz) :: a

  call fftw3d_f77_create_plan(there,nx,ny,nz, &
       FFTW_FORWARD,FFTW_ESTIMATE+FFTW_IN_PLACE);
  call fftwnd_f77_one(there,a,0)
  call fftwnd_f77_destroy_plan(there)

end subroutine fft_c2c_power
