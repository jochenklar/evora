!########################################################################!
!#                                                                      #!
!# startstop.f90                                                        #!
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

subroutine startevora
  use mpi
  use global
  implicit none

  call mpi_init(mpierror)
  call mpi_comm_size(mpi_comm_world,mpisize,mpierror)
  call mpi_comm_rank(mpi_comm_world,mpirank,mpierror)

  if (mpirank .eq. 0) then
     master = .true.
     notmaster = .false.
  else
     master = .false.
     notmaster = .true.
  endif

  if (mpirank .eq. mpisize-1) then
     last = .true.
     notlast = .false.
  else
     last = .false.
     notlast = .true.     
  endif

  if (mod(mpirank,2) .eq. 0) then
     even = .true.
     noteven = .false.
  else
     even = .false.
     noteven = .true.     
  endif

  if (mpisize .eq. 1) then
     parallel = .false.
  else
     parallel = .true.
  endif

end subroutine startevora

subroutine stopevora()
  use mpi
  use global
  implicit none

  call mpi_finalize(mpierror)
  stop

end subroutine stopevora 

subroutine cleanstop()
  use mpi
  use global
  implicit none

  write (*,'("abort: thread ",I0.3)') mpirank

  call mpi_abort(mpi_comm_world,666,mpierror)

end subroutine cleanstop
