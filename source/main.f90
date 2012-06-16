!########################################################################!
!#                                                                      #!
!# main.f90                                                             #!
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

program evora
  use global
  implicit none

  ! declare some variables
  logical :: go = .true.
  logical :: dump = .false.
  integer :: k=0,k1=0,k2=0,k3=0,k4=1
  character(127) :: name

  ! declare the main arrays
  real(8), dimension(:,:,:,:), allocatable :: q,w,evo
  real(8), dimension(:,:,:,:,:), allocatable :: s
  complex(8), dimension(:,:,:), allocatable :: phi
  real(8), dimension(:,:,:,:), allocatable :: buf1,buf2,buf3,buf4
  real(8), dimension(:,:), allocatable :: phibuf1,phibuf2,phibuf3,phibuf4

  ! initializes mpi if necessary
  call startevora

  ! read, broadcast and compute the parameters
  if (master) then
     call readparameter()
     call writeparameter()
  endif
  call bcastparameter()
  call compparameter()

  ! say hello
  if (master) call startupscreen() ! also performs some checks

  ! compute initial expansion factor and communicate
  if (supercomoving .and. master) call cosmostart()
  call bcasttime()

  ! compute units
  if (master) call compunits()
  call bcastunits()

  ! initialize gravitation
  if (gravity) call gravstart(phi)

  ! make domain decomposition
  call decomposition()

  ! allocate main arrays
  allocate(q(nq,lx:mx,ly:my,lz:mz))         ! conservative quantities
  allocate(w(nw,lx2:mx2,ly2:my2,lz2:mz2))   ! primitive quantities
  allocate(evo(nw,lx1:mx1,ly1:my1,lz1:mz1)) ! evolution term
  allocate(s(nw,3,lx1:mx1,ly1:my1,lz1:mz1)) ! slopes
  evo = 0.0
  s = 0.0

  ! build the initial conditions
  call ic(q)

  ! allocate array for gravitational potential
  if (gravity) allocate(phi(lx:mx,ly:my,lz:mz))

  ! allocate arrays buffer arrays for mpi
  if (parallel) allocate( &
          buf1(nw,lx2:mx2,ly2:my2,2), buf2(nw,lx2:mx2,ly2:my2,2), &
          buf3(nw,lx2:mx2,ly2:my2,2), buf4(nw,lx2:mx2,ly2:my2,2))
  if (gravity .and. parallel) allocate( &
          phibuf1(lx:mx,ly:my),phibuf2(lx:mx,ly:my), &
          phibuf3(lx:mx,ly:my),phibuf4(lx:mx,ly:my))

  ! initialize chemistry
  if (chemistry) call initchem(q)

  ! dump arrays into datafiles!
  call writedump(q,"in")

  ! first screen output
  call evolution(k2,q)

  ! first slices and cuts
  if (cut .or. slice .or. bov) call snapshot(k3,q)

  ! stop if iconly
  if (iconly) call stopevora

  ! here begins the timestep
  do while (go)

     ! compute primitive variables
     call primitive(w,q)

     ! compute boundary cells
     call bound(w)

     ! compute timestep
     dt = 1d99
     if (hydro)      call hydrostep(w)
     if (cooldt)     call coolstep(q,w)
     if (conduction) call thermalstep(q,w)
     if (parallel)   call reducestep        ! communicate the timestep

     if (master) then
        if (supercomoving) then
           call cosmostep()               ! compute timestep from expansion
           call compcosmo()              ! update cosmo
        endif

        call updatetime(go,dump,k4) ! update time and checks for end of simulation

        if (supercomoving) then
           call compunits()               ! update units
           if (background) call compbackground
           if (chemistry) call compcoef() ! update chemical coeficcients
        endif
     endif

     ! main communication step
     if (parallel) then
        call bcastlogical(go)
        call bcastlogical(dump)
        call bcasttime                          ! comm. of time and cosmo
        if (supercomoving) then
           call bcastunits()                    ! comm. of units
           if (chemistry) call bcastcoef()      ! comm. of chemical coef.
          endif
        call communicate(w,buf1,buf2,buf3,buf4) ! comm. of boundary
     endif

     ! compute the predictor step
     if (secondorder) call compevo(w,evo,s)

     ! this is the main hydro step
     if (hydro) call sweep(w,q,evo,s)

     if (gravity) then
        ! compute the gravitational potential via fft
        call gravpotential(q,phi)
        ! communicate potential boundary cells
        if (parallel) call gravcomm(phi,phibuf1,phibuf2,phibuf3,phibuf4)
        ! compute gravitational force and update conservative variables
        call gravintegration(q,phi,phibuf2,phibuf4)
     endif

     ! compute chemical evolution (and cooling)
     if (chemistry .or. tempfloor .or. jeansfloor) call chem(q)

     ! compute and apply thermal conduction
     if (conduction) then
        call primitive(w,q)
        call bound(w)
        if (parallel) call communicate(w,buf1,buf2,buf3,buf4)
        call thermalconduction(q,w)
     endif

     ! compute and append expansion terms in energy and entropy equation
     if (supercomoving .and. .not. atomic) call expansion(q)

     ! update timestep counter
     k = k + 1

     ! check for timestep
     if (supercomoving) then
        tsnap = cosmo * nsnap
     else
        tsnap = (t - tstart)/(tend - tstart) * nsnap
     endif

     if (master .and. tsnap .ge. dble(k1)) call screen(k1,k)
     if (tsnap .ge. dble(k2)) call evolution(k2,q)

     if (cut .or. slice .or. bov) then
        if (tsnap .ge. dble(k3)) call snapshot(k3,q)
     endif

     if (dump) then
        write (name,'(I0.3)') k4
        call writedump(q,trim(name))
        dump = .false.
        k4 = k4 + 1
     endif
  enddo

  ! final dump
  call writedump(q,"out")

  ! destroy fftw plans
  if (gravity) call gravstop()

  ! free memory
  deallocate(q,w,evo)
  if (gravity) deallocate(phi)
  if (parallel) deallocate(buf1,buf2,buf3,buf4)
  if (gravity .and. parallel) &
      deallocate(phibuf1,phibuf2,phibuf3,phibuf4)

  if (master) call endscreen

  call stopevora ! finalises mpi if necessary

end program evora
