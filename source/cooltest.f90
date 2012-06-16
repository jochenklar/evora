!########################################################################!
!#                                                                      #!
!# cooltest.f90                                                         #!
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

program cooltest
  use global
  implicit none

  real(8), parameter :: dt_chem = 0.01
  real(8), parameter :: dt_noncietime = 0.01
  real(8), parameter :: dt_cool = 0.01

  logical, parameter :: ciechem = .true. 
  !logical, parameter :: ciechem = .false.
  logical, parameter :: nonciechem = .true.
  !logical, parameter :: nonciechem = .false.
  logical, parameter :: noncietime = .true.
  !logical, parameter :: noncietime = .false.
  logical, parameter :: ciecool = .true. 
  !logical, parameter :: ciecool = .false. 
  logical, parameter :: nonciecool = .true.
  !logical, parameter :: nonciecool = .false.
  logical, parameter :: rates = .true.
  !logical, parameter :: rates = .false.

  character(80) :: rhostart

  real(8), parameter :: null=0.0 
  integer, parameter :: grid=10000
  real(8), parameter :: tempstart=100.0
  real(8), parameter :: temptime=1d6

  real(8), dimension(grid) :: temp
  real(8) :: expo,rho,press,oldpress

  real(8), dimension(ns) :: n,nstart
  real(8) :: temp1,lambda

  real(8), dimension(nchem) :: chemrate
  real(8), dimension(ncool) :: coolrate
  real(8), dimension(ncool) :: coolrateN

  real(8) :: hh, cc

  integer :: i

  call get_command_argument(2,rhostart)
  if (len_trim(rhostart) .eq. 0) then
     print *,"ERROR give command line argument!"
     print *,"USAGE bin/cooltest <parameter file> <density>"
     stop
  endif

  read (rhostart,*) rho

  call readparameter()
  call compparameter()
  redshift = zend
  cosmo = 1.0 / (redshift + 1.0)
  uv = 0
  call compunits()
  call compcoef()

  ! build temperature and pressure scale
  do i=1,grid
     expo = 2.0 + 6.0 * dble(i-1) / grid
     temp(i) = 10.0**expo
  enddo

  ! nuclear abundances  
  chiH = rhoH / mH
  chiHe = rhoHe / mHe

  nH  = chiH  * rho
  nHe = chiHe * rho

  substep = 0.01

  print *,"------------------------------------"
  print *,"evora - cooltest"
  print *,"------------------------------------"
  print *,"nH  =     ",nH*number0
  print *,"nHe =     ",nHe*number0
  print *,"rho =     ",rho
  print *,"n0  =     ",number0
  print *,"------------------------------------"
  print *,"dt chem = ",dt_chem
  print *,"dt cool = ",dt_noncietime
  print *,"dt cool = ",dt_cool
  print *,"------------------------------------"

  ! start values for q and n
  cooling = .false.
  call comprates(tempstart,chemrate,coolrate)
  call cieabundance(nstart,chemrate)

  if (background) then
     uv = j0
     call compbackground()
     call compcoef()
     print *,"uv = ",uv
     print *,"------------------------------------"
  endif

  ! --- testing chemical evolution into equilibrium --------------------

  dt = dt_chem ! in hubble units

  if (ciechem) then
     print *,"CIE..."
     open(unit=10,file="out/cie.chem")
     do i=1,grid
        n = nstart
        press = temp(i) * rho / temp0
        call cie(rho,press,n)  
        write (10,'(7E14.6)') temp(i),n*number0
     enddo
     close(10)
  endif

  if (nonciechem) then
     print *,"non-CIE..."
     open(unit=10,file="out/noncie.chem")
     do i=1,grid
        n = nstart
        press = temp(i) * rho / temp0
        call ncie(rho,press,n)     
        write (10,'(7E14.6)') temp(i),n*number0
     enddo
     close(10)
  endif

  ! --- non-equilibrium evolution with time ----------------------------

  if (noncietime) then
     print *,"non-CIE evolution with time..."
     open(unit=10,file="out/noncietime.chem")
     temp1 = temptime
     t = null
     n = nstart

     ! call temp2press(temp1,rho,n,press1)

     do i=1,grid     
        dt = dt_noncietime / dble(grid)
        t = t + dt
        call temp2press(temp1,rho,n,press)
        call ncie(rho,press,n)
        write (10,'(7E14.6)') t*time0,n*number0
     enddo
     close(10)
  endif

  ! --- testing cooling function ---------------------------------------

  dt = dt_cool
  cooling = .true.

  if (ciecool) then
     print *,"CIE with cooling..."
     open(unit=10,file="out/cie.cool")
     do i=1,grid
        call comprates(temp(i),chemrate,coolrate)
        call cieabundance(n,chemrate)
        call temp2press(temp(i),rho,n,press)
        oldpress = press
        call cie(rho,press,n)
        cc = max(oldpress-press,ratefloor)
        hh = max(press-oldpress,ratefloor)
        write (10,'(16E20.10)') temp(i),press,hh,cc,n*number0,nstart*number0
     enddo
     close(10)
  endif

  if (nonciecool) then
     print *,"non-CIE with cooling..."
     open(unit=10,file="out/noncie.cool")
     do i=1,grid
        call comprates(temp(i),chemrate,coolrate)
        call cieabundance(n,chemrate)
        call temp2press(temp(i),rho,n,press)
        oldpress = press
        call ncie(rho,press,n)
        cc = max(oldpress-press,ratefloor)
        hh = max(press-oldpress,ratefloor)
        write (10,'(16E20.10)') temp(i),press,hh,cc,n*number0,nstart*number0
     enddo
     close(10)
  endif

  ! --- plot the cooling function and the rates ------------------------

  if (rates) then
     print *,"rates and cooling function..."
     open(unit=50,file="out/cool.cool")
     open(unit=60,file="out/rate.chem")
     open(unit=70,file="out/rate.cool")
     do i=1,grid  
        call comprates(temp(i),chemrate,coolrate)
        call cieabundance(n,chemrate)
        call temp2press(temp(i),rho,n,press)

        write (60,'(10E14.6)') temp(i),chemrate
        write (70,'(14E14.6)') temp(i),coolrate

        if (cooling) then
           coolrateN(1)  = coolrate(1) * n(1)
           coolrateN(2)  = coolrate(2) * n(3)
           coolrateN(3)  = coolrate(3) * n(4)
           coolrateN(4)  = coolrate(4) * n(2)
           coolrateN(5)  = coolrate(5) * n(4)
           coolrateN(6)  = coolrate(6) * n(5)
           coolrateN(7)  = coolrate(7) * n(4)
           coolrateN(8)  = coolrate(8) * n(1)
           coolrateN(9)  = coolrate(9) * n(4)
           coolrateN(10) = coolrate(10) * (n(2) + n(4) + 4.0 * n(5))
           
           coolrateN(1:10) = coolrateN(1:10) * n(6)
        else
           coolrateN(1:10) = 0.0
        endif

        if (background .and. cooling) then
           coolrateN(11) = n(1) * coolrate(11)
           coolrateN(12) = n(3) * coolrate(12)
           coolrateN(13) = n(4) * coolrate(13)
        else
           coolrateN(11:13) = 0.0
        endif

        coolrateN = coolrateN/cool0/(nH*nH)
        coolrateN = max(coolrateN,ratefloor)

        ! compute cooling function
        lambda = sum(coolrateN(11:13)) - sum(coolrateN(1:10))
        write (50,'(17E14.6)') temp(i),coolrateN,abs(lambda), &
             abs(lambda)/(number0*number0),press
     enddo
     close(50)
     close(60)
     close(70)
  endif

  print *,"done!"
  print *,"------------------------------------"

end program cooltest
