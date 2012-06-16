!########################################################################!
!#                                                                      #!
!# units.f90                                                            #!
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

subroutine compunits()
  use global
  implicit none

  real(8) :: h,d0,l0,t0,m0,kappa_coef,mfp_coef

  if (supercomoving) then
    h  = hubble0 * 0.01
    d0 = omega0 * rhoc * h**2.0
    t0 = 1.0 / (hubble0 * kilometer / megapc)
    l0 = box * megapc / h
    m0 = atomicmass
  else
    h  = 1.0
    d0 = 1.0
    t0 = 1.0
    l0 = 1.0
    m0 = kb
  endif

  density0      = d0 / (cosmo**3.0)
  length0       = l0 * cosmo
  time0         = t0 * (cosmo**2.0)
  velocity0     = l0 / cosmo / t0
  temp0         = (l0**2.0) / (cosmo**2.0) / (t0**2.0) * m0 / kb
  temp0inv      = 1.0 / temp0
  press0        = d0 * l0*l0 / (t0*t0) / (cosmo**5.0)

  number0       = d0 / (cosmo**3.0) / m0

  chem0         = d0 * t0 / cosmo / m0
  photo0        = t0 * cosmo**2.0

  cool0         = g1 * (t0**3.0) * cosmo * d0 / (l0**2.0) / (m0**2.0) 
  heat0         = g1 * (t0**3.0) *(cosmo**4.0) / (l0**2.0) / m0

  kappa_coef    = 4.6d13 * 1d-20 / (37.8 / 40.0)
  kappa0        = (cosmo**5.0) * (t0**3.0) / d0 / (l0**4.0) * 0.5 * kappa_coef
  mfp_coef      = 0.023 * megapc * 1d-16 * 1d-3
  meanfreepath0 = (cosmo**2.0) * m0 / l0 / d0 * 0.5 * 4.2 * mfp_coef

  jeans0        = sqrt(pi * g * fb)
  if (jeansfloor) then
     jeansfloor0   = (dx / jeans0 / c_jeans)**2
  else
     jeansfloor0   = 0.0
  endif

end subroutine compunits

subroutine writeunits(name)
  use global
  implicit none

  integer :: ierror
  character(len=*) :: name
  character(len=24) :: filename

  ! write binary
  filename = "data/"//trim(name)//".units"
  open (unit=10, file=trim(filename),form="formatted",iostat=ierror)
  if (ierror .ne. 0) then
     print *,"ERROR file not found: ",filename
     stop
  endif

  write (10,'("t =             ",E16.8)') t
  write (10,'("cosmo =         ",E16.8)') cosmo
  write (10,'("redshift =      ",E16.8)') redshift
  write (10,'("length0 =       ",E16.8)') length0
  write (10,'("time0 =         ",E16.8)') time0
  write (10,'("density0 =      ",E16.8)') density0
  write (10,'("velocity0 =     ",E16.8)') velocity0
  write (10,'("temp0 =         ",E16.8)') temp0
  write (10,'("press0 =        ",E16.8)') press0
  write (10,'("number0 =       ",E16.8)') number0
  write (10,'("chem0 =         ",E16.8)') chem0
  write (10,'("photo0 =        ",E16.8)') photo0
  write (10,'("cool0 =         ",E16.8)') cool0
  write (10,'("heat0 =         ",E16.8)') heat0
  write (10,'("jeans0 =        ",E16.8)') jeans0
  write (10,'("jeansfloor0 =   ",E16.8)') jeansfloor0
  write (10,'("kappa0 =        ",E16.8)') kappa0
  write (10,'("meanfreepath0 = ",E16.8)') meanfreepath0
  
  close (10)

end subroutine writeunits

subroutine readunits(name)
  use global
  implicit none

  character(len=*) :: name

  character(len=127) :: filename

  integer :: ierror
  character(16) :: tmp

  filename = "data/"//trim(name)//".units"
  open (unit=10,file=trim(filename),form="formatted", &
       iostat=ierror,status="old")
  if (ierror .ne. 0) then
     print *,"ERROR file not found: ",trim(filename)
     stop
  endif

  read (10,'(A16,E16.8)') tmp,t
  read (10,'(A16,E16.8)') tmp,cosmo
  read (10,'(A16,E16.8)') tmp,redshift
  read (10,'(A16,E16.8)') tmp,length0
  read (10,'(A16,E16.8)') tmp,time0
  read (10,'(A16,E16.8)') tmp,density0
  read (10,'(A16,E16.8)') tmp,velocity0
  read (10,'(A16,E16.8)') tmp,temp0
  read (10,'(A16,E16.8)') tmp,press0
  read (10,'(A16,E16.8)') tmp,number0
  read (10,'(A16,E16.8)') tmp,chem0
  read (10,'(A16,E16.8)') tmp,photo0
  read (10,'(A16,E16.8)') tmp,cool0
  read (10,'(A16,E16.8)') tmp,heat0
  read (10,'(A16,E16.8)') tmp,jeans0
  read (10,'(A16,E16.8)') tmp,jeansfloor0
  read (10,'(A16,E16.8)') tmp,kappa0
  read (10,'(A16,E16.8)') tmp,meanfreepath0
  
  close (10)

end subroutine readunits
