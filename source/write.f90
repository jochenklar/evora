!########################################################################!
!#                                                                      #!
!# write.f90                                                            #!
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

subroutine writedump(q,name)
  use global
  implicit none

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q
  character(len=*) :: name

 
  character(len=127) :: filename
  character(len=4) :: node

  ! dump consevative quantities
  write (node,'(".",I0.3)') mpirank

  filename = name//node
  call writedata(q,nq,filename)

  if (master) call writeunits(name)

end subroutine writedump

subroutine writedata(a,nn,filename)
  use global
  implicit none

  real(8), dimension(nn,lx:mx,ly:my,lz:mz) :: a
  integer :: nn,ierror
  character(len=*) :: filename

  open (unit=10, file="data/"//filename,form="unformatted",iostat=ierror)
  if (ierror .ne. 0) then
     print *,"ERROR file not found: ",filename
     stop
  endif

  write (10) nn,nx,ny,nz,lz,mz
  write (10) a
  close (10)

end subroutine writedata

subroutine writebovhead(name,var,dir)
  use global
  implicit none

  character(len=*) :: name,var,dir

  character(len=127) :: filename1,filename2,varname
  integer :: ierror
  
  filename1 = trim(dir)//trim(name)//"."//trim(var)//".bov"
  filename2 = trim(name)//"."//trim(var)//".bovdata"

  if (trim(var) .eq. "dens") then
     varname = "Density"
  else if (trim(var) .eq. "temp") then
     varname = "Temperature"
  endif

  open (unit=10,file=trim(filename1), access="sequential", &
       form="formatted",iostat=ierror)
  call checkfile(ierror,filename1)

  if (supercomoving) then
     write (10,'("TIME: ",E10.2)') redshift
  else
     write (10,'("TIME: ",E10.2)') t
  endif
  write (10,'("DATA_FILE: ",A20)') filename2
  write (10,'("DATA_SIZE: ",3I5)') nx,ny,nz
  write (10,'("DATA_FORMAT: FLOAT")')
  write (10,'("VARIABLE:",A20)') trim(varname)
  write (10,'("DATA_ENDIAN: LITTLE")')
  write (10,'("CENTERING: ZONAL")')
  write (10,'("BRICK_ORIGIN: ",3E12.2)') xmin,ymin,zmin
  write (10,'("BRICK_SIZE: ",3E12.2)') xmax-xmin,ymax-ymin,zmax-zmin
  close (10)

end subroutine writebovhead

subroutine writecut(a,n,nvar,infile,name)
  use global
  implicit none

  real(8), dimension(nvar,1:n) :: a
  integer :: n,nvar
  character(len=*) :: infile,name

  integer :: ierror
  integer :: i
  character(len=127) :: filename
  character(len=9) :: form

  filename = "out/"//trim(infile)//"."//name//".cut"
  open (unit=25, file=trim(filename), &
       access="sequential",status="unknown", form="formatted",iostat=ierror)
  call checkfile(ierror,filename)
  write (25,'(A,A15,26A16)') "#","x","rho","rho u","rho v","rho w","E","S", &
       "rho_phys","u","v","w","p","T","xi_HI","xi_HII","xi_HeI","xi_HeII", &
       "xi_HeIII","xi_e","n_HI","n_HII","n_HeI","n_HeII","n_HeIII","n_e", &
       "lambda","gamma"
  write (form,'(A,I2,A6)') "(",nvar,"E16.6)"
  do i=1,n
     write (25,form) a(:,i)
  enddo
  close(25)
end subroutine writecut

subroutine writeslice(a,n1,n2,nvar,infile,name)
  use global
  implicit none

  real(8), dimension(nvar,1:n1,1:n2) :: a
  integer :: n1,n2,nvar
  character(len=*) :: infile,name

  integer :: ierror
  integer :: i1,i2
  character(len=127) :: filename
  character(len=9) :: form

  filename = "out/"//trim(infile)//"."//name//".slice"
  open (unit=25, file=trim(filename), &
       access="sequential",status="unknown", form="formatted",iostat=ierror)
  call checkfile(ierror,filename)
  write (25,'(A,A15,27A16)') "#","x","y","rho","rho u","rho v","rho w","E","S", &
       "rho_phys","u","v","w","p","T","xi_HI","xi_HII","xi_HeI","xi_HeII", &
       "xi_HeIII","xi_e","n_HI","n_HII","n_HeI","n_HeII","n_HeIII","n_e", &
       "lambda","gamma"
  write (form,'(A,I2,A6)') "(",nvar,"E16.6)"
  do i1=1,n1
     do i2=1,n2
        write (25,form) a(:,i1,i2)
     enddo
     write (25,*) " "
  enddo
  close(25)

end subroutine writeslice

subroutine checkfile(ierror,filename)
  implicit none
  
  integer :: ierror
  character(len=*) :: filename 

  if (ierror .ne. 0) then
     print *,"ERROR problem while opening :",trim(filename)
     stop
  endif

end subroutine checkfile

subroutine opendata(i,infile)
  use global
  implicit none

  integer, intent(in)      :: i
  character(*), intent(in) :: infile

  character(len=127) :: filename
  character(len=4) :: node

  integer :: nn,fnx,fny,fnz,ierror

  write(node,'(".",I0.3)') i

  filename = "data/"//trim(infile)//node

  print *,"processing file ",trim(filename)

  ! open file with conservative quantities
  open (unit=10, file=trim(filename),form="unformatted", &
       status="old",iostat=ierror)
  call checkfile(ierror,filename)
  read (10) nn,fnx,fny,fnz,lz,mz

  ! check
  if (nn .ne. nq) then
     print *,"ERROR not a proper data file ",filename,nn,"!=",nq
     stop
  endif
  if (fnx .ne. nx .or. fny .ne. ny .or. fnz .ne. nz) then
     print *,"ERROR dimensions do not match ",trim(filename)
     print *,fnx,fny,fnz,"!=",nx,ny,nz
     stop
  endif

end subroutine opendata

subroutine getstring(i,a)
  use global
  implicit none

  integer, intent(in)  :: i
  character(*), intent(out) :: a
  
  call get_command_argument(i,a)
  if (len_trim(a) .eq. 0) then
     print *,"ERROR give string as agument!"
     stop
  endif
end subroutine getstring
