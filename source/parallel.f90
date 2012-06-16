!########################################################################!
!#                                                                      #!
!# parallel.f90                                                         #!
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

subroutine decomposition()
  ! domain decomposition (given by fftw or slab decomposition)
  use mpi
  use global
  implicit none

  integer :: i

  lx  = 1      ; mx  = nx
  lx1 = lx - 1 ; lx2 = lx - 2 ; lx0 = lx + 1
  mx1 = mx + 1 ; mx2 = mx + 2 ; mx0 = mx - 1

  ly  = 1      ; my  = ny
  ly1 = ly - 1 ; ly2 = ly - 2 ; ly0 = ly + 1
  my1 = my + 1 ; my2 = my + 2 ; my0 = my - 1

  if (nz .lt. mpisize) then
     if (master) then
        print *,"ERROR: to few cells in z-direction for decomposition"
        print *,"nz = ",nz,"mpisize = ",mpisize
     endif
     call cleanstop
  endif

  if (gravity) then ! decomposition by fftw
     lz = local_x_start + 1
     mz = local_x_start + local_nx
  else ! slab decomposition z-wise
     mz = 0
     do i=0,mpisize-1
        if (mpirank .ge. i) then
           lz = mz + 1
           mz = lz + (nz - lz) / (mpisize - i)
        endif
     enddo
  endif

  lz1 = lz - 1 ; lz2 = lz - 2 ; lz0 = lz + 1
  mz1 = mz + 1 ; mz2 = mz + 2 ; mz0 = mz - 1

end subroutine decomposition

subroutine bcasttime() 
  use mpi
  use global
  implicit none

  integer, parameter :: nbuf=6
  real(8), dimension(nbuf) :: buf

  if (master) then
     buf(1) = t
     buf(2) = dt
     buf(3) = cosmo
     buf(4) = expand
     buf(5) = redshift
     buf(6) = uv
  endif
  
  call mpi_bcast(buf,nbuf,mpi_double_precision,0,mpi_comm_world,mpierror)

  if (notmaster) then
     t        = buf(1)
     dt       = buf(2)
     cosmo    = buf(3)
     expand   = buf(4)
     redshift = buf(5)
     uv       = buf(6)
  endif
end subroutine bcasttime

subroutine bcastunits()
  use mpi
  use global
  implicit none

  integer, parameter :: nbuf=16

  real(8), dimension(nbuf) :: buf

  if (master) then
     buf(1)  = density0
     buf(2)  = length0
     buf(3)  = time0
     buf(4)  = velocity0
     buf(5)  = number0
     buf(6)  = chem0
     buf(7)  = photo0
     buf(8)  = temp0
     buf(9)  = temp0inv
     buf(10) = press0
     buf(11) = heat0
     buf(12) = cool0
     buf(13) = jeans0
     buf(14) = jeansfloor0
     buf(15) = kappa0
     buf(16) = meanfreepath0
  endif
  
  call mpi_bcast(buf,nbuf,mpi_double_precision,0,mpi_comm_world,mpierror)

  if (notmaster) then
     density0      = buf(1)
     length0       = buf(2)
     time0         = buf(3)
     velocity0     = buf(4)
     number0       = buf(5)
     chem0         = buf(6)
     photo0        = buf(7)
     temp0         = buf(8)
     temp0inv      = buf(9)
     press0        = buf(10)
     heat0         = buf(11)
     cool0         = buf(12)
     jeans0        = buf(13)
     jeansfloor0   = buf(14)
     kappa0        = buf(15)
     meanfreepath0 = buf(16)
  endif
end subroutine bcastunits

subroutine bcastcoef()
  use mpi
  use global
  implicit none
  
  call mpi_bcast(coef,ncool*2,mpi_double_precision,0,mpi_comm_world,mpierror)
end subroutine bcastcoef

subroutine reducestep()
  use mpi
  use global
  implicit none

  real(8) :: newdt

  call mpi_reduce(dt,newdt,1,mpi_double_precision,mpi_min,0, &
       mpi_comm_world,mpierror)
  if (master) dt = newdt
end subroutine reducestep

subroutine reducesumdouble(a)
  use mpi
  use global
  implicit none

  real(8), intent(inout) :: a
  real(8) :: tmp

  call mpi_reduce(a,tmp,1,mpi_double_precision,mpi_sum,0, &
       mpi_comm_world,mpierror)
  if (master) a = tmp
end subroutine reducesumdouble

subroutine reducemaxdouble(a)
  use mpi
  use global
  implicit none

  real(8), intent(inout) :: a
  real(8) :: tmp

  call mpi_reduce(a,tmp,1,mpi_double_precision,mpi_max,0, &
       mpi_comm_world,mpierror)
  if (master) a = tmp
end subroutine reducemaxdouble

subroutine reducemindouble(a)
  use mpi
  use global
  implicit none

  real(8), intent(inout) :: a
  real(8) :: tmp

  call mpi_reduce(a,tmp,1,mpi_double_precision,mpi_min,0, &
       mpi_comm_world,mpierror)
  if (master) a = tmp
end subroutine reducemindouble

subroutine allreducesumdouble(a)
  use mpi
  use global
  implicit none

  real(8) :: a

  call mpi_allreduce(a,a,1,mpi_double_precision,mpi_sum, &
       mpi_comm_world,mpierror)

end subroutine allreducesumdouble

subroutine allreducesumint(a)
  use mpi
  use global
  implicit none

  integer :: a
  integer :: tmp

  call mpi_allreduce(a,tmp,1,mpi_integer,mpi_sum, &
       mpi_comm_world,mpierror)

  a = tmp

end subroutine allreducesumint

subroutine bcastdouble(a)
  use mpi
  use global
  implicit none

  real(8)  :: a

  call mpi_bcast(a,1,mpi_double_precision,0,mpi_comm_world,mpierror)

end subroutine bcastdouble

subroutine bcastlogical(a)
  use mpi
  use global
  implicit none

  logical  :: a

  call mpi_bcast(a,1,mpi_logical,0,mpi_comm_world,mpierror)

end subroutine bcastlogical

subroutine bcastinteger(a)
  use mpi
  use global
  implicit none

  integer :: a

  call mpi_bcast(a,1,mpi_integer,0,mpi_comm_world,mpierror)

end subroutine bcastinteger

subroutine bcastarray(a,b)

  use mpi
  use global
  implicit none

  integer :: b
  real(8), dimension(b) :: a

  call mpi_bcast(a,b,mpi_double_precision,0,mpi_comm_world,mpierror)

end subroutine bcastarray

subroutine reducedoublearray(a,ared,n)

  use mpi
  use global
  implicit none

  real(8), dimension(n) :: a,ared
  integer :: n

  call mpi_reduce(a,ared,n,mpi_double_precision,mpi_sum, &
       0,mpi_comm_world,mpierror)

end subroutine reducedoublearray

subroutine reduceintarray(a,ared,n)

  use mpi
  use global
  implicit none

  integer, dimension(n) :: a,ared
  integer :: n

  call mpi_reduce(a,ared,n,mpi_integer,mpi_sum, &
       0,mpi_comm_world,mpierror)

end subroutine reduceintarray

subroutine barrier()

  use mpi
  use global
  implicit none

  call mpi_barrier(mpi_comm_world,mpierror)

end subroutine barrier

subroutine communicate(w,buf1,buf2,buf3,buf4)
  ! communicates the ghost cells for the primitive variables
  ! in a non-blocking way
  use mpi
  use global
  implicit none

  real(8), dimension(nw,lx2:mx2,ly2:my2,lz2:mz2) :: w
  real(8), dimension(nw,lx2:mx2,ly2:my2,2) :: buf1,buf2,buf3,buf4

  integer :: count,req1,req2,req3,req4
  integer, dimension(MPI_STATUS_SIZE) :: status

  count = size(buf1)

  ! send up and receive down
  if (notlast) then
     buf1(1:nw,lx2:mx2,ly2:my2,1:2) = w(1:nw,lx2:mx2,ly2:my2,mz0:mz)
     call mpi_isend(buf1,count,mpi_double_precision, &
          mpirank+1,0,mpi_comm_world,req1,mpierror)
     call mpi_irecv(buf4,count,mpi_double_precision, &
          mpirank+1,1,mpi_comm_world,req4,mpierror)
  else if (last .and. boundtype .eq. 2) then
     buf1(1:nw,lx2:mx2,ly2:my2,1:2) = w(1:nw,lx2:mx2,ly2:my2,mz0:mz)
     call mpi_isend(buf1,count,mpi_double_precision, &
          0,2,mpi_comm_world,req1,mpierror)
     call mpi_irecv(buf4,count,mpi_double_precision, &
          0,3,mpi_comm_world,req4,mpierror)
  endif
  ! send down and receive up
  if (notmaster) then
     buf3(1:nw,lx2:mx2,ly2:my2,1:2) = w(1:nw,lx2:mx2,ly2:my2,lz:lz0)
     call mpi_isend(buf3,count,mpi_double_precision, &
          mpirank-1,1,mpi_comm_world,req3,mpierror)
     call mpi_irecv(buf2,count,mpi_double_precision, &
          mpirank-1,0,mpi_comm_world,req2,mpierror)
  else if (master .and. boundtype .eq. 2 ) then
     buf3(1:nw,lx2:mx2,ly2:my2,1:2) = w(1:nw,lx2:mx2,ly2:my2,lz:lz0)
     call mpi_isend(buf3,count,mpi_double_precision, &
          mpisize-1,3,mpi_comm_world,req3,mpierror)
     call mpi_irecv(buf2,count,mpi_double_precision, &
          mpisize-1,2,mpi_comm_world,req2,mpierror)
  endif
  ! wait and copy buffer to array

  if (notlast) then
     call mpi_wait(req1,status,mpierror)
     call mpi_wait(req4,status,mpierror)
     w(1:nw,lx2:mx2,ly2:my2,mz1:mz2) = buf4(1:nw,lx2:mx2,ly2:my2,1:2)
  else if (last .and. boundtype .eq. 2) then
     call mpi_wait(req1,status,mpierror)
     call mpi_wait(req4,status,mpierror)
     w(1:nw,lx2:mx2,ly2:my2,mz1:mz2) = buf4(1:nw,lx2:mx2,ly2:my2,1:2)
  endif
  if (notmaster) then
     call mpi_wait(req2,status,mpierror)
     call mpi_wait(req3,status,mpierror)
     w(1:nw,lx2:mx2,ly2:my2,lz2:lz1) = buf2(1:nw,lx2:mx2,ly2:my2,1:2)
  else if (master.and. boundtype .eq. 2) then
     call mpi_wait(req2,status,mpierror)
     call mpi_wait(req3,status,mpierror)
     w(1:nw,lx2:mx2,ly2:my2,lz2:lz1) = buf2(1:nw,lx2:mx2,ly2:my2,1:2)
  endif

end subroutine communicate

subroutine gravcomm(phi,phibuf1,phibuf2,phibuf3,phibuf4)
  ! communicates the ghost cells for the gravitational potential
  ! in a non-blocking way
  use mpi
  use global
  implicit none

  complex(8), dimension(lx:mx,ly:my,lz:mz) :: phi
  real(8), dimension(lx:mx,ly:my) :: phibuf1,phibuf2,phibuf3,phibuf4

  integer :: count,req1,req2,req3,req4
  integer, dimension(MPI_STATUS_SIZE) :: status

  count = size(phibuf1)

  ! send up and receive down
  if (notlast) then
     phibuf1(lx:mx,ly:my) = dble(phi(lx:mx,ly:my,mz))
     call mpi_isend(phibuf1,count,mpi_double_precision, &
          mpirank+1,0,mpi_comm_world,req1,mpierror)
     call mpi_irecv(phibuf4,count,mpi_double_precision, &
          mpirank+1,1,mpi_comm_world,req4,mpierror)
  else ! last
     phibuf1(lx:mx,ly:my) = dble(phi(lx:mx,ly:my,mz))
     call mpi_isend(phibuf1,count,mpi_double_precision, &
          0,2,mpi_comm_world,req1,mpierror)
     call mpi_irecv(phibuf4,count,mpi_double_precision, &
          0,3,mpi_comm_world,req4,mpierror)
  endif
  ! send down and receive up
  if (notmaster) then
     phibuf3(lx:mx,ly:my) = dble(phi(lx:mx,ly:my,lz))
     call mpi_isend(phibuf3,count,mpi_double_precision, &
          mpirank-1,1,mpi_comm_world,req3,mpierror)
     call mpi_irecv(phibuf2,count,mpi_double_precision, &
          mpirank-1,0,mpi_comm_world,req2,mpierror)
     
  else ! master
     phibuf3(lx:mx,ly:my) = dble(phi(lx:mx,ly:my,lz))
     call mpi_isend(phibuf3,count,mpi_double_precision, &
          mpisize-1,3,mpi_comm_world,req3,mpierror)
     call mpi_irecv(phibuf2,count,mpi_double_precision, &
          mpisize-1,2,mpi_comm_world,req2,mpierror)
  endif
  
  ! wait
  call mpi_wait(req1,status,mpierror)
  call mpi_wait(req2,status,mpierror)
  call mpi_wait(req3,status,mpierror)
  call mpi_wait(req4,status,mpierror)   

end subroutine gravcomm
