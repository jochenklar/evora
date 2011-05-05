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
