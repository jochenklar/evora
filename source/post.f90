program post
  use global
  implicit none

  character(80) :: mode,infile,option

  ! read mode
  call get_command_argument(1,mode)
  call get_command_argument(2,infile)
  call get_command_argument(3,option)
  if (len_trim(mode) .eq. 0 .or. len_trim(infile) .eq. 0) then
     print *,"ERROR to few arguments"
     print *,"USAGE bin/post <mode> <snapshot> <option>"
     stop
  endif

  ! check mode
  if (.not.(trim(mode) .eq. "cs" &
       .or. trim(mode) .eq. "bov" &
       .or. trim(mode) .eq. "phase" &
       .or. trim(mode) .eq. "ps")) then
     print *,"wrong mode"
     print *,"choose cs bov phase line"
     print *,"USAGE bin/post <mode> <snapshot> <option>"
     stop
  endif

  ! get the parameter
  call rereadparameter()
  call compparameter()
  call readunits(trim(infile))

  ! chemical stuff
  chiH = rhoH / mH
  chiHe = rhoHe / mHe
  if (background) call compbackground
  call compcoef()

  select case (trim(mode))
  case("cs")
     print *,"cuts and slices"
     call cuts_and_slices(infile)
  case("bov")
     print *,"bov"
     call bovs(infile,option)
  case("phase")
     print *,"phase diagram"
     call phasediagram(infile,option)
  case("ps")
     print *,"power spectrum"
     call powerspectrum(infile,option)
  end select

end program post

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! cuts and slices                                                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cuts_and_slices(infile)
  use global
  implicit none

  character(*) :: infile

  real(8), dimension(:,:,:,:), allocatable :: q
  character(80) :: filename
  integer :: j,k,ix,iy,iz,ierror

  integer :: i,i1,i2,ncut,nslice

  logical :: cut_x,cut_y,cut_z,cut_xy,cut_xyz,slice_xy,slice_xz,slice_yz
  real(8), dimension(:,:), allocatable :: cut_x_data,cut_y_data,cut_z_data,cut_xy_data,cut_xyz_data
  real(8), dimension(:,:,:), allocatable :: slice_xy_data,slice_xz_data,slice_yz_data

  real(8) :: xmid,ymid,zmid

  ! determine 'middle' cell
  xmid = max(nx/2,1) ; ymid = max(ny/2,1); zmid = max(nz/2,1);

  ! init logical variables
  cut_x=.false.;cut_y=.false.;cut_z=.false.;cut_xy=.false.;cut_xyz=.false.
  slice_xy=.false.;slice_xz=.false.;slice_yz=.false.

  ! decide what to do
  cut_x = .true.
  if (.not. onedimension) then
     cut_y = .true.
     slice_xy = .true.
     if (nx .eq. ny) cut_xy = .true.
     if (.not. twodimension) then
        cut_z = .true.
        slice_xz = .true.
        slice_yz = .true.
        if (nx .eq. ny .and. ny .eq. nz) cut_xyz = .true.
     endif
  endif

  ! allocate arrays and construct coordinates
  if (cut_x) then
     allocate(cut_x_data(nstate+1,nx))
     cut_x_data(1,:) = xmin + ((/(k,k=1,nx)/)-0.5)*dx
  endif
  if (cut_y) then
     allocate(cut_y_data(nstate+1,ny))
     cut_y_data(1,:) = ymin + ((/(k,k=1,ny)/)-0.5)*dy
  endif
  if (cut_z) then
     allocate(cut_z_data(nstate+1,nz))
     cut_z_data(1,:) = zmin + ((/(k,k=1,nz)/)-0.5)*dz
  endif
  if (cut_xy) then
     allocate(cut_xy_data(nstate+1,nx))
     cut_xy_data(1,:) = sqrt(2.0) * (xmin + ((/(k,k=1,nx)/)-0.5)*dx)
  endif
  if (cut_xyz) then
     allocate(cut_xyz_data(nstate+1,nx))
     cut_xyz_data(1,:) = sqrt(3.0) * (xmin + ((/(k,k=1,nx)/)-0.5)*dx)
  endif
  if (slice_xy) then
     allocate(slice_xy_data(nstate+2,nx,ny))
     do iy=1,ny
        do ix=1,nx
             slice_xy_data(1,ix,iy) = xmin + (ix-0.5) * dx
             slice_xy_data(2,ix,iy) = ymin + (iy-0.5) * dy
        enddo
     enddo
  endif
  if (slice_xz) then
     allocate(slice_xz_data(nstate+2,nx,nz))
     do iz=1,nz
        do ix=1,nx
           slice_xz_data(1,ix,iz) = xmin + (ix-0.5) * dx
           slice_xz_data(2,ix,iz) = zmin + (iz-0.5) * dz
        enddo
     enddo
  endif
  if (slice_yz) then
     allocate(slice_yz_data(nstate+2,ny,nz))
     do iz=1,nz
        do iy=1,ny
           slice_yz_data(1,iy,iz) = ymin + (iy-0.5) * dy
           slice_yz_data(2,iy,iz) = zmin + (iz-0.5) * dz
        enddo
     enddo
  endif

  ! gather data
  mz = 0 ; j = 0
  do while ((mz .eq. nz) .eqv. .false.)     
     ! open the datafile and read header
     call opendata(j,infile)

     ! allocate arrays
     allocate(q(nq,1:nx,1:ny,lz:mz))

     ! read array
     read (10) q
     close(10)

     ! loop through the grid
     do iz=lz,mz
        do iy=1,ny
           do ix=1,nx  
              ! cut_x
              if (cut_x .and. iy .eq. ymid .and. iz .eq. zmid) &
                 call compstate(q(:,ix,iy,iz),cut_x_data(2:nstate+1,ix))
              ! cut_y
              if (cut_y .and.  ix .eq. xmid .and. iz .eq. zmid) &
                 call compstate(q(:,ix,iy,iz),cut_y_data(2:nstate+1,iy))
              ! cut_z
              if (cut_z .and. ix .eq. xmid .and. iy .eq. ymid) &
                 call compstate(q(:,ix,iy,iz),cut_z_data(2:nstate+1,iz))
              ! cut_xy
              if (cut_xy .and. iz .eq. zmid .and. ix .eq. iy) &
                 call compstate(q(:,ix,iy,iz),cut_xy_data(2:nstate+1,ix))
              ! cut_xyz
              if (cut_xyz .and. ix .eq. iy .and. ix .eq. iz) &
                 call compstate(q(:,ix,iy,iz),cut_xyz_data(2:nstate+1,ix))
              ! slice_xy
              if (slice_xy .and. iz .eq. zmid) &
                 call compstate(q(:,ix,iy,iz),slice_xy_data(3:nstate+2,ix,iy))
              ! slice_xz
              if (slice_xz .and. iy .eq. ymid) &
                 call compstate(q(:,ix,iy,iz),slice_xz_data(3:nstate+2,ix,iz))
              ! slice_yz
              if (slice_yz .and. ix .eq. xmid) &
                 call compstate(q(:,ix,iy,iz),slice_yz_data(3:nstate+2,iy,iz))
           enddo
        enddo
     enddo

     ! get rid of the arrays
     deallocate(q)

     j = j + 1
  enddo

  ! write down cuts and slices
  if (cut_x)    call writecut(cut_x_data,nx,nstate+1,infile,"x")
  if (cut_y)    call writecut(cut_y_data,ny,nstate+1,infile,"y")
  if (cut_z)    call writecut(cut_z_data,nz,nstate+1,infile,"z")
  if (cut_xy)   call writecut(cut_xy_data,nx,nstate+1,infile,"xy")
  if (cut_xyz)  call writecut(cut_xy_data,nx,nstate+1,infile,"xyz")
  if (slice_xy) call writeslice(slice_xy_data,nx,ny,nstate+2,infile,"xy")
  if (slice_xz) call writeslice(slice_xz_data,nx,nz,nstate+2,infile,"xz")
  if (slice_yz) call writeslice(slice_yz_data,ny,nz,nstate+2,infile,"yz")

  ! clean up
  if (cut_x)    deallocate(cut_x_data)
  if (cut_y)    deallocate(cut_y_data)
  if (cut_z)    deallocate(cut_z_data)
  if (cut_xy)   deallocate(cut_xy_data)
  if (cut_xyz)  deallocate(cut_xyz_data)
  if (slice_xy) deallocate(slice_xy_data)
  if (slice_xz) deallocate(slice_xz_data)
  if (slice_yz) deallocate(slice_yz_data)

end subroutine cuts_and_slices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! bov                                                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bovs(infile,var)
  use global
  implicit none

  character(*) :: infile,var

  real(8), dimension(:,:,:,:), allocatable :: q
  real, dimension(:,:), allocatable :: tmp
  character(len=80) :: filename
  integer :: ierror,recordsize
  real(8) :: press,temp
  integer :: i,ix,iy,iz

  ! get one optional argument
  if ((trim(var) .eq. "dens" .or. trim(var) .eq. "temp") &
       .eqv. .false.) then
     print *,"ERROR 3nd argument has to be 'dens' or 'temp'!",var
     stop
  endif

  if (ny .eq. 1 .and. nz .eq. 1) then
     print *,"PROBLEM bov makes no sense in 1D!"
     stop
  endif

  recordsize = nx * ny * 4
  if (trim(var) .eq. "temp") allocate(tmp(1:nx,1:ny))

  ! writing the headerfile
  call writebovhead(infile,var,"out/")

  ! open direct-access file
  filename = "out/"//trim(infile)//"."//trim(var)//".bovdata"

  open (unit=20, file=trim(filename),form='unformatted', &
       access='direct', iostat=ierror, recl=recordsize)
  call checkfile(ierror,filename)

  ! gathering the data
  mz = 0
  i = 0
  do while ((mz .eq. nz) .eqv. .false.)    

     ! open the datafile and read header
     call opendata(i,infile)

     ! allocate arrays
     allocate(q(nq,1:nx,1:ny,lz:mz))

     ! read arrays
     read (10) q
     close(10)

     ! write slice after slice in file
     if (trim(var) .eq. "dens") then ! DENSITY
        do iz=lz,mz
           write (20,rec=iz) real(q(1,1:nx,1:ny,iz))
        enddo
     elseif (trim(var) .eq. "temp") then ! TEMPERATURE
        do iz=lz,mz
           do iy=1,ny
              do ix=1,nx
                 call pressure(q(:,ix,iy,iz),press)
                 call press2temp(press,q(1,ix,iy,iz),q(7:nq,ix,iy,iz),temp)
                 tmp(ix,iy) = real(temp)   
              enddo  
           enddo
           write (20,rec=iz) tmp(1:nx,1:ny)
        enddo
     endif

     deallocate(q)

     i = i + 1
  enddo
  close(20)

  if (trim(var) .eq. "temp") deallocate(tmp)

end subroutine bovs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! phasediagramm                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phasediagram(infile,option)
  use global
  implicit none

  character(*) :: infile,option

  real(8), dimension(:,:,:,:), allocatable :: q
  character(80) :: filename
  integer :: j,k,ix,iy,iz,ierror
  real(8) :: press,temp

  integer :: agrid,bgrid,ia,ib
  real(8), dimension(:,:,:), allocatable :: phase ! vol,mass
  real(8) :: aa,bb,adelta,bdelta,amin,bmin,amax,bmax
  real(8) :: totalmass,bin1,bin2,cum1,cum2
  real(8) :: photo,whim,cond,hot
  real(8), dimension(4,2) :: dave ! photo,condensed,whim,hot;vol,mass
  real(8), dimension(2)   :: overdens ! vol,mass

  if (len_trim(option) .eq. 0) then
     print *,"ERROR no. of gridpoints as option"
     stop
  endif
  read(option,*) agrid
  if (.not.(agrid .gt. 0 .and. agrid .le. 1024)) then
     print *,"ERROR no. of gridpoints makes no sense"
     stop
  endif

  bgrid = agrid 

  allocate(phase(agrid,bgrid,2))

  amin = -3.0 ; amax = 4.0 ; adelta = (amax - amin)/dble(agrid) ! dens
  bmin = -2.0 ; bmax = 8.0 ; bdelta = (bmax - bmin)/dble(bgrid) ! temp
  phase = 0.0 ; totalmass = 0.0

  ! gather data
  mz = 0 ; j = 0
  do while ((mz .eq. nz) .eqv. .false.)     
     ! open the datafile and read header
     call opendata(j,infile)
     
     ! allocate arrays
     allocate(q(nq,1:nx,1:ny,lz:mz))
     
     ! read array
     read (10) q
     close(10)

     ! loop through the grid
     do iz=lz,mz
        do iy=1,ny
           do ix=1,nx  
              totalmass = totalmass + q(1,ix,iy,iz)
              call pressure(q(:,ix,iy,iz),press)
              call press2temp(press,q(1,ix,iy,iz),q(7:nq,ix,iy,iz),temp)
              aa = (log10(q(1,ix,iy,iz)) - amin) / adelta
              bb = (log10(temp)          - bmin) / bdelta
              ia = floor(aa)
              ib = floor(bb)
              if (aa .ge. 1 .and. aa .le. agrid &
                   .and. bb .ge. 1 .and. bb .le. bgrid) then 
                 phase(ia,ib,1) = phase(ia,ib,1) + n3inv
                 phase(ia,ib,2) = phase(ia,ib,2) + q(1,ix,iy,iz)
              endif
           enddo
        enddo
     enddo

     ! get rid of the arrays
     deallocate(q)
     j = j + 1
  enddo

  ! finalize phase diagram
  phase(ia,ib,2) = phase(ia,ib,2) * density0 * length0**3.0 / solarmass * n3inv
 
  ! write phase diagramm
  filename = "out/"//trim(infile)//".phase"
  open(unit=10,file=trim(filename),form='formatted')
  do ia=1,agrid
     do ib=1,bgrid
        aa = amin + (ia+0.5)*adelta
        bb = bmin + (ib+0.5)*bdelta
        write(10,'(4E16.8)') aa,bb,phase(ia,ib,1),phase(ia,ib,2)
     enddo
     write(10,*)
  enddo
  close(10)

end subroutine phasediagram

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! powerspectrum                                                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine powerspectrum(infile,option)
  use global
  implicit none

  character(*) :: infile,option

  real(8), dimension(:,:,:,:), allocatable :: q
  character(80) :: filename
  integer :: i,j,ix,iy,iz,ierror

  integer, parameter :: nk=70
  real(8), parameter :: logstart=-2.05,logdelta=0.1

  real(8), dimension(nk) :: power
  integer, dimension(nk) :: number
  real(8) :: k,kx,ky,kz,kx2,ky2,kz2,kf

  complex(8), dimension(:,:,:), allocatable :: dens

  kf = 2.0 * pi / box

  ! allocate density array
  allocate(dens(nx,ny,nz))

  ! gather data
  mz = 0 ; j = 0
  do while ((mz .eq. nz) .eqv. .false.)     
     ! open the datafile and read header
     call opendata(j,infile)
     
     ! allocate arrays
     allocate(q(nq,1:nx,1:ny,lz:mz))
     
     ! read array
     read (10) q
     close(10)
     
     ! copy in density array
     dens(1:nx,1:ny,lz:mz) = q(1,:,:,:)

     ! get rid of the arrays
     deallocate(q)
     
     j = j + 1
  enddo

  ! fft the density field
  call fft_c2c_power(dens)
  dens = dens / dble(n3)

  power = 0.0
  number = 0.0

  ! obtain the powerspectrum
  do iz=1,nz
     call kvalue(iz,nz,nyqz,kz)
     kz2 = kz * kz
     do iy=1,ny
        call kvalue(iy,ny,nyqy,ky)
        ky2 = ky * ky
        do ix=1,nx
           call kvalue(ix,nx,nyqx,kx)
           kx2 = kx * kx
           k = sqrt(kx2 + ky2 + kz2) * kf
           if (k .gt. 0.0) then
              ! binning
              i = floor((log10(k) - logstart) / logdelta)
              if (i .ge. 1 .and. i .le. nk) then
                 power(i) = power(i) + abs(dens(ix,iy,iz))**2
                 number(i) = number(i) + 1
              endif
           endif
        enddo
     enddo
  enddo

  do i=1,nk
     if (number(i) .gt. 0) then
        power(i) = power(i) / number(i)
     endif
  enddo

  if (onedimension) then
     power = power * box
  else if (twodimension) then 
     power = power * box**2
  else
     power = power * box**3
  endif

  filename = "out/"//trim(infile)//".power"
  open (unit=25, file=trim(filename), &
       access="sequential",status="unknown", form="formatted",iostat=ierror)
  call checkfile(ierror,filename)
  write (25,'(A,A15,26A16)') "#","k","power"
  do i=1,nk
     k = 10**(logstart + (i+0.5) * logdelta)
     write (25,'(2E16.6)') k,power(i)
  enddo
  close(25)

  deallocate(dens)

end subroutine powerspectrum










