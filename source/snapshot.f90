subroutine startupscreen()
  use global
  implicit none

  open (unit=10, file="data/log", access="sequential", &
       status="unknown", position="append")

  write (* ,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write (10,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

  write (* ,*) "WELCOME TO EVORA"
  write (10,*) "WELCOME TO EVORA"

  write (* ,*) "todays menue:"
  write (10,*) "todays menue:"

  if (secondorder) then
     write (* ,*) '  second-order'
     write (10,*) '  second-order'
  else   
     write (* ,*) '  first-order'
     write (10,*) '  first-order'
  endif

  if (atomic) then
     write (* ,*) '  gamma forced to 5/3'
     write (10,*) '  gamma forced to 5/3'
  endif

  if (supercomoving) then
     write (* ,*) '  supercomoving coordinates'
     write (10,*) '  supercomoving coordinates'
  endif

  if (dualenergy) then
     write (* ,*) '  dual energy formalism'
     write (10,*) '  dual energy formalism'
  endif

  if (chemistry) then
     if (nq .ne. 12) then
        print *,"ERROR compile with -DCHEMFLAG for chemisty"
        call cleanstop
     endif

     if (noncie) then
        write (* ,*) '  non-CIE chemistry'
        write (10,*) '  non-CIE chemistry'
        if (nw .ne. 11 .or. nf .ne. 12) then
           print *,"ERROR compile with -DNONCIEFLAG for noncie"
           call cleanstop
        endif
     else
        write (* ,*) '  CIE chemistry'
        write (10,*) '  CIE chemistry'
     endif

     if (cooling) then
        write (* ,*) '  with cooling'
        write (10,*) '  with cooling'
     endif
  endif

  if (conduction) then
     write (* ,*) '  thermal conduction'
     write (10,*) '  thermal conduction'
  endif

  if (gravity) then
     write (* ,*) '  with gravitation'
     write (10,*) '  with gravitation'
     if (dx .ne. dy .or. dx .ne. dz) then
        print *,'ERROR non regular grid'
        call cleanstop()
     endif
     if (boundtype .ne. 2) then
        print *,'ERROR non periodic boundary'
        call cleanstop()
     endif
  endif

  write (*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write (10,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

  write (*,*) 'time integration ...'
  write (10,*) 'time integration ...'

  if (supercomoving) then           
     write (* ,*) & 
          "No.  k     time        timestep        a            z"
     write (10,*) & 
          "No.  k     time        timestep        a            z"
  else
     write (* ,*) "No.  k     time        timestep                percentage"
     write (10,*) "No.  k     time        timestep                percentage"
  endif 

  write (* ,*) "---------------------------------------------------------"
  write (10,*) "---------------------------------------------------------"

  close (10)

end subroutine startupscreen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine screen(k1,k)
  use global
  implicit none

  integer :: k1,k
  real(8) :: tpro

  open (unit=10, file="data/log", access="sequential", &
       status="unknown", position="append")
  if (supercomoving) then           
     write (*,'(I4,I6,E12.4,E12.4,F12.4,F12.2)') & 
          k1,k,t,dt,cosmo,redshift
     write (10,'(I4,I6,E12.4,E12.4,F12.4,F12.2)') & 
          k1,k,t,dt,cosmo,redshift
  else
     tpro = tsnap / nsnap * 100
     write (*,'(I4,I6,E12.4,E12.4,11X,F12.1,A1)') k1,k,t,dt,tpro,'%'
     write (10,'(I4,I6,E12.4,E12.4,11X,F12.1,A1)') k1,k,t,dt,tpro,'%'
  endif     
  close(10)

  k1 = k1 + 1

end subroutine screen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine evolution(k2,q)
  use global
  implicit none
 
  integer, intent(inout)                                :: k2
  real(8), intent(in) , dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(8) :: press,temp
  real(8) :: totalmass,totalenergy,totalentropy
  integer :: ix,iy,iz

  totalmass = sum(q(1,lx:mx,ly:my,lz:mz)) * volume
  totalenergy = sum(q(5,lx:nx,ly:my,lz:mz)) * volume
  totalentropy = sum(q(6,lx:nx,ly:my,lz:mz)) * volume

  if (parallel) then
     call reducesumdouble(totalmass)
     call reducesumdouble(totalenergy)
     call reducesumdouble(totalentropy)
  endif

  if (master) then
     open (unit=32, file="data/evo", access="sequential", &
          status="unknown", position="append")
     if (k2 .eq. 0) write (32,'(A,A15,23A16)') "#","t","dt","a","z",&
          "totalmass","totalenergy","totalentropy"
     write (32,'(24E16.6)') t,dt,cosmo,redshift, &
          totalmass,totalenergy,totalentropy
  endif

  call barrier()

  k2 = k2 + 1

end subroutine evolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine snapshot(k3,q)
  use global
  implicit none

  integer :: k3
  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q

  real(4), dimension(:,:), allocatable :: tmp

  real(8) :: press,temp,x,y,z,templ,tempr,ner,nel,j

  real(8), dimension(nchem) :: chemrate
  real(8), dimension(ncool) :: coolrate

  real(8) :: lambda,gamma

  character(len=127) :: name,filename

  integer :: x1,y1,z1,ierror,recordsize

  integer :: i,ix,iy,iz

  write(name,'(I0.3)') k3
   
  x1 = max(nx/2,1)
  y1 = max(ny/2,1)
  z1 = max(nz/2,1)

  recordsize = nx * ny * 4
  if (master .and. bov) then
     call writebovhead(name,"dens","snap/")
     call writebovhead(name,"temp","snap/")
  endif

  if (bov) allocate(tmp(1:nx,1:ny))

  ! loop over processes
  i = 0
  do 
     if (mpirank .eq. i) then
        
        if (cut) then
           ! x-cut
           if (lz .le. z1 .and. mz .ge. z1) then
              ! x-cut
              filename = "snap/"//trim(name)//".x"
              open (unit=20, file=trim(filename),form='formatted', &
                   access='sequential', iostat=ierror)
              call checkfile(ierror,filename)
              iz = z1
              iy = y1
              do ix=1,nx
                 call coordinates(ix,iy,iz,x,y,z)
                 call pressure(q(:,ix,iy,iz),press)
                 call press2temp(press,q(1,ix,iy,iz),q(7:nq,ix,iy,iz),temp)

                 write (20,'(5E16.6)') x,q(1,ix,iy,iz),temp, &
                      q(2,ix,iy,iz)/q(1,ix,iy,iz)*velocity0,press*press0
              enddo
              close(20)

              ! x0-cut
              if (onedimension .eqv. .false.) then
                 filename = "snap/"//trim(name)//".x0"
                 open (unit=20, file=trim(filename),form='formatted', &
                      access='sequential', iostat=ierror)
                 call checkfile(ierror,filename)
                 iz = z1
                 iy = 1
                 do ix=1,nx
                    call coordinates(ix,iy,iz,x,y,z)
                    call pressure(q(:,ix,iy,iz),press)
                    call press2temp(press,q(1,ix,iy,iz),q(7:nq,ix,iy,iz),temp)

                    write (20,'(5E16.6)') x,q(1,ix,iy,iz),temp, &
                         q(2,ix,iy,iz)/q(1,ix,iy,iz)*velocity0,press*press0
                 enddo
                 close(20)
              endif
           endif
        endif
        
        if (slice) then
           ! xy-slice
           if (lz .le. z1 .and. mz .ge. z1) then
              filename = "snap/"//trim(name)//".xy"
              open (unit=20, file=trim(filename),form='formatted', &
                   access='sequential', iostat=ierror)
              call checkfile(ierror,filename)
              iz = z1
              do ix=1,nx
                 do iy=1,ny
                    call coordinates(ix,iy,iz,x,y,z)
                    call pressure(q(:,ix,iy,iz),press)
                    call press2temp(press,q(1,ix,iy,iz),q(7:nq,ix,iy,iz),temp)
                    write (20,'(5E16.6)') x,y,q(1,ix,iy,iz),temp,press*press0
                 enddo
                 write (20,*)
              enddo
              close(20)
           endif

           ! xz-slice
           filename = "snap/"//trim(name)//".xz"
           open (unit=20, file=trim(filename),form='formatted', &
                access='sequential', iostat=ierror, position="append")
           call checkfile(ierror,filename)
           iy = y1
           do iz=lz,mz
              do ix=1,nx
                 call coordinates(ix,iy,iz,x,y,z)
                 call pressure(q(:,ix,iy,iz),press)
                 call press2temp(press,q(1,ix,iy,iz),q(7:nq,ix,iy,iz),temp)
                 write (20,'(5E16.6)') x,z,q(1,ix,iy,iz),temp,press*press0
              enddo
              write (20,*)
           enddo
           close(20)            
        endif
        
        if (bov) then
           ! density bov file
           filename = "snap/"//trim(name)//".dens.bovdata"
           open (unit=20, file=trim(filename),form='unformatted', &
                access='direct', iostat=ierror, recl=recordsize)
           call checkfile(ierror,filename)
           do iz=lz,mz
              write (20,rec=iz) real(q(1,1:nx,1:ny,iz))
           enddo
           close(20)
           ! temperature bov file
           filename = "snap/"//trim(name)//".temp.bovdata"
           open (unit=20, file=trim(filename),form='unformatted', &
                access='direct', iostat=ierror, recl=recordsize)
           call checkfile(ierror,filename)
           do iz=lz,mz
              do iy=1,ny
                 do ix=1,nx
                    call pressure(q(:,ix,iy,iz),press)
                    call press2temp(press,q(1,ix,iy,iz),q(7:nq,ix,iy,iz),temp)
                    tmp(ix,iy) = real(temp)
                 enddo
              enddo
              write (20,rec=iz) tmp
           enddo
           close(20)
        endif

     endif
     
     call barrier

     if (i .eq. mpisize) then
        k3 = k3 + 1
        return
     else
        i = i + 1
     endif

  enddo

  if (bov) deallocate(tmp)

end subroutine snapshot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine endscreen
  use global
  implicit none
  
  open (unit=10, file="data/log", access="sequential", &
       status="unknown", position="append")

  write (*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write (10,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

  write (* ,*) 'done!' 
  write (10,*) 'done!'

  write (*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write (10,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

  close (10)

end subroutine endscreen
