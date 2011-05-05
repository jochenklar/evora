subroutine readparameter()
  use global
  implicit none

  namelist /para/ iconly,hydro,secondorder,solver,limiter,cfl,g,atomic, &
       nx,ny,nz,xmin,xmax,ymin,ymax,zmin,zmax,boundtype, &
       cut,slice,bov,nsnap,gravity,tstart,tend, &
       supercomoving,c_cosmo,omega0,lambda0,fb,hubble0,box,zstart,zend, &
       chemistry,noncie,background,bg,j0,cooling,cooldt,c_cooldt, &
       tempfloor,temp_fl, &
       c_chem,substep,rhoH,rhoHe,jeansfloor,c_jeans, &
       conduction,c_cond,simplekappa,dualenergy,de_thresh,dump_redshift

  integer :: ierror

  ! set the parameter to default values
  ! hydro
  iconly = .false.
  hydro = .true.
  secondorder = .true. ; solver = 2 ; limiter = 1 
  cfl = 0.5 ; g = 1.4 ; atomic = .false.
  ! grid
  nx = 32 ; ny = 32 ; nz = 32
  xmin = 0.0 ; xmax = 1.0
  ymin = 0.0 ; ymax = 1.0
  zmin = 0.0 ; zmax = 1.0
  boundtype = 2
  ! output
  cut = .false. ; slice = .false. ; bov = .false.
  nsnap = 100 
  ! gravity
  gravity = .false.
  ! time
  tstart = 0.0 ; tend= 1.0
  ! cosmology
  supercomoving = .false. ; c_cosmo = 0.1
  omega0  = 1.0 ;  lambda0 = 0.0 ; hubble0 = 50.0
  fb = 1.0 ; box = 1.0
  zstart  = 99.0 ; zend    = 0.0
  ! chemistry, cooling and heating
  chemistry = .false. ; noncie = .false.
  cooling = .false. ; 
  background = .false. ; bg = 0 ; j0 = 1d-21
  c_chem = 0.1 ; substep = 0.1 
  rhoH = 0.76 ; rhoHe = 0.24
  ! global cooling timestep
  cooldt = .false. ; c_cooldt = 0.0
  ! temperature floor
  tempfloor = .false. ; temp_fl = 0.0 
  ! jeans floor
  jeansfloor = .false. ; c_jeans = 0.0
  ! thermal conduction
  conduction = .false. ; c_cond = 0.1; simplekappa = .false.
  ! dualenergy
  dualenergy = .false. ; de_thresh = 0.0
  ! dumps
  dump_redshift = -1.0

  ! read the parameter
  call get_command_argument(1,parameterfile)
  if (len_trim(parameterfile) .eq. 0) then
     print *,"ERROR give command line argument!"
     call cleanstop
  endif
  open (unit=10, file=parameterfile, status="old", iostat=ierror)
  if (ierror .ne. 0) then
     print *,"ERROR file not found: ",parameterfile
     call cleanstop()
  endif
  read (10,nml=para)
  close(10)

  ! some sanity
  if (chemistry .eqv. .false.) then
     background = .false.
     cooling = .false.
  end if
  if (cooling .eqv. .false.) cooldt = .false.

end subroutine readparameter

subroutine bcastparameter()
  use global
  implicit none

  ! distribute parameter
  if (parallel) then
     call bcastlogical(iconly)
     call bcastlogical(hydro)
     call bcastlogical(secondorder) 
     call bcastinteger(solver)     
     call bcastinteger(limiter)
     call bcastdouble(cfl)
     call bcastdouble(g)
     call bcastlogical(atomic)

     call bcastinteger(nx) ; call bcastinteger(ny) ; call bcastinteger(nz)
     call bcastdouble(xmin) ; call bcastdouble(xmax)
     call bcastdouble(ymin) ; call bcastdouble(ymax)
     call bcastdouble(zmin) ; call bcastdouble(zmax)
     call bcastinteger(boundtype)

     call bcastlogical(cut)
     call bcastlogical(slice)
     call bcastlogical(bov)
     call bcastinteger(nsnap)

     call bcastlogical(gravity)

     call bcastdouble(tstart) ; call bcastdouble(tend)

     call bcastlogical(supercomoving) 
     call bcastdouble(omega0)
     call bcastdouble(lambda0)
     call bcastdouble(hubble0)
     call bcastdouble(fb) ; call bcastdouble(box)
     call bcastdouble(zstart) ; call bcastdouble(zend)

     call bcastlogical(chemistry) ; call bcastlogical(noncie)
     call bcastlogical(cooling) ; call bcastlogical(background)
     call bcastinteger(bg) ; call bcastdouble(j0)

     call bcastdouble(c_chem); call bcastdouble(substep)
     call bcastdouble(rhoH) ; call bcastdouble(rhoHe)
     call bcastlogical(cooldt) ; call bcastdouble(c_cooldt)
     call bcastlogical(tempfloor) ; call bcastdouble(temp_fl)
     call bcastlogical(jeansfloor) ; call bcastdouble(c_jeans)

     call bcastlogical(conduction) ; call bcastdouble(c_cond)
     call bcastlogical(simplekappa)

     call bcastlogical(dualenergy) ; call bcastdouble(de_thresh)

     call bcastarray(dump_redshift,ndump)
  endif

end subroutine bcastparameter

subroutine compparameter()
  use global
  implicit none

  ! force gamma
  if (atomic) then
     g = 5.0 / 3.0
  endif
  ! one dimension ? two dimensions ?
  if (ny .eq. 1 .and. nz .eq. 1) then
     onedimension = .true.
  else
     onedimension = .false.
  endif
  if (nz .eq. 1) then
     twodimension = .true.
  else
     twodimension = .false.
  endif
  ! grid stuff 
  n3    = nx * ny * nz
  n3inv = 1.0 / n3
  nx1 = nx + 1
  ny1 = ny + 1
  nz1 = nz + 1
  nx2 = nx + 2
  ny2 = ny + 2
  nz2 = nz + 2
  dx = (xmax - xmin) / nx
  dy = (ymax - ymin) / ny
  dz = (zmax - zmin) / nz
  if (onedimension) then
     ! enforce regulagrid with 1D
     dy   = dx
     dz   = dx
     ymin = 0.0
     ymax = dy
     zmin = 0.0
     zmax = dz
     dy = (ymax - ymin) / ny
     dz = (zmax - zmin) / nz
  endif
  if (twodimension) then
     ! enforce regulagrid with 2D
     dz   = dx
     zmin = 0.0
     zmax = dz
     dz = (zmax - zmin) / nz
  endif
  dx2 = 0.5 * dx
  dy2 = 0.5 * dy
  dz2 = 0.5 * dz
  dxinv = 1.0 / dx
  dyinv = 1.0 / dy
  dzinv = 1.0 / dz
  dx2inv = 0.5 / dx
  dy2inv = 0.5 / dy
  dz2inv = 0.5 / dz
  x0 = xmin + dx / 2
  y0 = ymin + dy / 2
  z0 = zmin + dz / 2
  volume = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
  ! some gamma dependent things
  g1    = g - 1.0
  g1inv = 1.0 / g1
  g2    = g - 2.0
  gg1   = g / g1
  g35   = 3.0 * g - 5.0
  ! time
  t = tstart
  cosmo = 1.0
  !grav
  if (gravity) then
     if (supercomoving) then
        grav = 1.5 * omega0
     else
        grav = fourpig
     endif
  endif
  ! double energy
  if (dualenergy) then
     de_thresh1 = 1.0 - de_thresh
  else
     de_thresh1 = 1.0
  endif
  ! fundamental fractions
  third     = 1.0 / 3.0
  twothird  = 2.0 / 3.0
  fourthird = 4.0 / 3.0
  ! for the poison solver
  h2   = dx ** 2.0 / 4.0
  nyqx = max(nx/2,1)
  nyqy = max(ny/2,1)
  nyqz = max(nz/2,1)
  kfx  = dble(2.0 * pi / nx)
  kfy  = dble(2.0 * pi / ny)
  kfz  = dble(2.0 * pi / nz)
  kfx2 = kfx * 0.5
  kfy2 = kfy * 0.5
  kfz2 = kfz * 0.5

end subroutine compparameter

subroutine writeparameter()
  use global
  implicit none

  integer :: ierror

  open (unit=10, file="data/parameter", form="formatted", iostat=ierror)
  if (ierror .ne. 0) then
     print *,"ERROR file not found: data/parameter"
     call cleanstop()
  endif

  write (10,'("iconly =        ",L16)') iconly
  write (10,'("hydro =         ",L16)') hydro
  write (10,'("secondorder =   ",L16)') secondorder
  write (10,'("solver =        ",I16)') solver
  write (10,'("limiter =       ",I16)') limiter
  write (10,'("cfl =           ",E16.8)') cfl
  write (10,'("g =             ",E16.8)') g
  write (10,'("atomic =        ",L16)') atomic
  write (10,'("nx =            ",I16)') nx
  write (10,'("ny =            ",I16)') ny
  write (10,'("nz =            ",I16)') nz
  write (10,'("xmin =          ",E16.8)') xmin
  write (10,'("xmax =          ",E16.8)') xmax
  write (10,'("ymin =          ",E16.8)') ymin
  write (10,'("ymax =          ",E16.8)') ymax
  write (10,'("zmin =          ",E16.8)') zmin
  write (10,'("zmax =          ",E16.8)') zmax
  write (10,'("boundtype =     ",I16)') boundtype
  write (10,'("cut =           ",L16)') cut
  write (10,'("slice =         ",L16)') slice
  write (10,'("bov =           ",L16)') bov
  write (10,'("nsnap =         ",I16)') nsnap
  write (10,'("gravity =       ",L16)') gravity
  write (10,'("tstart =        ",E16.8)') tstart
  write (10,'("tend =          ",E16.8)') tend
  write (10,'("supercomoving = ",L16)') supercomoving
  write (10,'("c_cosmo =       ",E16.8)') c_cosmo
  write (10,'("omega0 =        ",E16.8)') omega0
  write (10,'("lambda0 =       ",E16.8)') lambda0
  write (10,'("fb =            ",E16.8)') fb
  write (10,'("hubble0 =       ",E16.8)') hubble0
  write (10,'("box =           ",E16.8)') box
  write (10,'("zstart =        ",E16.8)') zstart
  write (10,'("zend =          ",E16.8)') zend
  write (10,'("chemistry =     ",L16)') chemistry
  write (10,'("noncie =        ",L16)') noncie
  write (10,'("cooling =       ",L16)') cooling
  write (10,'("background =    ",L16)') background  
  write (10,'("bg =            ",I16)') bg
  write (10,'("j0 =            ",E16.8)') j0
  write (10,'("c_chem =        ",E16.8)') c_chem
  write (10,'("substep =       ",E16.8)') substep
  write (10,'("rhoH =          ",E16.8)') rhoH
  write (10,'("rhoHe =         ",E16.8)') rhoHe
  write (10,'("cooldt =        ",L16)') cooldt
  write (10,'("c_cooldt =      ",E16.8)') c_cooldt
  write (10,'("tempfloor =     ",L16)') tempfloor
  write (10,'("temp_fl =       ",E16.8)') temp_fl
  write (10,'("jeansfloor =    ",L16)') jeansfloor
  write (10,'("c_jeans =       ",E16.8)') c_jeans
  write (10,'("conduction =    ",L16)') conduction
  write (10,'("c_cond =        ",E16.8)') c_cond
  write (10,'("simplekappa =   ",L16)') simplekappa
  write (10,'("dualenergy =    ",L16)') dualenergy
  write (10,'("de_thresh =     ",E16.8)') de_thresh

  close(10)

end subroutine writeparameter

subroutine rereadparameter()
  use global
  implicit none

  integer :: ierror
  character(16) :: tmp

  open (unit=10, file="data/parameter",form="formatted", &
       status="old", iostat=ierror)
  if (ierror .ne. 0) then
     print *,"ERROR file not found: data/parameter"
     call cleanstop()
  endif

  read (10,'(A16,L16)') tmp,iconly
  read (10,'(A16,L16)') tmp,hydro
  read (10,'(A16,L16)') tmp,secondorder
  read (10,'(A16,I16)') tmp,solver
  read (10,'(A16,I16)') tmp,limiter
  read (10,'(A16,E16.8)') tmp,cfl
  read (10,'(A16,E16.8)') tmp,g
  read (10,'(A16,L16)') tmp,atomic
  read (10,'(A16,I16)') tmp,nx
  read (10,'(A16,I16)') tmp,ny
  read (10,'(A16,I16)') tmp,nz
  read (10,'(A16,E16.8)') tmp,xmin
  read (10,'(A16,E16.8)') tmp,xmax
  read (10,'(A16,E16.8)') tmp,ymin
  read (10,'(A16,E16.8)') tmp,ymax
  read (10,'(A16,E16.8)') tmp,zmin
  read (10,'(A16,E16.8)') tmp,zmax
  read (10,'(A16,I16)') tmp,boundtype
  read (10,'(A16,L16)') tmp,cut
  read (10,'(A16,L16)') tmp,slice
  read (10,'(A16,L16)') tmp,bov
  read (10,'(A16,I16)') tmp,nsnap
  read (10,'(A16,L16)') tmp,gravity
  read (10,'(A16,E16.8)') tmp,tstart
  read (10,'(A16,E16.8)') tmp,tend
  read (10,'(A16,L16)') tmp,supercomoving
  read (10,'(A16,E16.8)') tmp,c_cosmo
  read (10,'(A16,E16.8)') tmp,omega0
  read (10,'(A16,E16.8)') tmp,lambda0
  read (10,'(A16,E16.8)') tmp,fb
  read (10,'(A16,E16.8)') tmp,hubble0
  read (10,'(A16,E16.8)') tmp,box
  read (10,'(A16,E16.8)') tmp,zstart
  read (10,'(A16,E16.8)') tmp,zend
  read (10,'(A16,L16)') tmp,chemistry
  read (10,'(A16,L16)') tmp,noncie
  read (10,'(A16,L16)') tmp,cooling
  read (10,'(A16,L16)') tmp,background
  read (10,'(A16,I16)') tmp,bg
  read (10,'(A16,E16.8)') tmp,j0
  read (10,'(A16,E16.8)') tmp,c_chem
  read (10,'(A16,E16.8)') tmp,substep
  read (10,'(A16,E16.8)') tmp,rhoH
  read (10,'(A16,E16.8)') tmp,rhoHe
  read (10,'(A16,L16)') tmp,cooldt
  read (10,'(A16,E16.8)') tmp,c_cooldt
  read (10,'(A16,L16)') tmp,tempfloor
  read (10,'(A16,E16.8)') tmp,temp_fl
  read (10,'(A16,L16)') tmp,jeansfloor
  read (10,'(A16,E16.8)') tmp,c_jeans
  read (10,'(A16,L16)') tmp,conduction
  read (10,'(A16,E16.8)') tmp,c_cond
  read (10,'(A16,L16)') tmp,conduction
  read (10,'(A16,L16)') tmp,dualenergy
  read (10,'(A16,E16.8)') tmp,de_thresh

  close(10)

end subroutine rereadparameter

