module global
  implicit none 

  ! number of variables
#ifdef NONCIE
  integer, parameter :: nq=12, nf=12, nw=11, nstate=26
#elif CIE
  integer, parameter :: nq=12, nf= 6, nw= 5, nstate=26
#else
  integer, parameter :: nq= 6, nf= 6, nw= 5, nstate=12
#endif

  ! some (almost) global variables, to be read from parameter.txt
  ! or to be calculated by subroutine readparameter

  ! general 
  logical :: iconly 

  ! hydro
  logical :: secondorder,atomic,hydro
  integer :: solver,limiter
  real(8) :: g,cfl,g1,g1inv,g2,gg1,g35

  ! gravity
  logical :: gravity
  real(8) :: grav
  integer :: nyqx,nyqy,nyqz
  real(8) :: h2,kfx,kfy,kfz,kfx2,kfy2,kfz2

  ! cosmology
  logical :: supercomoving
  real(8) :: c_cosmo,omega0,lambda0,fb,box,hubble0

  ! dual energy
  logical :: dualenergy
  real(8) :: de_thresh,de_thresh1

  ! chemistry
  logical :: chemistry,noncie,background
  logical :: cooling,cooldt,tempfloor,jeansfloor
  real(8) :: c_chem,c_cooldt,rhoH,rhoHe,substep,temp_fl,j0,c_jeans
  integer :: bg

  ! thermal conduction
  logical :: conduction,simplekappa
  real(8) :: c_cond

  ! grid
  integer :: nx,ny,nz
  real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
  integer :: boundtype
  integer :: n3,nx1,ny1,nz1,nx2,ny2,nz2
  integer :: lx,ly,lz,lx1,ly1,lz1,lx2,ly2,lz2,lx0,ly0,lz0
  integer :: mx,my,mz,mx1,my1,mz1,mx2,my2,mz2,mx0,my0,mz0
  real(8) :: dx,dy,dz,dx2,dy2,dz2,dxinv,dyinv,dzinv,dx2inv,dy2inv,dz2inv
  real(8) :: x0,y0,z0,n3inv,volume
  logical :: onedimension,twodimension

  ! time
  real(8) :: tstart,tend,t,dt,tsnap
  real(8) :: zstart,zend,expand,cosmo,oldcosmo,cosmodot,redshift

  ! output
  logical :: cut,slice,bov
  integer :: nsnap

  ! units
  real(8) :: length0, time0, density0
  real(8) :: velocity0,temp0,temp0inv,press0,number0,chem0,photo0,cool0,heat0
  real(8) :: jeans0,jeansfloor0,kappa0,meanfreepath0

  ! physical constants and mathematical stuff
  real(8), parameter :: pi =         3.1415927
  real(8), parameter :: lightspeed = 2.9979250d+10
  real(8), parameter :: gravconst =  6.67259d-8
  real(8), parameter :: fourpig =    8.38553911d-7
  real(8), parameter :: rhoc =       1.88d-29
  real(8), parameter :: kb =         1.3806200d-16
  real(8), parameter :: me =         9.109d-28
  real(8), parameter :: charge =     4.8032068d-10
  real(8), parameter :: avogadro =   6.0221418d23
  real(8), parameter :: atomicmass = 1.660538782d-24
  real(8), parameter :: solarmass =  1.988435d33
  real(8), parameter :: mH =         1.00794
  real(8), parameter :: mHe =        4.002602 
  real(8), parameter :: megapc =     3.08568025d24
  real(8), parameter :: kilometer =  1d5
  real(8), parameter :: year =       31557600
  real(8), parameter :: gff =        1.5
  real(8), parameter :: sigma_lya =  4.5d-18
  real(8), parameter :: sigma_th =   6.65d-25


  real(8) :: third,twothird,fourthird
 
  ! chemical stuff
  integer, parameter :: ns=6,ns1=5
  integer, parameter :: nchem=9,ncool=13
  real(8), parameter :: nfloor = 1.0d-40,ratefloor = 1.0d-40

  real(8), parameter :: cie_tol = 1d-10
  integer, parameter :: cie_nit = 1000 
  real(8), parameter :: noncie_tol = 1d-10
  integer, parameter :: noncie_nit = 10
  real(8), parameter :: cool_tol = 1d-10
  integer, parameter :: cool_nit = 1000

  real(8) :: nH,nHe,chiH,chiHe,uv
  real(8), dimension(ncool,2) :: coef

  ! snapshots
  integer, parameter :: ndump = 100
  real(8), dimension(ndump)  :: dump_redshift

  ! fftw 
  integer(8) :: there, back ! plans
  integer :: local_nx,local_x_start
  integer :: local_ny_after_transpose,local_y_start_after_transpose
  integer :: total_local_size

  ! mpi stuff
  logical :: parallel
  logical :: master,notmaster,last,notlast,even,noteven
  integer :: mpirank,mpisize,mpierror,mpistatus

  ! parameter file 
  character(len=30) :: parameterfile

  ! fftw flags
  integer, parameter :: FFTW_FORWARD=-1, FFTW_BACKWARD=1
  integer, parameter :: FFTW_REAL_TO_COMPLEX=-1, FFTW_COMPLEX_TO_REAL=1
  integer, parameter :: FFTW_ESTIMATE=0, FFTW_MEASURE=1
  integer, parameter :: FFTW_OUT_OF_PLACE=0, FFTW_IN_PLACE=8
  integer, parameter :: FFTW_USE_WISDOM=16, FFTW_THREADSAFE=128
  integer, parameter :: FFTW_NORMAL_ORDER=0,FFTW_TRANSPOSED_ORDER=1
  integer, parameter :: FFTW_SCRAMBLED_INPUT=8192,FFTW_SCRAMBLED_OUTPUT=16384

end module global
