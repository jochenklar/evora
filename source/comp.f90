!########################################################################!
!#                                                                      #!
!# comp.f90                                                             #!
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

subroutine q2w(q,w)
  use global
  implicit none
  
  real(8), dimension(nq), intent(in)  :: q
  real(8), dimension(nw), intent(out) :: w

  real(8) :: densinv
  integer :: i

  densinv = 1.0 / q(1)
  w(1) = q(1)
  w(2:4) = q(2:4) * densinv
  call pressure(q(:),w(5))
  do i=6,nw
     w(i) = q(i+1) * densinv
  enddo
end subroutine q2w

subroutine pressure(q,press)
  ! computes the p from given conservative variables
  ! using energy density or entropy (dual energy formalism)
  use global
  implicit none

  real(8), intent(inout), dimension(nq) :: q
  real(8), intent(out)                  :: press

  real(8) :: m2,kin,ratio

  m2 = (q(2)**2) + (q(3)**2) + (q(4)**2)
  kin = 0.5 * m2 / q(1)
  ratio = kin / q(5) 

  if (dualenergy) then
     if (ratio .gt. de_thresh1) then
        ! S - system
        press = q(6) * (q(1)**g1)
        q(5) = kin + press * g1inv
     else
        ! E - System
        press = g1 * (q(5) - kin)
        q(6) = press / (q(1)**g1)
     endif
  else
     press = g1 * (q(5) - kin)
  endif

  if (press .le. 0.0) then
     print *,'ERROR negative pressure occured while computing'
     print *,'primitive variables!'
     print *,'time:', t
     print *,'redshift:', redshift
     print *,de_thresh1
     print *,kin,q(5),ratio
     print *,press,q(6)
     stop
  endif
  
end subroutine pressure

subroutine press2temp(press,rho,n,temp)
  ! computes the T from given p,rho,n
  use global
  implicit none

  real(8), intent(in)                :: press,rho
  real(8), intent(in), dimension(ns) :: n
  real(8), intent(out)               :: temp

  if (chemistry) then
     temp = press / sum(n) * temp0
  else
     temp = press / rho * temp0
  endif

end subroutine press2temp

subroutine temp2press(temp,rho,n,press)
  ! computes the p from given T,rho,n
  use global
  implicit none

  real(8), intent(in)                :: temp,rho
  real(8), intent(in), dimension(ns) :: n
  real(8), intent(out)               :: press

  if (chemistry) then
     press = temp * sum(n) * temp0inv
  else
     press = temp * rho * temp0inv
  endif

end subroutine temp2press

subroutine w2n(w,n)
  ! computes n from w
  use global
  implicit none

  real(8), intent(in),  dimension(nw) :: w
  real(8), intent(out), dimension(ns) :: n

  integer :: i

  do i=6,nw
     n(i-5) = w(i) * w(1) 
  enddo

end subroutine w2n

subroutine jeanslength(rho,press,jeans)
  ! computes the ratio between cell size and the local jeans lenght
  ! in box units
  use global
  implicit none

  real(8), intent(in)  :: rho,press
  real(8), intent(out) :: jeans

  jeans = jeans0 * sqrt(press) / rho 

end subroutine jeanslength

subroutine coordinates(ix,iy,iz,x,y,z)
  ! computes the comooving coordinates of a cell in box units
  use global
  implicit none

  integer :: ix,iy,iz
  real(8) :: x,y,z

  x = x0 + dble(ix-1) * dx;
  y = y0 + dble(iy-1) * dy;
  z = z0 + dble(iz-1) * dz;

end subroutine coordinates

subroutine statistic(a)
  ! computes the mean and the variance of an array a(i,j,k)
  use global
  implicit none

  real(8), dimension(lx:mx,ly:my,lz:mz) :: a
  real(8) :: mean,sigma
  integer :: ix,iy,iz

  ! get the mean value
  mean = 0
  do iz=lz,mz
     do iy=ly,my
        do ix=lx,mx
           mean = mean + a(ix,iy,iz)
        enddo
     enddo
  enddo  
  ! mpi stuff here 
  mean = mean/n3

  ! get the variance
  sigma = 0
  do iz=lz,mz
     do iy=ly,my
        do ix=lx,mx
           sigma = sigma + abs(mean - a(ix,iy,iz))
        enddo
     enddo
  enddo    
  ! mpi stuff here
  sigma = sigma/n3
  
  write(*,'("mean = ",F8.3," variance = ",F8.3)') mean,sigma

end subroutine statistic

subroutine compstate(q,state)
  ! computes the state vector containing physical and dimensionless
  ! quantities mainly for ascii-output
  use global
  implicit none

  real(8), dimension(nq) :: q
  real(8), dimension(nstate) :: state

  real(8), dimension(nchem) :: chemrate
  real(8), dimension(ncool) :: coolrate

  real(8) :: press,temp
  integer :: i

  call pressure(q,press)
  call press2temp(press,q(1),q(7:nq),temp)

  state(1)    = q(1) / fb
  state(2:6)  = q(2:6)
  state(7)    = state(1) * density0
  state(8:10) = q(2:4) / state(1) * velocity0
  state(11)   = press * press0
  state(12)   = temp

  if (chemistry) then
     do i=7,nq
        state(i+6) = q(i) / state(1)
        state(i+12) = q(i) * number0   
     enddo
     if (cooling) then
        nH  = chiH  * q(1)
        nHe = chiHe * q(1)
        call comprates(temp,chemrate,coolrate)
        call cieabundance(q(7:nq),chemrate)
        call cool(q(7:nq),coolrate,state(nstate-1),state(nstate))
     endif
  endif
end subroutine compstate

subroutine kvalue(i,n,nyq,k)
  implicit none

  integer, intent(in)  :: i,n,nyq
  real(8), intent(out) :: k

  if (i <= nyq) then 
     k = dble(i-1)
  else 
     k = - dble(n - (i-1))
  endif
end subroutine kvalue
