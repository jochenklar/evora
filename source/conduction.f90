!########################################################################!
!#                                                                      #!
!# conduction.f90                                                       #!
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

subroutine thermalconduction(q,w) 
  use global
  implicit none

  real(8), dimension(nq,lx:mx,ly:my,lz:mz) :: q
  real(8), dimension(nw,lx2:mx2,ly2:my2,lz2:mz2) :: w

  ! left and right temperature as seen by the interface
  real(8) :: templ,tempr,nel,ner 

  ! thermal fluxes as seen by the cell center
  real(8) :: jl,jr

  real(8), dimension(ns) :: n

  real(8) :: delta

  integer :: ix,iy,iz

  ! xsweep
  if (nx .ne. 1) then
     delta = dt * dxinv
     do iz=lz,mz
        do iy=ly,my
           call w2n(w(:,lx1,iy,iz),n)
           call press2temp(w(5,lx1,iy,iz),w(1,lx1,iy,iz),n,templ)
           nel = n(6)
           call w2n(w(:,lx,iy,iz),n)
           call press2temp(w(5,lx,iy,iz),w(1,lx,iy,iz),n,tempr)
           ner = n(6)
           call thermalflux(templ,tempr,nel,ner,jr,dxinv)
           do ix=lx,mx
              jl = jr
              templ = tempr
              nel   = ner
              call w2n(w(:,ix+1,iy,iz),n)
              call press2temp(w(5,ix+1,iy,iz),w(1,ix+1,iy,iz),n,tempr)
              ner = n(6)
              call thermalflux(templ,tempr,nel,ner,jr,dxinv)
              call thermalintegration(jl,jr,q(:,ix,iy,iz),delta)
           enddo
        enddo
     enddo
  endif

  ! ysweep
  if (ny .ne. 1) then
     delta = dt * dyinv
     do iz=lz,mz
        do ix=lx,mx
           call w2n(w(:,ix,ly1,iz),n)
           call press2temp(w(5,ix,ly1,iz),w(1,ix,ly1,iz),n,templ)
           nel = n(6)
           call w2n(w(:,ix,ly,iz),n)
           call press2temp(w(5,ix,ly,iz),w(1,ix,ly,iz),n,tempr)
           ner = n(6)
           call thermalflux(templ,tempr,nel,ner,jr,dyinv)
           do iy=ly,my
              jl = jr
              templ = tempr
              nel   = ner
              call w2n(w(:,ix,iy+1,iz),n)
              call press2temp(w(5,ix,iy+1,iz),w(1,ix,iy+1,iz),n,tempr)
              ner = n(6)
              call thermalflux(templ,tempr,nel,ner,jr,dyinv)
              call thermalintegration(jl,jr,q(:,ix,iy,iz),delta)
           enddo
        enddo
     enddo
  endif

  ! zsweep
  if (nz .ne. 1) then
     delta = dt * dzinv
     do iy=ly,my
        do ix=lx,mx
           call w2n(w(:,ix,iy,lz1),n)
           call press2temp(w(5,ix,iy,lz1),w(1,ix,iy,lz1),n,templ)
           nel = n(6)
           call w2n(w(:,ix,iy,lz),n)
           call press2temp(w(5,ix,iy,lz),w(1,ix,iy,lz),n,tempr)
           ner = n(6)
           call thermalflux(templ,tempr,nel,ner,jr,dzinv)
           do iz=lz,mz
              jl = jr
              templ = tempr
              nel   = ner
              call w2n(w(:,ix,iy,iz+1),n)
              call press2temp(w(5,ix,iy,iz+1),w(1,ix,iy,iz+1),n,tempr)
              ner = n(6)
              call thermalflux(templ,tempr,nel,ner,jr,dzinv)
              call thermalintegration(jl,jr,q(:,ix,iy,iz),delta)
           enddo
        enddo
     enddo
  endif

end subroutine thermalconduction

subroutine thermalflux(templ,tempr,nel,ner,j,deltainv)
   use global
  implicit none

  real(8), intent(in)  :: templ,tempr,nel,ner
  real(8), intent(out) :: j
  real(8), intent(in)  :: deltainv

  real(8) :: kappa,kappa_eff,meanfreepath,dtempdx

  if (simplekappa) then
     kappa_eff = 1.0 / g1
  else
     kappa = kappa0 * ((templ**2.5)+(tempr**2.5))
     meanfreepath = meanfreepath0 * ((templ / nel)+(tempr / ner))
     kappa_eff = kappa / (1.0 + meanfreepath * abs(dtempdx))
  endif
  
  dtempdx = (tempr - templ) * deltainv
  
  j = - kappa_eff * dtempdx

end subroutine thermalflux

subroutine thermalintegration(jl,jr,q,delta)
  use global
  implicit none
  
  real(8), intent(in) :: jl,jr
  real(8), intent(inout), dimension(nq) :: q
  real(8), intent(in) :: delta

  real(8) :: de

  de = delta * (jl - jr)

  q(5) = q(5) + de
  q(6) = q(6) + g1 * de / (q(1)**g1)

end subroutine thermalintegration
