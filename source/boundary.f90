!########################################################################!
!#                                                                      #!
!# boundary.f90                                                         #!
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

subroutine bound(w)
  use global
  implicit none

  real(8), dimension(nw,lx2:mx2,ly2:my2,lz2:mz2) :: w

  ! xboundary
  if (nx .ne. 1) then
     if (lx .eq. 1) then ! lower xboundary
        if (boundtype .eq. 1 .or. boundtype .eq. 3) then ! zero gradient
           w(:,lx2,:,:) = w(:,lx0,:,:)
           w(:,lx1,:,:) = w(:,lx,:,:)
        else if (boundtype .eq. 2) then ! periodic
           w(:,lx2,:,:) = w(:,mx0,:,:)
           w(:,lx1,:,:) = w(:,mx,:,:)
        else if (boundtype .eq. 3) then ! reflecting
           w(1   ,lx2,:,:) =   w(1   ,lx0,:,:)
           w(2:4 ,lx2,:,:) = - w(2:4 ,lx0,:,:)
           w(5:nw,lx2,:,:) =   w(5:nw,lx0,:,:)
           w(1   ,lx1,:,:) =   w(1   ,lx ,:,:)
           w(2:4 ,lx1,:,:) = - w(2:4 ,lx ,:,:)
           w(5:nw,lx1,:,:) =   w(5:nw,lx ,:,:)
        endif
     endif
     if (mx .eq. nx) then ! upper xboundary
        if (boundtype .eq. 1) then ! zero gradient
           w(:,mx2,:,:) = w(:,mx0,:,:)
           w(:,mx1,:,:) = w(:,mx,:,:)
        else if (boundtype .eq. 2) then ! periodic
           w(:,mx2,:,:) = w(:,lx0,:,:)
           w(:,mx1,:,:) = w(:,lx,:,:)
        else if (boundtype .eq. 3) then ! reflecting
           w(1   ,mx1,:,:) =   w(1   ,mx ,:,:)
           w(2:4 ,mx1,:,:) = - w(2:4  ,mx ,:,:)
           w(5:nw,mx1,:,:) =   w(5:nw,mx ,:,:)
           w(1   ,mx2,:,:) =   w(1   ,mx0,:,:)
           w(2:4 ,mx2,:,:) = - w(2:4 ,mx0,:,:)
           w(5:nw,mx2,:,:) =   w(5:nw,mx0,:,:)
        endif
     endif
  endif

  ! yboundary
  if (ny .ne. 1) then
     if (ly .eq. 1) then ! lower yboundary
        if (boundtype .eq. 1 .or. boundtype .eq. 3) then ! zero gradient
           w(:,:,ly2,:) = w(:,:,ly0,:)
           w(:,:,ly1,:) = w(:,:,ly,:) 
        else if (boundtype .eq. 2) then ! periodic
           w(:,:,ly2,:) = w(:,:,my0,:)
           w(:,:,ly1,:) = w(:,:,my,:) 
        endif
     endif
     if (my .eq. ny) then ! upper yboundary
        if (boundtype .eq. 1 .or. boundtype .eq. 3) then ! zero gradient
           w(:,:,my2,:) = w(:,:,my0,:)
           w(:,:,my1,:) = w(:,:,my,:)  
        else if (boundtype .eq. 2) then ! periodic
           w(:,:,my2,:) = w(:,:,ly0,:)
           w(:,:,my1,:) = w(:,:,ly,:) 
        endif
     endif
  endif

  ! zboundary
  if (nz .ne. 1) then 
     if (lz .eq. 1) then ! lower zboundary
        if (boundtype .eq. 1 .or. boundtype .eq. 3) then ! zero gradient
           ! compute free-flow boundary conditions
           w(:,:,:,lz2) = w(:,:,:,lz0)
           w(:,:,:,lz1) = w(:,:,:,lz)
        else if (boundtype .eq. 2 .and. mpisize .eq. 1) then ! periodic
           w(:,:,:,lz2) = w(:,:,:,mz0)
           w(:,:,:,lz1) = w(:,:,:,mz) 
        endif
     endif
     if (mz .eq. nz) then ! upper zboundary
        if (boundtype .eq. 1 .or. boundtype .eq. 3) then
           w(:,:,:,mz2) = w(:,:,:,mz0)
           w(:,:,:,mz1) = w(:,:,:,mz)
        else if (boundtype .eq. 2 .and. mpisize .eq. 1) then ! periodic
           w(:,:,:,mz2) = w(:,:,:,lz0)
           w(:,:,:,mz1) = w(:,:,:,lz)  
        endif
     endif
  endif

end subroutine bound

