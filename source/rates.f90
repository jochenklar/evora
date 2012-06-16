!########################################################################!
!#                                                                      #!
!# rates.f90                                                            #!
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

subroutine compcoef()
  use global
  implicit none

  ! collisional ionisation
  ! HI
  coef(1,1)  = chem0 * 5.85d-11
  coef(1,2)  = cool0 * 1.27d-21
  ! HeI
  coef(2,1)  = chem0 * 2.38d-11
  coef(2,2)  = cool0 * 9.38d-22
  ! HeII
  coef(3,1)  = chem0 * 5.68d-12
  coef(3,2)  = cool0 * 4.95d-22

  ! recombination
  ! HII
  coef(4,1)  = chem0 * 8.4d-11
  coef(4,2)  = cool0 * 8.7d-27
  ! HeI
  coef(5,1)  = chem0 * 1.5d-10
  coef(5,2)  = cool0 * 1.55d-26
  ! HeII 
  coef(6,1)  = chem0 * 3.36d-10
  coef(6,2)  = cool0 * 3.48d-26

  ! dielectric recombination
  coef(7,1)  = chem0 * 1.9d-3
  coef(7,2)  = cool0 * 1.24d-13

  ! collisional exitation cooling
  ! HI
  coef(8,2)  = cool0 * 7.5d-19
  ! HeII
  coef(9,2)  = cool0 * 5.54d-17
  
  ! bremsstrahlung
  coef(10,2) = cool0 * 1.42d-27 * gff

  ! photo ionisation
  ! HI
  coef(11,1) = photo0 * 2.54d8 * uv
  coef(11,2) = heat0  * 2.54d8 * 7.75d-12  * uv

  ! HeI
  coef(12,1) = photo0 * 2.49d8  * uv
  coef(12,2) = heat0  * 2.49d8 * 2.19d-11  * uv

  ! HeII
  coef(13,1) = photo0 * 1.60d7  * uv
  coef(13,2) = heat0  * 1.60d7 * 3.10d-11  * uv

end subroutine compcoef

subroutine comprates(temp,chemrate,coolrate)
  use global
  implicit none

  real(8) :: temp
  real(8), dimension(nchem) :: chemrate
  real(8), dimension(ncool) :: coolrate
  real(8) :: temproot,tempinv,temp15
  real(8) :: fec,fc1,e1,e2,e3,fr,fr1,fr2,edr

  ! some things we dont want to compute double
  temproot = sqrt(temp)
  tempinv = 1.0 / temp
  temp15 = temp**(-1.5)
 
  fec = 1.0 / (1.0 + sqrt(temp * 1d-5))
  fc1 = temproot * fec
 
  e1 = exp(-157809.1 * tempinv)
  e2 = exp(-286335.4 * tempinv)
  e3 = exp(-631515.0 * tempinv)

  fr = ((temp * 1d-3)**(-0.2)) / (1.0 + ((temp * 1d-6)**0.7))
  fr1 = fr / temproot
  fr2 = fr * temproot

  edr = exp((-470000.0) * tempinv) * (1.0 + 0.3 * (-94000.0) * tempinv)

  ! collisional ionisation 
  chemrate(1) = coef(1,1) * fc1 * e1                          ! HI
  chemrate(2) = coef(2,1) * fc1 * e2                          ! HeI
  chemrate(3) = coef(3,1) * fc1 * e3                          ! HeII
  ! recombination
  chemrate(4) = coef(4,1) * fr1                               ! HII
  chemrate(5) =  coef(5,1) * (temp**(-0.6353)) &                ! HeII
       + coef(7,1) * temp15 * edr
  chemrate(6) = coef(6,1) * fr1                               ! HeIII

  if (background) then
     chemrate(7) = coef(11,1)
     chemrate(8) = coef(12,1)
     chemrate(9) = coef(13,1)
  else
     chemrate(7:9) = 0.0
  endif

  if (cooling) then
     ! collisional ionisation 
     coolrate(1) = coef(1,2) * fc1 * e1                       ! HI
     coolrate(2) = coef(2,2) * fc1 * e2                       ! HeI
     coolrate(3) = coef(3,2) * fc1 * e3                       ! HeII
     ! recombination
     coolrate(4) = coef(4,2) * fr2                            ! HII
     coolrate(5) = coef(5,2) * (temp**0.3647)                 ! HeII
     coolrate(6) = coef(6,2) * fr2                            ! HeIII
     ! dielectric recombination
     coolrate(7) = coef(7,2) * temp15 * edr                   
     ! collisional exitation cooling
     coolrate(8) = coef(8,2) * fec * exp(-118348.0 * tempinv) ! HI
     coolrate(9) = coef(9,2) * (temp**(-0.397)) * fec &         ! HeII
          * exp(-473638.0 * tempinv)
     ! bremsstrahlung
     coolrate(10) = coef(10,2) * temproot

     if (background) then
        coolrate(11) = coef(11,2)
        coolrate(12) = coef(12,2)
        coolrate(13) = coef(13,2)
     else
        coolrate(11) = 0.0
        coolrate(12) = 0.0
        coolrate(13) = 0.0
     endif
  endif
   
end subroutine comprates

subroutine cool(n,coolrate,lambda,gamma)
  use global
  implicit none

  real(8), dimension(ns),    intent(in) :: n
  real(8), dimension(ncool), intent(in) :: coolrate
  real(8), intent(out)                  :: lambda,gamma

  ! compute proper cooling rate
  if (cooling) then
     lambda = coolrate(1)  * n(1) &
          + coolrate(2) * n(3) &
          + coolrate(3)  * n(4) &
          + coolrate(4)  * n(2) &
          + coolrate(5)  * n(4) &
          + coolrate(6)  * n(5) &
          + coolrate(7)  * n(4) &
          + coolrate(8)  * n(1) &
          + coolrate(9)  * n(4) &
          + coolrate(10) * (n(2) + n(4) + 4.0 * n(5))

     lambda = lambda * n(6)
     
     if (background) then
        gamma = n(1) * coolrate(11) &
             + n(3) * coolrate(12) &
             + n(4) * coolrate(13)
     else
        gamma = 0.0
     endif
  endif

end subroutine cool

