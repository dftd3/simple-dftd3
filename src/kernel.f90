! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

!> Reciprocal space representation of the damped dispersion potential.
!>
!> All damping functions of the rational family reduce the pair potential to a sum
!> of terms of the shape
!>
!>    phi(r) = prefactor * r**m / (r**alpha + radius**alpha)
!>
!> with even m and even alpha. Since the poles of the denominator are known, the
!> three-dimensional Fourier transform
!>
!>    phi_hat(k) = 4*pi/k * integral_0^infinity r*sin(k*r)*phi(r) dr
!>
!> follows from contour integration over the alpha/2 poles in the upper half plane
!> and decays exponentially with exp(-k*radius*sin(pi/alpha)).
module dftd3_fourier_kernel
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: fourier_term, max_fourier_terms
   public :: is_supported_term, get_fourier_transform, get_potential_zero
   public :: get_reciprocal_cutoff

   !> Maximum number of terms needed to represent a damped pair potential
   integer, parameter :: max_fourier_terms = 2


   !> Single term of the damped pair potential
   type :: fourier_term

      !> Scaling of the term
      real(wp) :: prefactor = 0.0_wp

      !> Power of the numerator
      integer :: m = 0

      !> Power of the denominator
      integer :: alpha = 0

      !> Damping radius entering the denominator
      real(wp) :: radius = 0.0_wp

   end type fourier_term


contains


!> Check whether the closed form of the Fourier transform is applicable
pure function is_supported_term(term) result(supported)

   !> Term of the damped pair potential
   type(fourier_term), intent(in) :: term

   !> Whether the term can be transformed analytically
   logical :: supported

   supported = term%alpha > 0 .and. term%m >= 0 &
      & .and. mod(term%alpha, 2) == 0 .and. mod(term%m, 2) == 0 &
      & .and. term%m + 3 < term%alpha .and. term%radius > 0.0_wp

end function is_supported_term


!> Fourier transform of a single term and its derivative w.r.t. the wave number
pure subroutine get_fourier_transform(term, kval, phi, dphi)

   !> Term of the damped pair potential
   type(fourier_term), intent(in) :: term

   !> Wave number
   real(wp), intent(in) :: kval

   !> Fourier transform of the term
   real(wp), intent(out) :: phi

   !> Derivative of the Fourier transform w.r.t. the wave number
   real(wp), intent(out) :: dphi

   integer :: ipole, npole, nexp
   real(wp) :: pre, theta, ct, st, arg, damp, acc, dacc, krad

   phi = 0.0_wp
   dphi = 0.0_wp
   if (abs(term%prefactor) <= 0.0_wp) return

   nexp = term%m - term%alpha + 2
   pre = 4.0_wp*pi*pi / term%alpha * term%prefactor * term%radius**nexp

   if (kval <= 0.0_wp) then
      phi = pre * term%radius / sin(pi*(term%m + 3)/real(term%alpha, wp))
      return
   end if

   npole = term%alpha / 2
   krad = kval * term%radius
   acc = 0.0_wp
   dacc = 0.0_wp
   do ipole = 0, npole - 1
      theta = pi * (2*ipole + 1) / real(term%alpha, wp)
      ct = cos(theta)
      st = sin(theta)
      arg = 0.5_wp*pi + nexp*theta + krad*ct
      damp = exp(-krad*st)
      acc = acc + sin(arg) * damp
      dacc = dacc + cos(arg + theta) * damp
   end do

   phi = pre * acc / kval
   dphi = (pre * term%radius * dacc - phi) / kval

end subroutine get_fourier_transform


!> Value of the damped pair potential at vanishing distance
pure function get_potential_zero(term) result(phi)

   !> Term of the damped pair potential
   type(fourier_term), intent(in) :: term

   !> Value of the term at the origin
   real(wp) :: phi

   if (term%m > 0) then
      phi = 0.0_wp
   else
      phi = term%prefactor / term%radius**term%alpha
   end if

end function get_potential_zero


!> Wave number beyond which the term is damped below the requested tolerance
pure function get_reciprocal_cutoff(term, tolerance) result(kcut)

   !> Term of the damped pair potential
   type(fourier_term), intent(in) :: term

   !> Requested relative accuracy
   real(wp), intent(in) :: tolerance

   !> Reciprocal space cutoff
   real(wp) :: kcut

   kcut = -log(tolerance) / (term%radius * sin(pi/real(term%alpha, wp)))

end function get_reciprocal_cutoff


end module dftd3_fourier_kernel
