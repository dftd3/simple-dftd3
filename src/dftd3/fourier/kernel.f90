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
   public :: fourier_kernel, new_fourier_kernel
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


   !> Wave-number-independent coefficients of the Fourier transform
   type :: fourier_kernel
      private
      real(wp) :: prefactor = 0.0_wp, radius = 0.0_wp, zero = 0.0_wp
      real(wp), allocatable :: angle(:), cosine(:), sine(:), phase(:)
   end type fourier_kernel

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


!> Prepare a term once before evaluating it at many wave numbers
pure elemental function new_fourier_kernel(term) result(kernel)

   type(fourier_term), intent(in) :: term
   type(fourier_kernel) :: kernel
   integer :: ipole, nexp

   if (abs(term%prefactor) <= 0.0_wp) return

   nexp = term%m - term%alpha + 2
   kernel%prefactor = 4.0_wp*pi*pi / term%alpha * term%prefactor * term%radius**nexp
   kernel%radius = term%radius
   kernel%zero = kernel%prefactor * term%radius / sin(pi*(term%m + 3)/real(term%alpha, wp))
   kernel%angle = [(pi * (2*ipole + 1) / real(term%alpha, wp), ipole=0, term%alpha/2-1)]
   kernel%cosine = cos(kernel%angle)
   kernel%sine = sin(kernel%angle)
   kernel%phase = 0.5_wp*pi + nexp*kernel%angle

end function new_fourier_kernel


!> Fourier transform and its derivative w.r.t. the wave number
pure subroutine get_fourier_transform(kernel, kval, phi, dphi, grad)

   !> Prepared coefficients of the damped pair potential
   type(fourier_kernel), intent(in) :: kernel

   !> Wave number
   real(wp), intent(in) :: kval

   !> Fourier transform of the term
   real(wp), intent(out) :: phi

   !> Derivative of the Fourier transform w.r.t. the wave number
   real(wp), intent(out) :: dphi

   !> Whether the derivative is required
   logical, intent(in) :: grad

   integer :: ipole
   real(wp) :: arg, damp, acc, dacc, krad

   phi = 0.0_wp
   dphi = 0.0_wp
   if (abs(kernel%prefactor) <= 0.0_wp) return

   if (kval <= 0.0_wp) then
      phi = kernel%zero
      return
   end if

   krad = kval * kernel%radius
   acc = 0.0_wp
   dacc = 0.0_wp
   do ipole = 1, size(kernel%angle)
      arg = kernel%phase(ipole) + krad*kernel%cosine(ipole)
      damp = exp(-krad*kernel%sine(ipole))
      acc = acc + sin(arg) * damp
      if (grad) dacc = dacc + cos(arg + kernel%angle(ipole)) * damp
   end do

   phi = kernel%prefactor * acc / kval
   if (grad) dphi = (kernel%prefactor * kernel%radius * dacc - phi) / kval

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
