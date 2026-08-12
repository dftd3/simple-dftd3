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

module dftd3_ncoord
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   implicit none
   private

   public :: get_coordination_number, add_coordination_number_derivs
   public :: add_coordination_number_hessian

   !> Steepness of counting function
   real(wp), parameter :: default_kcn = 16.0_wp

contains


!> Wrapper for geometric fractional coordination number
!> with standard exponential counting function.
subroutine get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out), optional :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out), optional :: dcndL(:, :, :)

   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: error

   call new_ncoord(ncoord, mol, cn_count%exp, error, &
      & kcn=default_kcn, cutoff=cutoff, rcov=rcov)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)

end subroutine get_coordination_number


subroutine add_coordination_number_derivs(mol, trans, cutoff, rcov, dEdcn, gradient, sigma)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Derivative of expression with respect to the coordination number
   real(wp), intent(in) :: dEdcn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates
   real(wp), intent(inout) :: gradient(:, :)

   !> Derivative of the CN with respect to strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: error

   call new_ncoord(ncoord, mol, cn_count%exp, error, &
      & kcn=default_kcn, cutoff=cutoff, rcov=rcov)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call ncoord%add_coordination_number_derivs(mol, trans, dEdcn, gradient, sigma)

end subroutine add_coordination_number_derivs


!> Add the second derivative of the exponential coordination number contracted
!> with the derivative of the energy w.r.t. the coordination number.
!>
!> The counting function is reproduced here rather than taken from mctc-lib,
!> which only provides derivatives up to first order.
subroutine add_coordination_number_hessian(mol, trans, cutoff, rcov, dEdcn, hessian)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Derivative of expression with respect to the coordination number
   real(wp), intent(in) :: dEdcn(:)

   !> Second derivative of the energy w.r.t. the Cartesian coordinates
   real(wp), intent(inout) :: hessian(:, :)

   integer :: iat, jat, izp, jzp, itr, ic, jc, ii, jj
   real(wp) :: vec(3), r2, r1, cutoff2, rc, expterm, cf, dcf, d2cf
   real(wp) :: dEdcnij, block(3, 3), tmp

   cutoff2 = cutoff*cutoff

   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         rc = rcov(izp) + rcov(jzp)
         dEdcnij = dEdcn(iat) + dEdcn(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            expterm = exp(-default_kcn*(rc/r1 - 1.0_wp))
            cf = 1.0_wp/(1.0_wp + expterm)
            ! first and second derivative of the counting function w.r.t. r
            tmp = cf * (1.0_wp - cf)
            dcf = -tmp * default_kcn * rc / r2
            d2cf = tmp * (1.0_wp - 2.0_wp*cf) * (default_kcn*rc)**2 / (r2*r2) &
               & + 2.0_wp * tmp * default_kcn * rc / (r2*r1)

            do ic = 1, 3
               do jc = 1, 3
                  block(ic, jc) = dEdcnij * ( &
                     & d2cf * vec(ic) * vec(jc) / r2 &
                     & - dcf * vec(ic) * vec(jc) / (r2*r1))
               end do
               block(ic, ic) = block(ic, ic) + dEdcnij * dcf / r1
            end do

            do ic = 1, 3
               ii = 3*(iat - 1) + ic
               do jc = 1, 3
                  jj = 3*(jat - 1) + jc
                  hessian(ii, 3*(iat - 1) + jc) = hessian(ii, 3*(iat - 1) + jc) + block(ic, jc)
                  hessian(3*(jat - 1) + ic, jj) = hessian(3*(jat - 1) + ic, jj) + block(ic, jc)
                  hessian(ii, jj) = hessian(ii, jj) - block(ic, jc)
                  hessian(3*(jat - 1) + ic, 3*(iat - 1) + jc) = &
                     & hessian(3*(jat - 1) + ic, 3*(iat - 1) + jc) - block(ic, jc)
               end do
            end do
         end do
      end do
   end do

end subroutine add_coordination_number_hessian


end module dftd3_ncoord
