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

!> Translation of accuracy targets into summation parameters.
module dftd3_accuracy
   use dftd3_cutoff, only : realspace_cutoff, get_lattice_points
   use dftd3_model, only : d3_model
   use dftd3_ncoord, only : get_coordination_number
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3
   implicit none
   private

   public :: get_realspace_cutoff

   !> Margin on the estimated truncation error to cover the neglected C8 tail
   real(wp), parameter :: cutoff_safety = 1.25_wp

   !> Shortest cutoff for which the asymptotic tail estimate is reliable
   real(wp), parameter :: cutoff_minimum = 25.0_wp


contains


!> Realspace cutoffs reproducing the two-body dispersion energy within a
!> requested accuracy per atom.
!>
!> Beyond the damping region the pair sum is dominated by the bare C6 tail.
!> Integrating it over the d periodic dimensions gives the truncation error
!>
!>    dE/atom = rho * <C6> * S(d) / (2 * (6 - d) * R**(6 - d))
!>
!> with the number density rho of the periodic cell, the mean pair coefficient
!> <C6> and the surface S(d) of the d-dimensional unit sphere. Inverting this
!> relation yields the cutoff. The estimate assumes a scaling factor s6 of one
!> and neglects the faster decaying C8 tail, both covered by a safety margin.
!>
!> Only the two-body cutoff is derived, the coordination number cutoff is not an
!> accuracy parameter because the counting function has a non-vanishing limit,
!> and the three-body cutoff is left at its default.
function get_realspace_cutoff(mol, disp, accuracy) result(cutoff)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Requested accuracy of the two-body dispersion energy in Hartree per atom
   real(wp), intent(in) :: accuracy

   !> Realspace cutoffs
   type(realspace_cutoff) :: cutoff

   integer :: mref, ndim, iat, izp, isp, jsp, iref, jref, jat
   real(wp) :: measure, c6mean, amp, r2
   real(wp), allocatable :: cn(:), gwvec(:, :), wsum(:, :), lattr(:, :)

   cutoff = realspace_cutoff()
   if (accuracy <= 0.0_wp) return

   mref = maxval(disp%ref)
   allocate(cn(mol%nat), gwvec(mref, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, cn)
   call disp%weight_references(mol, cn, gwvec)

   ! contract the reference coefficients with the accumulated weights per species
   allocate(wsum(mref, mol%nid), source=0.0_wp)
   do iat = 1, mol%nat
      if (disp%ghost(iat)) cycle
      izp = mol%id(iat)
      wsum(:, izp) = wsum(:, izp) + gwvec(:, iat)
   end do
   c6mean = 0.0_wp
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         do iref = 1, disp%ref(isp)
            do jref = 1, disp%ref(jsp)
               c6mean = c6mean + disp%c6(iref, jref, isp, jsp) &
                  & * wsum(iref, isp) * wsum(jref, jsp)
            end do
         end do
      end do
   end do
   c6mean = c6mean / mol%nat**2

   ndim = count(mol%periodic)
   measure = get_periodic_measure(mol, ndim)

   select case(ndim)
   case(3)
      amp = 2.0_wp*pi/3.0_wp * mol%nat/measure * c6mean
      cutoff%disp2 = (cutoff_safety * amp / accuracy)**(1.0_wp/3.0_wp)
   case(2)
      amp = 0.25_wp*pi * mol%nat/measure * c6mean
      cutoff%disp2 = (cutoff_safety * amp / accuracy)**(0.25_wp)
   case(1)
      amp = 0.2_wp * mol%nat/measure * c6mean
      cutoff%disp2 = (cutoff_safety * amp / accuracy)**(0.2_wp)
   case default
      ! a finite system is summed exactly once all pairs are covered
      r2 = 0.0_wp
      do iat = 1, mol%nat
         do jat = 1, iat - 1
            r2 = max(r2, sum((mol%xyz(:, iat) - mol%xyz(:, jat))**2))
         end do
      end do
      cutoff%disp2 = sqrt(r2)
   end select

   cutoff%disp2 = max(cutoff%disp2, cutoff_minimum)
   cutoff%disp3 = min(cutoff%disp3, cutoff%disp2)

end function get_realspace_cutoff


!> Measure of the periodic cell, a volume, an area or a length
pure function get_periodic_measure(mol, ndim) result(measure)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Number of periodic dimensions
   integer, intent(in) :: ndim

   !> Measure of the periodic cell
   real(wp) :: measure

   integer :: ii, idx(3)
   real(wp) :: avec(3), bvec(3)

   measure = 1.0_wp
   if (ndim < 1) return

   idx(:ndim) = pack([(ii, ii = 1, 3)], mol%periodic)

   select case(ndim)
   case(3)
      measure = abs(matdet_3x3(mol%lattice))
   case(2)
      ! area from the Gram determinant of the two periodic lattice vectors
      avec(:) = mol%lattice(:, idx(1))
      bvec(:) = mol%lattice(:, idx(2))
      measure = sqrt(max(sum(avec**2)*sum(bvec**2) - dot_product(avec, bvec)**2, &
         & 0.0_wp))
   case default
      measure = norm2(mol%lattice(:, idx(1)))
   end select

end function get_periodic_measure


end module dftd3_accuracy
