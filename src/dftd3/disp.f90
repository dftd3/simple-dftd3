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

module dftd3_disp
   use dftd3_blas, only : d3_gemv
   use dftd3_cutoff, only : realspace_cutoff, get_lattice_points
   use dftd3_damping, only : damping_param
   use dftd3_data, only : get_covalent_rad
   use dftd3_model, only : d3_model
   use dftd3_ncoord, only : get_coordination_number
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   implicit none
   private

   public :: get_dispersion


contains


subroutine get_dispersion(mol, disp, param, cutoff, energy, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion energy
   real(wp), intent(out) :: energy

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   logical :: grad
   integer :: mref
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: dEdcn(:), energies(:)
   real(wp), allocatable :: lattr(:, :)

   mref = maxval(disp%ref)
   grad = present(gradient).and.present(sigma)

   allocate(cn(mol%nat))
   if (grad) allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, cn, dcndr, dcndL)

   allocate(gwvec(mref, mol%nat))
   if (grad) allocate(gwdcn(mref, mol%nat))
   call disp%weight_references(mol, cn, gwvec, gwdcn)

   allocate(c6(mol%nat, mol%nat))
   if (grad) allocate(dc6dcn(mol%nat, mol%nat))
   call disp%get_atomic_c6(mol, gwvec, gwdcn, c6, dc6dcn)

   allocate(energies(mol%nat))
   energies(:) = 0.0_wp
   if (grad) then
      allocate(dEdcn(mol%nat))
      dEdcn(:) = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
   end if
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call param%get_dispersion2(mol, lattr, cutoff%disp2, disp%rvdw, disp%r4r2, c6, dc6dcn, &
      & energies, dEdcn, gradient, sigma)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_dispersion3(mol, lattr, cutoff%disp3, disp%rvdw, disp%r4r2, c6, dc6dcn, &
      & energies, dEdcn, gradient, sigma)
   if (grad) then
      call d3_gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
      call d3_gemv(dcndL, dEdcn, sigma, beta=1.0_wp)
   end if

   energy = sum(energies)

end subroutine get_dispersion


end module dftd3_disp
