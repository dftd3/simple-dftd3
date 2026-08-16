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
   use, intrinsic :: iso_fortran_env, only : error_unit
   use dftd3_cutoff, only : realspace_cutoff, get_lattice_points
   use dftd3_damping, only : damping_param, get_dispersion2_hessian
   use dftd3_model, only : d3_model
   use dftd3_ncoord, only : get_coordination_number, get_partitioned_coordination_number, &
      & add_coordination_number_derivs, add_coordination_number_hessian
   use dftd3_partition, only : work_partition, work_reducer
   use mctc_data, only : get_covalent_rad
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   implicit none
   private

   public :: get_dispersion, get_pairwise_dispersion


   !> Calculate dispersion energy
   interface get_dispersion
      module procedure :: get_dispersion_atomic
      module procedure :: get_dispersion_scalar
      module procedure :: get_dispersion_atomic_v2
      module procedure :: get_dispersion_scalar_v2
   end interface get_dispersion

   !> Calculate pairwise representation of the dispersion energy
   interface get_pairwise_dispersion
      module procedure :: get_pairwise_dispersion
      module procedure :: get_pairwise_dispersion_v2
   end interface get_pairwise_dispersion

contains


!> Calculate atom-resolved dispersion energies.
!>
!> The dispersion model and the damping parameters have to agree on the summation
!> technique, an inconsistent setup is reported in the error handler.
subroutine get_dispersion_atomic_v2(error, mol, disp, param, cutoff, energies, &
      & gradient, sigma, hessian, partition, reducer)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion energy
   real(wp), intent(out) :: energies(:)

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   !> Dispersion hessian
   real(wp), intent(out), contiguous, optional :: hessian(:, :)

   !> Work partition of the interaction loops, see dftd3_partition
   type(work_partition), intent(in), optional :: partition

   !> Communication backend, enables partitioning the coordination number
   class(work_reducer), intent(in), optional :: reducer

   logical :: grad, hess, ewald
   integer :: mref
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :)
   real(wp), allocatable :: dEdcn(:)
   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)

   mref = maxval(disp%ref)
   grad = present(gradient) .or. present(sigma)
   hess = present(hessian)

   ! a low-rank model under three-dimensional boundary conditions asks for the
   ! reciprocal space summation, which the damping function has to support
   ewald = allocated(disp%lowrank) .and. all(mol%periodic)
   if (ewald .and. .not.param%supports_ewald()) then
      call fatal_error(error, "Damping function does not support Ewald summation")
      return
   end if

   ! the low-rank expansion has no second derivatives, and for a periodic model
   ! they would mix the reciprocal space energy with a real space curvature
   if (allocated(disp%lowrank) .and. hess) then
      call fatal_error(error, "Hessian is only available for the real space summation")
      return
   end if

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   if (present(reducer)) then
      call get_partitioned_coordination_number(mol, lattr, cutoff%cn, disp%rcov, cn, &
         & partition)
      call reducer%reduce(cn, error)
      if (allocated(error)) return
   else
      call get_partitioned_coordination_number(mol, lattr, cutoff%cn, disp%rcov, cn)
   end if

   allocate(gwvec(mref, mol%nat))
   if (grad) allocate(gwdcn(mref, mol%nat))
   call disp%weight_references(mol, cn, gwvec, gwdcn)

   allocate(c6(mol%nat, mol%nat))
   if (grad) allocate(dc6dcn(mol%nat, mol%nat))
   call disp%get_atomic_c6(mol, gwvec, gwdcn, c6, dc6dcn)

   energies(:) = 0.0_wp
   if (grad) then
      allocate(dEdcn(mol%nat))
      dEdcn(:) = 0.0_wp
      if (present(gradient)) gradient(:, :) = 0.0_wp
      if (present(sigma)) sigma(:, :) = 0.0_wp
      allocate(gradient_local(3, mol%nat), source=0.0_wp)
      allocate(sigma_local(3, 3), source=0.0_wp)
   end if
   if (ewald) then
      call param%get_dispersion2_ewald(mol, disp, gwvec, gwdcn, energies, dEdcn, &
         & gradient_local, sigma_local, error, partition)
      if (allocated(error)) return
   else
      call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
      call param%get_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, &
         & disp%rvdw, disp%r4r2, c6, dc6dcn, energies, dEdcn, gradient_local, &
         & sigma_local, partition)
   end if
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_dispersion3(mol, lattr, cutoff%disp3, cutoff%width3, &
      & disp%rvdw, disp%r4r2, c6, dc6dcn, energies, dEdcn, gradient_local, &
      & sigma_local, partition)
   if (grad) then
      call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
      if (present(reducer)) then
         call reducer%reduce(dEdcn, error)
         if (allocated(error)) return
         call add_coordination_number_derivs(mol, lattr, cutoff%cn, disp%rcov, dEdcn, &
            & gradient_local, sigma_local, partition)
      else
         call add_coordination_number_derivs(mol, lattr, cutoff%cn, disp%rcov, dEdcn, &
            & gradient_local, sigma_local)
      end if
      if (present(gradient)) gradient(:, :) = gradient_local(:, :)
      if (present(sigma)) sigma(:, :) = sigma_local(:, :)
   end if
   if (hess) then
      call get_dispersion_hessian(mol, disp, param, cutoff, hessian, partition)
   end if

end subroutine get_dispersion_atomic_v2


!> Calculate atom-resolved dispersion energies.
!>
!> deprecated: removed with the v2 API, use the interface with error handling
subroutine get_dispersion_atomic(mol, disp, param, cutoff, energies, gradient, sigma, hessian)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion energy
   real(wp), intent(out) :: energies(:)

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   !> Dispersion hessian
   real(wp), intent(out), contiguous, optional :: hessian(:, :)

   type(error_type), allocatable :: error

   call get_dispersion_atomic_v2(error, mol, disp, param, cutoff, energies, &
      & gradient, sigma, hessian)

   ! this interface cannot propagate the inconsistent setup to the caller
   if (allocated(error)) then
      write(error_unit, '("[Fatal]", 1x, a)') error%message
      error stop
   end if

end subroutine get_dispersion_atomic


!> Analytical second derivatives of the dispersion energy w.r.t. the coordinates.
!>
!> The energy depends on the coordinates directly and through the coordination
!> number. Both contributions are accumulated separately and the coordination
!> number part is contracted with dCN/dR afterwards.
subroutine get_dispersion_hessian(mol, disp, param, cutoff, hessian, partition)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion hessian
   real(wp), intent(out) :: hessian(:, :)

   !> Work partition of the interaction loops
   type(work_partition), intent(in), optional :: partition

   integer :: mref, nat, ndim, iat, ic, kat
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: gwvec(:, :), gwdcn(:, :), gwd2cn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), d2c6dcn2(:, :), d2c6dcnij(:, :)
   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: dEdcn(:), dEdcndr(:, :), dEdcndcn(:, :), dr(:, :)

   nat = mol%nat
   ndim = 3*nat
   mref = maxval(disp%ref)

   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, cn, dcndr, dcndL)

   allocate(gwvec(mref, nat), gwdcn(mref, nat), gwd2cn(mref, nat))
   call disp%weight_references(mol, cn, gwvec, gwdcn, gwd2cn)

   allocate(c6(nat, nat), dc6dcn(nat, nat), d2c6dcn2(nat, nat), d2c6dcnij(nat, nat))
   call disp%get_atomic_c6(mol, gwvec, gwdcn, c6, dc6dcn, gwd2cn, d2c6dcn2, d2c6dcnij)

   hessian(:, :) = 0.0_wp
   allocate(dEdcn(nat), source=0.0_wp)
   allocate(dEdcndr(ndim, nat), source=0.0_wp)
   allocate(dEdcndcn(nat, nat), source=0.0_wp)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
   call get_dispersion2_hessian(param, mol, lattr, cutoff%disp2, cutoff%width2, &
      & disp%rvdw, disp%r4r2, c6, dc6dcn, d2c6dcn2, d2c6dcnij, hessian, dEdcn, &
      & dEdcndr, dEdcndcn, partition)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call param%get_dispersion3_hessian(mol, lattr, cutoff%disp3, &
      & cutoff%width3, disp%rvdw, disp%r4r2, c6, dc6dcn, d2c6dcn2, d2c6dcnij, &
      & hessian, dEdcn, dEdcndr, dEdcndcn, partition)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call add_coordination_number_hessian(mol, lattr, cutoff%cn, disp%rcov, dEdcn, hessian)

   allocate(dr(ndim, nat))
   do kat = 1, nat
      do iat = 1, nat
         do ic = 1, 3
            dr(3*(iat - 1) + ic, kat) = dcndr(ic, iat, kat)
         end do
      end do
   end do

   hessian(:, :) = hessian + matmul(dEdcndr, transpose(dr)) &
      & + matmul(dr, transpose(dEdcndr)) &
      & + matmul(dr, matmul(dEdcndcn, transpose(dr)))

end subroutine get_dispersion_hessian


!> Calculate scalar dispersion energy.
!>
!> deprecated: removed with the v2 API, use the interface with error handling
subroutine get_dispersion_scalar(mol, disp, param, cutoff, energy, gradient, sigma, hessian)

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

   !> Dispersion hessian
   real(wp), intent(out), contiguous, optional :: hessian(:, :)

   real(wp), allocatable :: energies(:)

   allocate(energies(mol%nat))

   call get_dispersion_atomic(mol, disp, param, cutoff, energies, gradient, sigma, hessian)

   energy = sum(energies)

end subroutine get_dispersion_scalar


!> Calculate scalar dispersion energy, reporting an inconsistent setup.
subroutine get_dispersion_scalar_v2(error, mol, disp, param, cutoff, energy, &
      & gradient, sigma, hessian, partition, reducer)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

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

   !> Dispersion hessian
   real(wp), intent(out), contiguous, optional :: hessian(:, :)

   !> Work partition of the interaction loops, see dftd3_partition
   type(work_partition), intent(in), optional :: partition

   !> Communication backend, enables partitioning the coordination number
   class(work_reducer), intent(in), optional :: reducer

   real(wp), allocatable :: energies(:)

   allocate(energies(mol%nat))
   energy = 0.0_wp

   call get_dispersion_atomic_v2(error, mol, disp, param, cutoff, energies, &
      & gradient, sigma, hessian, partition, reducer)
   if (allocated(error)) return

   energy = sum(energies)

end subroutine get_dispersion_scalar_v2


!> Calculate the pairwise representation
subroutine get_pairwise_dispersion_v2(error, mol, disp, param, cutoff, energy2, energy3)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Pairwise representation of additive dispersion energy
   real(wp), intent(out) :: energy2(:, :)

   !> Pairwise representation of non-additive dispersion energy
   real(wp), intent(out) :: energy3(:, :)

   integer :: mref
   real(wp), allocatable :: cn(:), gwvec(:, :), c6(:, :), lattr(:, :)

   if (allocated(disp%lowrank) .and. all(mol%periodic)) then
      call fatal_error(error, "Pairwise analysis is only available for the "//&
         & "real space summation")
      return
   end if

   mref = maxval(disp%ref)

   allocate(cn(mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, cn)

   allocate(gwvec(mref, mol%nat))
   call disp%weight_references(mol, cn, gwvec)

   allocate(c6(mol%nat, mol%nat))
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   energy2(:, :) = 0.0_wp
   energy3(:, :) = 0.0_wp
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
    call param%get_pairwise_dispersion2(mol, lattr, cutoff%disp2, cutoff%width2, &
       & disp%rvdw, disp%r4r2, c6, energy2)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
    call param%get_pairwise_dispersion3(mol, lattr, cutoff%disp3, cutoff%width3, &
       & disp%rvdw, disp%r4r2, c6, energy3)

end subroutine get_pairwise_dispersion_v2


!> Wrapper to handle the evaluation of pairwise representation of the dispersion energy
!>
!> deprecated: removed with the v2 API, use the interface with error handling
subroutine get_pairwise_dispersion(mol, disp, param, cutoff, energy2, energy3)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Pairwise representation of additive dispersion energy
   real(wp), intent(out) :: energy2(:, :)

   !> Pairwise representation of non-additive dispersion energy
   real(wp), intent(out) :: energy3(:, :)

   type(error_type), allocatable :: error

   call get_pairwise_dispersion_v2(error, mol, disp, param, cutoff, energy2, energy3)

   if (allocated(error)) then
      write(error_unit, '("[Fatal]", 1x, a)') error%message
      error stop
   end if

end subroutine get_pairwise_dispersion


end module dftd3_disp
