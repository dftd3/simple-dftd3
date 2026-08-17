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

module test_partition
   use dftd3, only : get_dispersion, realspace_cutoff, d3_model, new_d3_model, &
      & d3_lowrank_config, damping_param, d3_param, &
      & rational_damping_param, new_rational_damping, &
      & zero_damping_param, new_zero_damping, &
      & mzero_damping_param, new_mzero_damping, &
      & optimizedpower_damping_param, new_optimizedpower_damping, &
      & cso_damping_param, new_cso_damping, &
      & z_damping_param, new_z_damping, &
      & work_partition, new_work_partition, serial_work_partition, work_reducer, &
      & get_lattice_points
   use dftd3_gcp, only : gcp_param, get_gcp_param, get_geometric_counterpoise, &
      & get_geometric_counterpoise_hessian
   use dftd3_ncoord, only : get_partitioned_coordination_number, &
      & add_coordination_number_derivs
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   use mstore, only : get_structure
   implicit none
   private

   public :: collect_partition

   integer, parameter :: nparts = 3
   integer, parameter :: ndamping = 6
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   !> Summing the parts regroups the additions, this is well below a missing pair
   real(wp), parameter :: cnthr = 1.0e-10_wp

   type(d3_param), parameter :: pbe0_d3bj = d3_param(&
      & s6=1.0_wp, s8=1.2177_wp, a1=0.4145_wp, a2=4.8593_wp, s9=1.0_wp, alp=14.0_wp)

   !> Reducer over a single part, exercises the reduction path without MPI
   type, extends(work_reducer) :: identity_reducer
   contains
      procedure :: reduce => identity_reduce
   end type identity_reducer


contains


!> Collect all exported unit tests
subroutine collect_partition(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("invalid partition", test_invalid), &
      & new_unittest("serial partition", test_serial), &
      & new_unittest("coordination number", test_coordination_number), &
      & new_unittest("c6 on demand", test_c6_on_demand), &
      & new_unittest("reducer", test_reducer), &
      & new_unittest("dispersion", test_dispersion_partitioned), &
      & new_unittest("hessian", test_hessian_partitioned), &
      & new_unittest("ewald", test_ewald_partitioned), &
      & new_unittest("counterpoise", test_gcp_partitioned) &
      & ]

end subroutine collect_partition


subroutine test_invalid(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(work_partition) :: partition
   type(error_type), allocatable :: partition_error

   call new_work_partition(partition_error, partition, -1, nparts)
   if (.not.allocated(partition_error)) then
      call test_failed(error, "Negative part index did not return an error")
      return
   end if

   call new_work_partition(partition_error, partition, nparts, nparts)
   if (.not.allocated(partition_error)) then
      call test_failed(error, "Out of range part index did not return an error")
      return
   end if

   call new_work_partition(partition_error, partition, 0, 0)
   if (.not.allocated(partition_error)) then
      call test_failed(error, "Vanishing number of parts did not return an error")
   end if

end subroutine test_invalid


!> Passing the serial partition must be identical to omitting the argument
subroutine test_serial(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   type(rational_damping_param) :: param
   type(gcp_param) :: gcp
   real(wp) :: energy, serial_energy
   real(wp) :: sigma(3, 3), serial_sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), serial_gradient(:, :)

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   call new_rational_damping(param, pbe0_d3bj)

   allocate(gradient(3, mol%nat), serial_gradient(3, mol%nat))
   call get_dispersion(error, mol, d3, param, realspace_cutoff(), energy, gradient, sigma)
   if (allocated(error)) return
   call get_dispersion(error, mol, d3, param, realspace_cutoff(), serial_energy, &
      & serial_gradient, serial_sigma, partition=serial_work_partition)
   if (allocated(error)) return

   call check(error, serial_energy, energy, thr=thr)
   if (allocated(error)) return

   if (any(abs(serial_gradient - gradient) > thr) .or. &
      & any(abs(serial_sigma - sigma) > thr)) then
      call test_failed(error, "Serial dispersion partition does not match")
      return
   end if

   call get_gcp_param(gcp, mol, method="hf3c")

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call get_geometric_counterpoise(mol, gcp, realspace_cutoff(), energy, gradient, sigma)

   serial_energy = 0.0_wp
   serial_gradient(:, :) = 0.0_wp
   serial_sigma(:, :) = 0.0_wp
   call get_geometric_counterpoise(mol, gcp, realspace_cutoff(), serial_energy, &
      & serial_gradient, serial_sigma, serial_work_partition)

   call check(error, serial_energy, energy, thr=thr)
   if (allocated(error)) return

   if (any(abs(serial_gradient - gradient) > thr) .or. &
      & any(abs(serial_sigma - sigma) > thr)) then
      call test_failed(error, "Serial counterpoise partition does not match")
   end if

end subroutine test_serial


!> The local coordination number loops must reproduce mctc-lib and must sum up
!> over the parts
subroutine test_coordination_number(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_coordination_number(error, "MB16-43", "01")
   if (allocated(error)) return

   call check_coordination_number(error, "X23", "acetic")

end subroutine test_coordination_number


subroutine check_coordination_number(error, set, id)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Structure set and identifier
   character(len=*), intent(in) :: set, id

   type(structure_type) :: mol
   type(d3_model) :: d3
   type(realspace_cutoff) :: cutoff
   type(work_partition) :: partition
   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: ncoord_error
   real(wp), allocatable :: lattr(:, :), cn(:), cn_ref(:), cn_part(:), dEdcn(:)
   real(wp), allocatable :: gradient(:, :), gradient_ref(:, :)
   real(wp) :: sigma(3, 3), sigma_ref(3, 3)
   integer :: ipart, iat

   call get_structure(mol, set, id)
   call new_d3_model(d3, mol)
   cutoff = realspace_cutoff()
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)

   ! steepness of the counting function, private in dftd3_ncoord
   call new_ncoord(ncoord, mol, cn_count%exp, ncoord_error, kcn=16.0_wp, &
      & cutoff=cutoff%cn, rcov=d3%rcov)
   if (allocated(ncoord_error)) then
      call test_failed(error, ncoord_error%message)
      return
   end if

   allocate(cn(mol%nat), cn_ref(mol%nat), cn_part(mol%nat), dEdcn(mol%nat))
   call ncoord%get_coordination_number(mol, lattr, cn_ref)
   call get_partitioned_coordination_number(mol, lattr, cutoff%cn, d3%rcov, cn)
   if (any(abs(cn - cn_ref) > cnthr)) then
      call test_failed(error, "Coordination number does not match mctc-lib for "//set)
      return
   end if

   cn(:) = 0.0_wp
   do ipart = 0, nparts - 1
      call new_work_partition(error, partition, ipart, nparts)
      if (allocated(error)) return
      call get_partitioned_coordination_number(mol, lattr, cutoff%cn, d3%rcov, &
         & cn_part, partition)
      cn(:) = cn + cn_part
   end do
   if (any(abs(cn - cn_ref) > cnthr)) then
      call test_failed(error, "Partitioned coordination number does not sum up for "//set)
      return
   end if

   do iat = 1, mol%nat
      dEdcn(iat) = 0.01_wp*iat
   end do
   allocate(gradient(3, mol%nat), gradient_ref(3, mol%nat), source=0.0_wp)
   sigma(:, :) = 0.0_wp
   sigma_ref(:, :) = 0.0_wp

   call ncoord%add_coordination_number_derivs(mol, lattr, dEdcn, gradient_ref, sigma_ref)
   do ipart = 0, nparts - 1
      call new_work_partition(error, partition, ipart, nparts)
      if (allocated(error)) return
      call add_coordination_number_derivs(mol, lattr, cutoff%cn, d3%rcov, dEdcn, &
         & gradient, sigma, partition)
   end do
   if (any(abs(gradient - gradient_ref) > cnthr) .or. &
      & any(abs(sigma - sigma_ref) > cnthr)) then
      call test_failed(error, "Partitioned CN derivatives do not sum up for "//set)
   end if

end subroutine check_coordination_number


!> Without the three-body term the C6 coefficients are only evaluated for the
!> pairs a part consumes, summing the parts has to restore the complete matrix
subroutine test_c6_on_demand(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   type(d3_param) :: input
   type(rational_damping_param) :: param
   type(realspace_cutoff) :: cutoff
   type(work_partition) :: partition
   real(wp), allocatable :: lattr(:, :), cn(:), gwvec(:, :), gwdcn(:, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), c6_ref(:, :), dc6dcn_ref(:, :)
   real(wp), allocatable :: c6_sum(:, :), dc6dcn_sum(:, :)
   real(wp), allocatable :: gradient(:, :), part_gradient(:, :), sum_gradient(:, :)
   real(wp) :: energy, part_energy, sum_energy
   real(wp) :: sigma(3, 3), part_sigma(3, 3), sum_sigma(3, 3)
   integer :: mref, ipart

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   cutoff = realspace_cutoff()
   mref = maxval(d3%ref)

   allocate(cn(mol%nat), gwvec(mref, mol%nat), gwdcn(mref, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
   call get_partitioned_coordination_number(mol, lattr, cutoff%cn, d3%rcov, cn)
   call d3%weight_references(mol, cn, gwvec, gwdcn)

   allocate(c6(mol%nat, mol%nat), dc6dcn(mol%nat, mol%nat))
   allocate(c6_ref(mol%nat, mol%nat), dc6dcn_ref(mol%nat, mol%nat))
   allocate(c6_sum(mol%nat, mol%nat), dc6dcn_sum(mol%nat, mol%nat), source=0.0_wp)

   call d3%get_atomic_c6(mol, gwvec, gwdcn, c6_ref, dc6dcn_ref)
   do ipart = 0, nparts - 1
      call new_work_partition(error, partition, ipart, nparts)
      if (allocated(error)) return
      call d3%get_atomic_c6(mol, gwvec, gwdcn, c6, dc6dcn, partition=partition)
      c6_sum(:, :) = c6_sum + c6
      dc6dcn_sum(:, :) = dc6dcn_sum + dc6dcn
   end do

   if (any(abs(c6_sum - c6_ref) > thr) .or. &
      & any(abs(dc6dcn_sum - dc6dcn_ref) > thr)) then
      call test_failed(error, "Partitioned C6 coefficients do not sum up")
      return
   end if

   ! the two-body only calculation must not notice the missing coefficients
   input = pbe0_d3bj
   input%s9 = 0.0_wp
   call new_rational_damping(param, input)

   allocate(gradient(3, mol%nat), part_gradient(3, mol%nat), sum_gradient(3, mol%nat))
   call get_dispersion(error, mol, d3, param, cutoff, energy, gradient, sigma)
   if (allocated(error)) return

   sum_energy = 0.0_wp
   sum_gradient(:, :) = 0.0_wp
   sum_sigma(:, :) = 0.0_wp
   do ipart = 0, nparts - 1
      call new_work_partition(error, partition, ipart, nparts)
      if (allocated(error)) return
      call get_dispersion(error, mol, d3, param, cutoff, part_energy, part_gradient, &
         & part_sigma, partition=partition)
      if (allocated(error)) return
      sum_energy = sum_energy + part_energy
      sum_gradient(:, :) = sum_gradient + part_gradient
      sum_sigma(:, :) = sum_sigma + part_sigma
   end do

   call check(error, sum_energy, energy, thr=thr)
   if (allocated(error)) return

   if (any(abs(sum_gradient - gradient) > thr) .or. &
      & any(abs(sum_sigma - sigma) > thr)) then
      call test_failed(error, "Partitioned two-body dispersion does not sum up")
   end if

end subroutine test_c6_on_demand


!> Reducing over a single part must not change the result
subroutine test_reducer(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   type(rational_damping_param) :: param
   type(identity_reducer) :: reducer
   real(wp) :: energy, ref_energy, sigma(3, 3), ref_sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), ref_gradient(:, :)

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   call new_rational_damping(param, pbe0_d3bj)
   allocate(gradient(3, mol%nat), ref_gradient(3, mol%nat))

   call get_dispersion(error, mol, d3, param, realspace_cutoff(), ref_energy, &
      & ref_gradient, ref_sigma)
   if (allocated(error)) return

   call get_dispersion(error, mol, d3, param, realspace_cutoff(), energy, &
      & gradient, sigma, partition=serial_work_partition, reducer=reducer)
   if (allocated(error)) return

   call check(error, energy, ref_energy, thr=thr)
   if (allocated(error)) return

   if (any(abs(gradient - ref_gradient) > thr) .or. &
      & any(abs(sigma - ref_sigma) > thr)) then
      call test_failed(error, "Reduced serial partition does not match")
   end if

end subroutine test_reducer


subroutine identity_reduce(self, val, error)

   !> Communication backend
   class(identity_reducer), intent(in) :: self

   !> Values to sum over all parts
   real(wp), intent(inout), contiguous :: val(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

end subroutine identity_reduce


!> Damping parameters covering every implementation of the partitioned loops
subroutine make_damping(param, idamp)

   !> Damping parameters
   class(damping_param), allocatable, intent(out) :: param

   !> Selector of the damping function
   integer, intent(in) :: idamp

   type(rational_damping_param) :: rational
   type(zero_damping_param) :: zero
   type(mzero_damping_param) :: mzero
   type(optimizedpower_damping_param) :: optimizedpower
   type(cso_damping_param) :: cso
   type(z_damping_param) :: z

   select case(idamp)
   case(1)
      call new_rational_damping(rational, d3_param(&
         & s6=1.0_wp, s8=1.2177_wp, a1=0.4145_wp, a2=4.8593_wp, s9=1.0_wp, alp=14.0_wp))
      param = rational
   case(2)
      call new_zero_damping(zero, d3_param(&
         & s6=1.0_wp, s9=1.0_wp, alp=14.0_wp, rs8=1.0_wp, rs6=1.094_wp, s8=1.682_wp))
      param = zero
   case(3)
      call new_mzero_damping(mzero, d3_param(&
         & s6=1.0_wp, s9=1.0_wp, alp=14.0_wp, rs8=1.0_wp, rs6=2.340218_wp, s8=0.0_wp, &
         & bet=0.129434_wp))
      param = mzero
   case(4)
      call new_optimizedpower_damping(optimizedpower, d3_param(&
         & s6=0.97388_wp, s9=1.0_wp, alp=14.0_wp, bet=6.0_wp, a1=0.150_wp, s8=0.0_wp, &
         & a2=4.25_wp))
      param = optimizedpower
   case(5)
      call new_cso_damping(cso, d3_param(&
         & s6=1.0_wp, s9=1.0_wp, alp=14.0_wp, a1=0.86_wp, a2=2.5_wp, rs6=0.0_wp, rs8=6.25_wp))
      param = cso
   case(6)
      call new_z_damping(z, d3_param(&
         & s6=1.0_wp, s9=1.0_wp, alp=14.0_wp, a1=200770.0_wp, a2=0.0_wp, rs6=0.0_wp, &
         & rs8=6.25_wp))
      param = z
   case default
      continue
   end select

end subroutine make_damping


!> Summing all parts must reproduce the complete result for every damping function,
!> both for the energy-only and the derivative code path
subroutine test_dispersion_partitioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   class(damping_param), allocatable :: param
   type(work_partition) :: partition
   integer :: idamp, part
   real(wp) :: energy, part_energy, sum_energy
   real(wp) :: sigma(3, 3), part_sigma(3, 3), sum_sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), part_gradient(:, :), sum_gradient(:, :)
   real(wp), allocatable :: energies(:), part_energies(:), sum_energies(:)

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)

   allocate(gradient(3, mol%nat), part_gradient(3, mol%nat), sum_gradient(3, mol%nat))
   allocate(energies(mol%nat), part_energies(mol%nat), sum_energies(mol%nat))

   do idamp = 1, ndamping
      call make_damping(param, idamp)

      ! atom-resolved energy without derivatives
      call get_dispersion(error, mol, d3, param, realspace_cutoff(), energies)
      if (allocated(error)) return

      sum_energies(:) = 0.0_wp
      do part = 0, nparts - 1
         call new_work_partition(error, partition, part, nparts)
         if (allocated(error)) return
         call get_dispersion(error, mol, d3, param, realspace_cutoff(), part_energies, &
            & partition=partition)
         if (allocated(error)) return
         sum_energies(:) = sum_energies + part_energies
      end do

      if (any(abs(sum_energies - energies) > thr)) then
         call test_failed(error, "Partitioned atomic dispersion energy does not match")
         return
      end if

      call get_dispersion(error, mol, d3, param, realspace_cutoff(), energy, gradient, sigma)
      if (allocated(error)) return

      sum_energy = 0.0_wp
      sum_gradient(:, :) = 0.0_wp
      sum_sigma(:, :) = 0.0_wp
      do part = 0, nparts - 1
         call new_work_partition(error, partition, part, nparts)
         if (allocated(error)) return
         call get_dispersion(error, mol, d3, param, realspace_cutoff(), part_energy, &
            & part_gradient, part_sigma, partition=partition)
         if (allocated(error)) return
         sum_energy = sum_energy + part_energy
         sum_gradient(:, :) = sum_gradient + part_gradient
         sum_sigma(:, :) = sum_sigma + part_sigma
      end do

      call check(error, sum_energy, energy, thr=thr)
      if (allocated(error)) return

      if (any(abs(sum_gradient - gradient) > thr)) then
         call test_failed(error, "Partitioned dispersion gradient does not match")
         return
      end if

      if (any(abs(sum_sigma - sigma) > thr)) then
         call test_failed(error, "Partitioned dispersion virial does not match")
         return
      end if
   end do

end subroutine test_dispersion_partitioned


subroutine test_hessian_partitioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   class(damping_param), allocatable :: param
   type(work_partition) :: partition
   integer :: idamp, part, ndim
   real(wp) :: energy
   real(wp), allocatable :: hessian(:, :), part_hessian(:, :), sum_hessian(:, :)

   call get_structure(mol, "MB16-43", "02")
   call new_d3_model(d3, mol)

   ndim = 3*mol%nat
   allocate(hessian(ndim, ndim), part_hessian(ndim, ndim), sum_hessian(ndim, ndim))

   do idamp = 1, ndamping
      call make_damping(param, idamp)

      call get_dispersion(error, mol, d3, param, realspace_cutoff(), energy, hessian=hessian)
      if (allocated(error)) return

      sum_hessian(:, :) = 0.0_wp
      do part = 0, nparts - 1
         call new_work_partition(error, partition, part, nparts)
         if (allocated(error)) return
         call get_dispersion(error, mol, d3, param, realspace_cutoff(), energy, &
            & hessian=part_hessian, partition=partition)
         if (allocated(error)) return
         sum_hessian(:, :) = sum_hessian + part_hessian
      end do

      if (any(abs(sum_hessian - hessian) > thr)) then
         call test_failed(error, "Partitioned dispersion hessian does not match")
         return
      end if
   end do

end subroutine test_hessian_partitioned


!> The reciprocal space summation is partitioned over the lattice points
subroutine test_ewald_partitioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   type(rational_damping_param) :: param
   type(work_partition) :: partition
   integer :: part
   real(wp) :: energy, part_energy, sum_energy
   real(wp) :: sigma(3, 3), part_sigma(3, 3), sum_sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), part_gradient(:, :), sum_gradient(:, :)

   call get_structure(mol, "X23", "acetic")
   call new_d3_model(d3, mol, lowrank=d3_lowrank_config(tolerance=1.0e-6_wp))
   call new_rational_damping(param, pbe0_d3bj)

   allocate(gradient(3, mol%nat), part_gradient(3, mol%nat), sum_gradient(3, mol%nat))
   call get_dispersion(error, mol, d3, param, realspace_cutoff(), energy, gradient, sigma)
   if (allocated(error)) return

   ! an explicit serial partition takes the same branches as a distributed one
   call get_dispersion(error, mol, d3, param, realspace_cutoff(), part_energy, &
      & part_gradient, part_sigma, partition=serial_work_partition)
   if (allocated(error)) return
   call check(error, part_energy, energy, thr=thr)
   if (allocated(error)) return

   sum_energy = 0.0_wp
   sum_gradient(:, :) = 0.0_wp
   sum_sigma(:, :) = 0.0_wp
   do part = 0, nparts - 1
      call new_work_partition(error, partition, part, nparts)
      if (allocated(error)) return
      call get_dispersion(error, mol, d3, param, realspace_cutoff(), part_energy, &
         & part_gradient, part_sigma, partition=partition)
      if (allocated(error)) return
      sum_energy = sum_energy + part_energy
      sum_gradient(:, :) = sum_gradient + part_gradient
      sum_sigma(:, :) = sum_sigma + part_sigma
   end do

   call check(error, sum_energy, energy, thr=thr)
   if (allocated(error)) return

   if (any(abs(sum_gradient - gradient) > thr)) then
      call test_failed(error, "Partitioned Ewald gradient does not match")
      return
   end if

   if (any(abs(sum_sigma - sigma) > thr)) then
      call test_failed(error, "Partitioned Ewald virial does not match")
   end if

end subroutine test_ewald_partitioned


subroutine test_gcp_partitioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(gcp_param) :: param
   type(work_partition) :: partition
   integer :: part, ndim
   real(wp) :: energy, part_energy, sum_energy
   real(wp) :: sigma(3, 3), part_sigma(3, 3), sum_sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), part_gradient(:, :), sum_gradient(:, :)
   real(wp), allocatable :: hessian(:, :), part_hessian(:, :), sum_hessian(:, :)

   call get_structure(mol, "MB16-43", "01")
   call get_gcp_param(param, mol, method="hf3c")

   ndim = 3*mol%nat
   allocate(gradient(3, mol%nat), part_gradient(3, mol%nat), sum_gradient(3, mol%nat))
   allocate(hessian(ndim, ndim), part_hessian(ndim, ndim), sum_hessian(ndim, ndim))

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call get_geometric_counterpoise(mol, param, realspace_cutoff(), energy, gradient, sigma)
   call get_geometric_counterpoise_hessian(mol, param, realspace_cutoff(), hessian)

   sum_energy = 0.0_wp
   sum_gradient(:, :) = 0.0_wp
   sum_sigma(:, :) = 0.0_wp
   sum_hessian(:, :) = 0.0_wp
   do part = 0, nparts - 1
      call new_work_partition(error, partition, part, nparts)
      if (allocated(error)) return
      part_energy = 0.0_wp
      part_gradient(:, :) = 0.0_wp
      part_sigma(:, :) = 0.0_wp
      call get_geometric_counterpoise(mol, param, realspace_cutoff(), part_energy, &
         & part_gradient, part_sigma, partition)
      call get_geometric_counterpoise_hessian(mol, param, realspace_cutoff(), &
         & part_hessian, partition)
      sum_energy = sum_energy + part_energy
      sum_gradient(:, :) = sum_gradient + part_gradient
      sum_sigma(:, :) = sum_sigma + part_sigma
      sum_hessian(:, :) = sum_hessian + part_hessian
   end do

   call check(error, sum_energy, energy, thr=thr)
   if (allocated(error)) return

   if (any(abs(sum_gradient - gradient) > thr)) then
      call test_failed(error, "Partitioned counterpoise gradient does not match")
      return
   end if

   if (any(abs(sum_sigma - sigma) > thr)) then
      call test_failed(error, "Partitioned counterpoise virial does not match")
      return
   end if

   if (any(abs(sum_hessian - hessian) > thr)) then
      call test_failed(error, "Partitioned counterpoise hessian does not match")
      return
   end if

   ! the energy-only path uses a separate set of loops
   energy = 0.0_wp
   call get_geometric_counterpoise(mol, param, realspace_cutoff(), energy)

   sum_energy = 0.0_wp
   do part = 0, nparts - 1
      call new_work_partition(error, partition, part, nparts)
      if (allocated(error)) return
      part_energy = 0.0_wp
      call get_geometric_counterpoise(mol, param, realspace_cutoff(), part_energy, &
         & partition=partition)
      sum_energy = sum_energy + part_energy
   end do

   call check(error, sum_energy, energy, thr=thr)

end subroutine test_gcp_partitioned


end module test_partition
