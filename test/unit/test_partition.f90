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
      & d3_lowrank_config, rational_damping_param, new_rational_damping, d3_param, &
      & work_partition, new_work_partition, serial_work_partition
   use dftd3_gcp, only : gcp_param, get_gcp_param, get_geometric_counterpoise, &
      & get_geometric_counterpoise_hessian
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   implicit none
   private

   public :: collect_partition

   integer, parameter :: nparts = 3
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)

   type(d3_param), parameter :: pbe0_d3bj = d3_param(&
      & s6=1.0_wp, s8=1.2177_wp, a1=0.4145_wp, a2=4.8593_wp, s9=1.0_wp, alp=14.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_partition(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("invalid partition", test_invalid), &
      & new_unittest("serial partition", test_serial), &
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


!> Summing all parts must reproduce the serial energy, gradient and virial
subroutine test_dispersion_partitioned(error)

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

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   call new_rational_damping(param, pbe0_d3bj)

   allocate(gradient(3, mol%nat), part_gradient(3, mol%nat), sum_gradient(3, mol%nat))
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
   end if

end subroutine test_dispersion_partitioned


subroutine test_hessian_partitioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d3_model) :: d3
   type(rational_damping_param) :: param
   type(work_partition) :: partition
   integer :: part, ndim
   real(wp) :: energy
   real(wp), allocatable :: hessian(:, :), part_hessian(:, :), sum_hessian(:, :)

   call get_structure(mol, "MB16-43", "02")
   call new_d3_model(d3, mol)
   call new_rational_damping(param, pbe0_d3bj)

   ndim = 3*mol%nat
   allocate(hessian(ndim, ndim), part_hessian(ndim, ndim), sum_hessian(ndim, ndim))
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
   end if

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
   end if

end subroutine test_gcp_partitioned


end module test_partition
