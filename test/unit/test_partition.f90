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
   integer, parameter :: ndamping = 6
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
