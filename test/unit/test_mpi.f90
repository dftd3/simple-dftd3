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

module test_mpi
   use dftd3, only : dftd3_has_mpi, dftd3_has_feature, new_mpi_work_partition, &
      & get_dispersion_mpi, work_partition, realspace_cutoff, d3_model, &
      & new_d3_model, d3_param, rational_damping_param, new_rational_damping
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   implicit none
   private

   public :: collect_mpi

   !> Placeholder handle, no communicator is usable without a call to MPI_Init
   integer, parameter :: dummy_comm = 0

   type(d3_param), parameter :: pbe0_d3bj = d3_param(&
      & s6=1.0_wp, s8=1.2177_wp, a1=0.4145_wp, a2=4.8593_wp, s9=1.0_wp, alp=14.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_mpi(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("feature query", test_feature), &
      & new_unittest("uninitialized partition", test_uninitialized_partition), &
      & new_unittest("uninitialized dispersion", test_uninitialized_dispersion) &
      & ]

end subroutine collect_mpi


subroutine test_feature(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check(error, dftd3_has_feature("mpi"), dftd3_has_mpi)
   if (allocated(error)) return

   call check(error, .not.dftd3_has_feature("this-is-not-a-feature"))

end subroutine test_feature


!> Without MPI support the missing feature is reported, with MPI support the
!> uninitialized library is, in neither case may a partition come back
subroutine test_uninitialized_partition(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(error_type), allocatable :: mpi_error
   type(work_partition) :: partition

   call new_mpi_work_partition(mpi_error, partition, dummy_comm)
   if (.not.allocated(mpi_error)) then
      call test_failed(error, "Work partition created without initialized MPI")
      return
   end if

   call check(error, mpi_error%message, expected_message())

end subroutine test_uninitialized_partition


subroutine test_uninitialized_dispersion(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(error_type), allocatable :: mpi_error
   type(structure_type) :: mol
   type(d3_model) :: d3
   type(rational_damping_param) :: param
   real(wp) :: energy

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   call new_rational_damping(param, pbe0_d3bj)

   call get_dispersion_mpi(mpi_error, mol, d3, param, realspace_cutoff(), &
      & dummy_comm, energy)
   if (.not.allocated(mpi_error)) then
      call test_failed(error, "Dispersion evaluated without initialized MPI")
      return
   end if

   call check(error, mpi_error%message, expected_message())

end subroutine test_uninitialized_dispersion


pure function expected_message() result(message)
   character(len=:), allocatable :: message

   if (dftd3_has_mpi) then
      message = "MPI is not initialized"
   else
      message = "s-dftd3 was built without MPI support"
   end if
end function expected_message


end module test_mpi
