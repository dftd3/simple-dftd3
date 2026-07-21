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

module test_output
   use dftd3_output, only : turbomole_gradient, turbomole_gradlatt
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   implicit none
   private

   public :: collect_output

contains


!> Collect all exported unit tests
subroutine collect_output(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("turbomole-gradient-corrupt", test_gradient_corrupt), &
      & new_unittest("turbomole-gradlatt-corrupt", test_gradlatt_corrupt) &
      & ]

end subroutine collect_output


!> A corrupted (non-numeric) gradient entry in an existing $grad file must be
!> detected and reported via stat/=0 instead of being silently ignored.
subroutine test_gradient_corrupt(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 2
   integer, parameter :: num(nat) = [1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.4_wp], &
      & shape(xyz))
   character(len=*), parameter :: fname = "test_output_gradient_corrupt.tmp"
   type(structure_type) :: mol
   real(wp) :: gradient(3, nat)
   integer :: stat, io

   call new(mol, num, xyz)
   gradient = 0.0_wp

   open(newunit=io, file=fname, status="replace", action="write")
   write(io, '("$grad")')
   write(io, '(2x,"cycle =",1x,i6,4x,"SCF energy =",f18.11,3x,"|dE/dxyz| =",f10.6)') &
      & 1, 0.0_wp, 0.0_wp
   write(io, "(3(F20.14,2x),4x,a2)") mol%xyz(1,1), mol%xyz(2,1), mol%xyz(3,1), "H "
   write(io, "(3(F20.14,2x),4x,a2)") mol%xyz(1,2), mol%xyz(2,2), mol%xyz(3,2), "H "
   write(io, '("   xxx   yyy   zzz")')
   write(io, "(3D22.13)") 0.0_wp, 0.0_wp, 0.0_wp
   write(io, '("$end")')
   close(io)

   call turbomole_gradient(mol, fname, 0.0_wp, gradient, stat)

   open(newunit=io, file=fname, status="old")
   close(io, status="delete")

   call check(error, stat, 1)
   if (allocated(error)) return

end subroutine test_gradient_corrupt


!> A corrupted (non-numeric) gradient entry in an existing $gradlatt file must
!> be detected and reported via stat/=0 instead of being silently ignored.
subroutine test_gradlatt_corrupt(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 1
   integer, parameter :: num(nat) = [6]
   real(wp), parameter :: xyz(3, nat) = reshape([0.0_wp, 0.0_wp, 0.0_wp], shape(xyz))
   real(wp), parameter :: lattice(3, 3) = reshape([&
      & 10.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 10.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 10.0_wp], &
      & shape(lattice))
   character(len=*), parameter :: fname = "test_output_gradlatt_corrupt.tmp"
   type(structure_type) :: mol
   real(wp) :: sigma(3, 3)
   integer :: stat, io

   call new(mol, num, xyz, periodic=[.true.], lattice=lattice)
   sigma = 0.0_wp

   open(newunit=io, file=fname, status="replace", action="write")
   write(io, '("$gradlatt")')
   write(io, '(2x,"cycle =",1x,i6,4x,"SCF energy =",f18.11,3x,"|dE/dlatt| =",f10.6)') &
      & 1, 0.0_wp, 0.0_wp
   write(io, "(3(F20.14,2x))") mol%lattice(1,1), mol%lattice(2,1), mol%lattice(3,1)
   write(io, "(3(F20.14,2x))") mol%lattice(1,2), mol%lattice(2,2), mol%lattice(3,2)
   write(io, "(3(F20.14,2x))") mol%lattice(1,3), mol%lattice(2,3), mol%lattice(3,3)
   write(io, '("   xxx   yyy   zzz")')
   write(io, "(3D22.13)") 0.0_wp, 0.0_wp, 0.0_wp
   write(io, "(3D22.13)") 0.0_wp, 0.0_wp, 0.0_wp
   write(io, '("$end")')
   close(io)

   call turbomole_gradlatt(mol, fname, 0.0_wp, sigma, stat)

   open(newunit=io, file=fname, status="old")
   close(io, status="delete")

   call check(error, stat, 1)
   if (allocated(error)) return

end subroutine test_gradlatt_corrupt

end module test_output
