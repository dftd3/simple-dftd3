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

!> Check that a calculation distributed over MPI ranks reproduces the serial
!> result, this driver is only built with MPI support and run under mpiexec
program mpi_tester
   use, intrinsic :: iso_fortran_env, only : error_unit
   use dftd3, only : get_dispersion, get_dispersion_mpi, new_mpi_work_partition, &
      & work_partition, realspace_cutoff, d3_model, new_d3_model, d3_param, &
      & rational_damping_param, new_rational_damping, dftd3_has_mpi, &
      & d3_lowrank_config
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mpi_f08, only : MPI_COMM_WORLD, MPI_Init, MPI_Finalize, MPI_Abort, MPI_Comm_rank
   use mstore, only : get_structure
   implicit none

   real(wp), parameter :: thr = 1.0e-9_wp
   type(d3_param), parameter :: pbe0_d3bj = d3_param(&
      & s6=1.0_wp, s8=1.2177_wp, a1=0.4145_wp, a2=4.8593_wp, s9=1.0_wp, alp=14.0_wp)

   type(error_type), allocatable :: error
   type(structure_type) :: mol
   type(d3_model) :: d3
   type(d3_param) :: input
   type(rational_damping_param) :: param
   type(work_partition) :: partition
   real(wp) :: energy, energy_ref, sigma(3, 3), sigma_ref(3, 3)
   real(wp), allocatable :: gradient(:, :), gradient_ref(:, :)
   integer :: stat, rank, idamp

   call MPI_Init(stat)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, stat)

   if (.not.dftd3_has_mpi) call fatal("Library reports MPI support as unavailable")

   ! the library takes plain integer handles from every language binding
   call new_mpi_work_partition(error, partition, MPI_COMM_WORLD%mpi_val)
   if (allocated(error)) call fatal(error%message)

   call get_structure(mol, "MB16-43", "01")
   call new_d3_model(d3, mol)
   allocate(gradient(3, mol%nat), gradient_ref(3, mol%nat))

   ! without the three-body term the C6 coefficients are partitioned as well
   do idamp = 1, 2
      input = pbe0_d3bj
      input%s9 = merge(0.0_wp, 1.0_wp, idamp == 1)
      call new_rational_damping(param, input)

      call check_distributed("real space")
   end do

   ! the reciprocal space summation hands out whole terms of the low-rank
   ! expansion, on a particle mesh and over the reciprocal lattice alike
   call get_structure(mol, "X23", "acetic")
   deallocate(gradient, gradient_ref)
   allocate(gradient(3, mol%nat), gradient_ref(3, mol%nat))
   call new_rational_damping(param, pbe0_d3bj)

   call new_d3_model(d3, mol, lowrank=d3_lowrank_config(tolerance=1.0e-6_wp))
   call check_distributed("particle mesh")

   call new_d3_model(d3, mol, lowrank=d3_lowrank_config(tolerance=1.0e-6_wp, mesh=-1))
   call check_distributed("reciprocal lattice")

   if (rank == 0) write(error_unit, "(a)") "# distributed result matches serial result"

   call MPI_Finalize(stat)

contains


   subroutine check_distributed(label)
      character(len=*), intent(in) :: label

      call get_dispersion(error, mol, d3, param, realspace_cutoff(), energy_ref, &
         & gradient_ref, sigma_ref)
      if (allocated(error)) call fatal(error%message)

      call get_dispersion_mpi(error, mol, d3, param, realspace_cutoff(), &
         & MPI_COMM_WORLD%mpi_val, energy, gradient, sigma)
      if (allocated(error)) call fatal(error%message)

      call expect(abs(energy - energy_ref) < thr, label//" energy")
      call expect(all(abs(gradient - gradient_ref) < thr), label//" gradient")
      call expect(all(abs(sigma - sigma_ref) < thr), label//" virial")
   end subroutine check_distributed

   subroutine expect(cond, label)
      logical, intent(in) :: cond
      character(len=*), intent(in) :: label
      if (.not.cond) call fatal("Distributed "//label//" differs from serial result")
   end subroutine expect

   subroutine fatal(message)
      character(len=*), intent(in) :: message
      write(error_unit, "(2a)") "Error: ", message
      call MPI_Abort(MPI_COMM_WORLD, 1, stat)
   end subroutine fatal

end program mpi_tester
