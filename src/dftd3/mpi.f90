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

!> MPI aware entry points of the dispersion calculation.
!>
!> Every rank holds the complete structure, evaluates the share of the
!> interaction loops assigned to its rank and the results are reduced over the
!> communicator. All routines report the missing feature in the error handler if
!> s-dftd3 was built without MPI support, see dftd3_feature.
module dftd3_mpi
   use dftd3_cutoff, only : realspace_cutoff
   use dftd3_damping, only : damping_param
   use dftd3_disp, only : get_dispersion
   use dftd3_model, only : d3_model
   use dftd3_mpi_utils, only : get_mpi_comm_info, mpi_allreduce_sum
   use dftd3_partition, only : work_partition, new_work_partition
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   implicit none
   private

   public :: new_mpi_work_partition, get_dispersion_mpi


   !> Calculate dispersion energy distributed over an MPI communicator
   interface get_dispersion_mpi
      module procedure :: get_dispersion_mpi_atomic
      module procedure :: get_dispersion_mpi_scalar
   end interface get_dispersion_mpi


contains


!> Create the work partition of the calling rank in a communicator
subroutine new_mpi_work_partition(error, partition, comm)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> New work partition
   type(work_partition), intent(out) :: partition

   !> MPI communicator handle, users of mpi_f08 pass comm%mpi_val
   integer, intent(in) :: comm

   integer :: rank, nranks

   call get_mpi_comm_info(comm, rank, nranks, error)
   if (allocated(error)) return

   call new_work_partition(error, partition, rank, nranks)

end subroutine new_mpi_work_partition


!> Calculate atom-resolved dispersion energies distributed over a communicator
subroutine get_dispersion_mpi_atomic(error, mol, disp, param, cutoff, comm, &
      & energies, gradient, sigma, hessian)

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

   !> MPI communicator handle, users of mpi_f08 pass comm%mpi_val
   integer, intent(in) :: comm

   !> Dispersion energy
   real(wp), intent(out), contiguous :: energies(:)

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   !> Dispersion hessian
   real(wp), intent(out), contiguous, optional :: hessian(:, :)

   type(work_partition) :: partition

   call new_mpi_work_partition(error, partition, comm)
   if (allocated(error)) return

   call get_dispersion(error, mol, disp, param, cutoff, energies, gradient, &
      & sigma, hessian, partition)
   if (allocated(error)) return

   call mpi_allreduce_sum(comm, energies, error)
   if (allocated(error)) return

   if (present(gradient)) then
      call mpi_allreduce_sum(comm, gradient, error)
      if (allocated(error)) return
   end if

   if (present(sigma)) then
      call mpi_allreduce_sum(comm, sigma, error)
      if (allocated(error)) return
   end if

   if (present(hessian)) then
      call mpi_allreduce_sum(comm, hessian, error)
      if (allocated(error)) return
   end if

end subroutine get_dispersion_mpi_atomic


!> Calculate scalar dispersion energy distributed over a communicator
subroutine get_dispersion_mpi_scalar(error, mol, disp, param, cutoff, comm, &
      & energy, gradient, sigma, hessian)

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

   !> MPI communicator handle, users of mpi_f08 pass comm%mpi_val
   integer, intent(in) :: comm

   !> Dispersion energy
   real(wp), intent(out) :: energy

   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   !> Dispersion hessian
   real(wp), intent(out), contiguous, optional :: hessian(:, :)

   real(wp), allocatable :: energies(:)

   energy = 0.0_wp
   allocate(energies(mol%nat))

   call get_dispersion_mpi_atomic(error, mol, disp, param, cutoff, comm, energies, &
      & gradient, sigma, hessian)
   if (allocated(error)) return

   energy = sum(energies)

end subroutine get_dispersion_mpi_scalar


end module dftd3_mpi
