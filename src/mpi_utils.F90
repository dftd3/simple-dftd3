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

!> Thin wrappers around the MPI calls used by s-dftd3.
!>
!> This is the only place besides dftd3_feature that is preprocessed. Without
!> MPI support every wrapper reports the missing feature in the error handler.
module dftd3_mpi_utils
#ifdef WITH_MPI
   use mpi_f08, only : MPI_Comm, MPI_Comm_rank, MPI_Comm_size, MPI_Initialized, &
      & MPI_Allreduce, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_IN_PLACE, MPI_SUCCESS
#endif
   use mctc_env, only : wp, error_type, fatal_error
   implicit none
   private

   public :: get_mpi_comm_info, mpi_allreduce_sum


   !> Sum a real quantity over all ranks of a communicator, in place
   interface mpi_allreduce_sum
      module procedure :: allreduce_sum_rank1
      module procedure :: allreduce_sum_rank2
   end interface mpi_allreduce_sum

   character(len=*), parameter :: unavailable = &
      & "s-dftd3 was built without MPI support"


contains


!> Determine rank and size of a communicator
subroutine get_mpi_comm_info(comm, rank, nranks, error)

   !> MPI communicator handle, users of mpi_f08 pass comm%mpi_val
   integer, intent(in) :: comm

   !> Zero-based rank of the calling process
   integer, intent(out) :: rank

   !> Number of processes in the communicator
   integer, intent(out) :: nranks

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

#ifdef WITH_MPI
   type(MPI_Comm) :: mpi_comm
   integer :: stat
   logical :: initialized

   rank = 0
   nranks = 1

   call MPI_Initialized(initialized, stat)
   if (stat /= MPI_SUCCESS .or. .not.initialized) then
      call fatal_error(error, "MPI is not initialized")
      return
   end if

   ! the integer handle is the portable currency between the language bindings
   mpi_comm%mpi_val = comm

   call MPI_Comm_rank(mpi_comm, rank, stat)
   if (stat /= MPI_SUCCESS) then
      call fatal_error(error, "Could not determine rank of MPI communicator")
      return
   end if

   call MPI_Comm_size(mpi_comm, nranks, stat)
   if (stat /= MPI_SUCCESS) then
      call fatal_error(error, "Could not determine size of MPI communicator")
      return
   end if
#else
   rank = 0
   nranks = 1
   call fatal_error(error, unavailable)
#endif

end subroutine get_mpi_comm_info


!> Sum an array over all ranks of a communicator
subroutine allreduce_sum_rank1(comm, val, error)

   !> MPI communicator handle
   integer, intent(in) :: comm

   !> Values to reduce
   real(wp), intent(inout), contiguous :: val(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

#ifdef WITH_MPI
   type(MPI_Comm) :: mpi_comm
   integer :: stat

   mpi_comm%mpi_val = comm
   call MPI_Allreduce(MPI_IN_PLACE, val, size(val), MPI_DOUBLE_PRECISION, MPI_SUM, &
      & mpi_comm, stat)
   if (stat /= MPI_SUCCESS) then
      call fatal_error(error, "Could not reduce results over MPI communicator")
      return
   end if
#else
   call fatal_error(error, unavailable)
#endif

end subroutine allreduce_sum_rank1


!> Sum a matrix over all ranks of a communicator
subroutine allreduce_sum_rank2(comm, val, error)

   !> MPI communicator handle
   integer, intent(in) :: comm

   !> Values to reduce
   real(wp), intent(inout), contiguous, target :: val(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), pointer :: ptr(:)

   ptr(1:size(val)) => val
   call allreduce_sum_rank1(comm, ptr, error)

end subroutine allreduce_sum_rank2


end module dftd3_mpi_utils
