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

!> Query the optional features s-dftd3 was built with
module dftd3_feature
   implicit none
   private

   public :: dftd3_has_mpi, dftd3_has_feature


   !> Whether distributed memory parallelisation with MPI is available
#ifdef WITH_MPI
   logical, parameter :: dftd3_has_mpi = .true.
#else
   logical, parameter :: dftd3_has_mpi = .false.
#endif


contains


!> Query an optional feature by name at runtime, unknown features are not available
pure function dftd3_has_feature(feature) result(has)

   !> Name of the feature, currently only "mpi" is defined
   character(len=*), intent(in) :: feature

   !> Whether the feature is available in this build
   logical :: has

   select case(feature)
   case default
      has = .false.
   case("mpi")
      has = dftd3_has_mpi
   end select

end function dftd3_has_feature


end module dftd3_feature
