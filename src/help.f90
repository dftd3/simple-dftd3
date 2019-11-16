! This file is part of s-dftd3.
!
! Copyright (C) 2019 Sebastian Ehlert
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module d3mod_help
   implicit none
   public :: helpmsg

contains

subroutine helpmsg(unit)
   integer, intent(in) :: unit
   include 's-dftd3-version.fh'
   write(unit, '(a,1x,a,1x,a)') &
      "Usage:", name, "[options] <input>"
end subroutine helpmsg

end module d3mod_help
