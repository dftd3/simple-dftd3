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

module d3def_environment
   use iso_fortran_env
   implicit none
   public :: d3_environment
   private


   type :: d3_error
      integer :: error_code = 0
      character(len=:), allocatable :: message
      character(len=:), allocatable :: origin
   end type d3_error

   type :: d3_environment
      private
      integer, public :: unit = output_unit
      integer :: error_code = 0
      integer :: error_count = 0
      type(d3_error), allocatable :: message(:)
   contains
      generic :: new => new_default
      generic :: reset => new_default
      procedure, private :: new_default => env_new_default
      procedure :: get_error => env_get_error
      procedure :: get_finalize => env_get_finalize
      procedure :: set_error => env_set_error
      procedure :: set_finalize => env_set_finalize
      procedure :: checkpoint => env_checkpoint
      procedure, private :: push_back => env_push_back
      procedure, private :: resize => env_resize
      procedure, private :: terminate => env_terminate
      procedure :: destroy => env_destroy
      final :: env_finalizer
   end type d3_environment


   interface len
      module procedure :: env_length
   end interface

   interface size
      module procedure :: env_size
   end interface

contains


subroutine env_new_default(self)
   class(d3_environment), intent(out) :: self
end subroutine env_new_default

integer pure elemental function env_length(self) result(length)
   class(d3_environment), intent(in) :: self
   if (allocated(self%message)) then
      length = self%error_count
   else
      length = 0
   endif
end function env_length

integer pure elemental function env_size(self) result(length)
   class(d3_environment), intent(in) :: self
   if (allocated(self%message)) then
      length = size(self%message)
   else
      length = 0
   endif
end function env_size


logical pure elemental function env_get_error(self) result(error)
   class(d3_environment), intent(in) :: self
   error = self%error_code > 0
end function env_get_error

logical pure elemental function env_get_finalize(self) result(error)
   class(d3_environment), intent(in) :: self
   error = self%error_code /= 0
end function env_get_finalize

subroutine env_set_error(self, message, origin, error_code, terminate)
   class(d3_environment), intent(inout) :: self
   character(len=*), intent(in) :: message
   character(len=*), intent(in), optional :: origin
   integer, intent(in), optional :: error_code
   logical, intent(in), optional :: terminate
   if (present(error_code)) then
      self%error_code = error_code
   else
      self%error_code = 128
   endif
   call self%push_back(message)
   if (present(terminate)) then
      if (terminate) call self%checkpoint("forced termination")
   endif
end subroutine env_set_error

subroutine env_set_finalize(self, message, origin, error_code, terminate)
   class(d3_environment), intent(inout) :: self
   character(len=*), intent(in) :: message
   character(len=*), intent(in), optional :: origin
   integer, intent(in), optional :: error_code
   logical, intent(in), optional :: terminate
   if (present(error_code)) then
      self%error_code = error_code
   else
      self%error_code = -1
   endif
   call self%push_back(message)
   if (present(terminate)) then
      if (terminate) call self%checkpoint("forced termination")
   endif
end subroutine env_set_finalize

subroutine env_checkpoint(self, message)
   class(d3_environment), intent(inout) :: self
   character(len=*), intent(in) :: message
   integer :: ierr
   if (self%get_error()) then
      write(self%unit, '(a)')
      write(self%unit, '("ERROR !>",1x,a,1x,"<!")') message
      do ierr = 1, len(self)
         write(self%unit, '("-->",1x,a)') self%message(ierr)%message
      enddo
      write(self%unit, '(a)')
   endif
   if (self%get_finalize()) call self%terminate
end subroutine env_checkpoint


pure subroutine env_push_back(self, message, origin, error_code)
   class(d3_environment), intent(inout) :: self
   character(len=*), intent(in) :: message
   character(len=*), intent(in), optional :: origin
   integer, intent(in), optional :: error_code
   integer :: pos, last, length
   last = len(self)
   if (last >= size(self)) call self%resize
   self%error_count = last+1
   self%message(last+1)%message = message
end subroutine env_push_back

pure subroutine env_resize(self, n)
   class(d3_environment), intent(inout) :: self
   integer, intent(in), optional :: n
   type(d3_error), allocatable :: message(:)
   integer :: length, current_length
   current_length = size(self)
   if (current_length > 0) then
      if (present(n)) then
         if (n <= current_length) return
         length = n
      else
         length = current_length + current_length/2 + 1
      endif
      allocate(message(length))
      message(:current_length) = self%message(:current_length)
      deallocate(self%message)
      call move_alloc(message, self%message)
   else
      if (present(n)) then
         length = n
      else
         length = 64
      endif
      allocate(self%message(length))
   endif
end subroutine env_resize

subroutine env_terminate(self)
   class(d3_environment), intent(inout) :: self
   if (self%get_error()) then
      error stop "abnormal termination of s-dftd3"
   else
      stop "normal termination of s-dftd3"
   endif
end subroutine env_terminate


pure subroutine env_destroy(self)
   class(d3_environment), intent(inout) :: self
end subroutine env_destroy

pure subroutine env_finalizer(self)
   type(d3_environment), intent(inout) :: self
   call self%destroy
end subroutine env_finalizer


end module d3def_environment

