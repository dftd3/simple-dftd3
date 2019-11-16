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

module d3def_argparser
   use iso_fortran_env
   implicit none
   public :: d3_argparser
   public :: size, len, file, write(formatted)
   private


   type :: d3_argument
      character(len=:), allocatable :: string
      logical :: option = .false.
      logical :: unread = .true.
      logical :: isfile = .false.
   contains
      procedure :: new => arg_new
      procedure :: to_string => arg_to_string
      procedure :: get_value => arg_get_value
      procedure, private :: to_int16 => arg_to_int16
      procedure, private :: to_int32 => arg_to_int32
      procedure, private :: to_int64 => arg_to_int64
      procedure, private :: to_real32 => arg_to_real32
      procedure, private :: to_real64 => arg_to_real64
      procedure, private :: to_logical => arg_to_logical
      procedure :: destroy => arg_destroy
   end type d3_argument

   type :: d3_argparser
      type(d3_argument), allocatable :: args(:)
   contains
      generic :: new => new_default
      generic :: reset => new_default
      procedure, private :: new_default => ap_new_default
      procedure :: has_option => ap_has_option
      procedure :: get_argument => ap_get_argument
      procedure :: get_option => ap_get_option
      procedure :: get_file => ap_get_filename
      procedure :: destroy => ap_destroy
      final :: ap_finalizer
   end type d3_argparser


   interface size
      module procedure ap_get_size
   end interface size

   interface file
      module procedure ap_get_filecount
   end interface file

   interface len
      module procedure ap_get_length
   end interface len

   interface write(formatted)
      module procedure ap_write_formatted
   end interface write(formatted)


contains


subroutine ap_new_default(self)
   class(d3_argparser), intent(out) :: self
   integer :: iarg, nargs
   logical :: getopts
   nargs = command_argument_count()
   getopts = .true.
   allocate(self%args(nargs))
   do iarg = 1, nargs
      call self%args(iarg)%new(iarg, getopts)
   enddo
end subroutine ap_new_default

subroutine arg_new(self, iarg, getopts)
   class(d3_argument), intent(out) :: self
   integer, intent(in) :: iarg
   logical, intent(inout) :: getopts
   integer :: length, iopt
   ! get command line argument
   call get_command_argument(iarg, length=length)
   allocate(character(len=length) :: self%string)
   call get_command_argument(iarg, self%string)
   ! check if the argument provided is a file
   inquire(file=self%string, exist=self%isfile)
   ! check if we got two dashes
   if (self%string == '--' .and. getopts) then
      getopts = .false.
      self%isfile = .false.
   endif
   if (getopts) then
      iopt = index(self%string, '--')
      self%option = iopt == 1
   endif
end subroutine arg_new


logical function ap_has_option(self, option) result(found)
   class(d3_argparser), intent(inout) :: self
   character(len=*), intent(in) :: option
   integer :: iarg
   found = .false.
   do iarg = 1, size(self)
      associate( arg => self%args(iarg) )
         if (arg%option .and. arg%unread) then
            found = arg%string == '--'//option
            if (found) then
               arg%unread = .false.
               exit
            endif
         endif
      end associate
   enddo
end function ap_has_option

subroutine ap_get_argument(self, argument)
   class(d3_argparser), intent(inout) :: self
   character(len=:), allocatable, intent(out) :: argument
   integer :: iarg, error
   logical :: found_local
   found_local = .false.
   do iarg = 1, size(self)
      associate( arg => self%args(iarg) )
         if (arg%unread) then
            arg%unread = .false.
            argument = arg%string
            exit
         endif
      end associate
   enddo
end subroutine ap_get_argument

subroutine ap_get_option(self, option, value, found)
   class(d3_argparser), intent(inout) :: self
   character(len=*), intent(in) :: option
   class(*), intent(out) :: value
   logical, intent(out), optional :: found
   integer :: iarg, error
   logical :: found_local
   found_local = .false.
   do iarg = 1, size(self)-1
      associate( arg => self%args(iarg), sec => self%args(iarg+1) )
         if (arg%option .and. arg%unread .and. sec%unread) then
            found_local = arg%string == '--'//option
            if (found_local) then
               call sec%get_value(value, error)
               if (error /= 0) then
                  arg%unread = .false.
                  sec%unread = .false.
               else
                  found_local = .false.
               endif
               exit
            endif
         endif
      end associate
   enddo
   if (present(found)) found = found_local
end subroutine ap_get_option

subroutine ap_get_filename(self, filename)
   class(d3_argparser), intent(inout) :: self
   character(len=:), allocatable, intent(out) :: filename
   integer :: iarg
   do iarg = 1, size(self)
      associate( arg => self%args(iarg) )
         if (arg%isfile .and. arg%unread) then
            arg%unread = .false.
            filename = arg%string
            exit
         endif
      end associate
   enddo
end subroutine ap_get_filename


integer pure elemental function ap_get_size(self) result(length)
   class(d3_argparser), intent(in) :: self
   if (allocated(self%args)) then
      length = size(self%args)
   else
      length = 0
   endif
end function ap_get_size

integer pure elemental function ap_get_length(self) result(length)
   class(d3_argparser), intent(in) :: self
   if (allocated(self%args)) then
      length = count(self%args%unread)
   else
      length = 0
   endif
end function ap_get_length

integer pure elemental function ap_get_filecount(self) result(length)
   class(d3_argparser), intent(in) :: self
   if (allocated(self%args)) then
      length = count(self%args%unread .and. self%args%isfile)
   else
      length = 0
   endif
end function ap_get_filecount

subroutine ap_write_formatted(self, unit, iotype, v_list, iostat, iomsg)
   class(d3_argparser), intent(in) :: self
   integer, intent(in) :: unit
   character(len=*), intent(in) :: iotype
   integer, intent(in) :: v_list(:)
   integer, intent(out) :: iostat
   character(len=*), intent(inout) :: iomsg
   character(len=:), allocatable :: string
   integer :: iarg
   do iarg = 1, size(self)
      call self%args(iarg)%to_string(string)
      write(unit, '(i0,":",1x,a)') iarg, string
   enddo
end subroutine ap_write_formatted

subroutine arg_to_string(self, string)
   class(d3_argument), intent(in) :: self
   character(len=:), allocatable, intent(out) :: string
   string = self%string
   if (self%isfile) string = string // ' (file exist)'
   if (self%option) string = string // ' (is option)'
end subroutine arg_to_string

subroutine arg_get_value(self, value, error)
   class(d3_argument), intent(in) :: self
   class(*), intent(out) :: value
   integer, intent(out), optional :: error
   integer :: stat
   select type(value)
   type is (real(real32)); call self%to_real32(value, stat)
   type is (real(real64)); call self%to_real64(value, stat)
   type is (integer(int16)); call self%to_int16(value, stat)
   type is (integer(int32)); call self%to_int32(value, stat)
   type is (integer(int64)); call self%to_int64(value, stat)
   type is (logical); call self%to_logical(value, stat)
   class default
      stat = 128
   end select
   if (present(error)) error = stat
end subroutine arg_get_value

subroutine arg_to_int16(self, value, error)
   class(d3_argument), intent(in) :: self
   integer(int16), intent(out) :: value
   integer, intent(out), optional :: error
   integer :: stat
   read(self%string, *, iostat=stat) value
   if (present(error)) error = stat
end subroutine arg_to_int16

subroutine arg_to_int32(self, value, error)
   class(d3_argument), intent(in) :: self
   integer(int32), intent(out) :: value
   integer, intent(out), optional :: error
   integer :: stat
   read(self%string, *, iostat=stat) value
   if (present(error)) error = stat
end subroutine arg_to_int32

subroutine arg_to_int64(self, value, error)
   class(d3_argument), intent(in) :: self
   integer(int64), intent(out) :: value
   integer, intent(out), optional :: error
   integer :: stat
   read(self%string, *, iostat=stat) value
   if (present(error)) error = stat
end subroutine arg_to_int64

subroutine arg_to_real32(self, value, error)
   class(d3_argument), intent(in) :: self
   real(real32), intent(out) :: value
   integer, intent(out), optional :: error
   integer :: stat
   read(self%string, *, iostat=stat) value
   if (present(error)) error = stat
end subroutine arg_to_real32

subroutine arg_to_real64(self, value, error)
   class(d3_argument), intent(in) :: self
   real(real64), intent(out) :: value
   integer, intent(out), optional :: error
   integer :: stat
   read(self%string, *, iostat=stat) value
   if (present(error)) error = stat
end subroutine arg_to_real64

subroutine arg_to_logical(self, value, error)
   class(d3_argument), intent(in) :: self
   logical, intent(out) :: value
   integer, intent(out), optional :: error
   integer :: stat
   select case(self%string)
   case('Y','y','Yes','yes','T','t','true','True','1')
      stat = .true.
      value = .true.
   case('N','n','No','no','F','f','false','False','0')
      stat = .true.
      value = .false.
   case default
      stat = .false.
   end select
   if (present(error)) error = stat
end subroutine arg_to_logical


subroutine arg_destroy(self)
   class(d3_argument), intent(inout) :: self
   if (allocated(self%string)) deallocate(self%string)
end subroutine arg_destroy

subroutine ap_destroy(self)
   class(d3_argparser), intent(inout) :: self
   integer :: iarg
   do iarg = 1, size(self)
      call self%args(iarg)%destroy
   enddo
   if (allocated(self%args)) deallocate(self%args)
end subroutine ap_destroy

subroutine ap_finalizer(self)
   type(d3_argparser), intent(inout) :: self
   call self%destroy
end subroutine ap_finalizer


end module d3def_argparser
