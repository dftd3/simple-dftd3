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

program sdftd3_prog
   use d3def_environment
   use d3def_molecule
   use d3def_options
   use d3def_results
   use d3mod_header
   use d3mod_main
   use d3mod_utils_filetype
   implicit none
   type(d3_environment) :: env
   type(d3_options) :: opt
   type(d3_molecule) :: mol
   type(d3_results) :: res
   character(len=:), allocatable :: filename
   logical :: exist
   integer :: ifile

   call env%new

   call read_command_line_arguments(env, filename, opt)
   call env%checkpoint("reading command line arguments")

   call sdftd3_header(env%unit)

   inquire(file=filename, exist=exist)
   if (exist) then
      open(file=filename, newunit=ifile)
      call mol%read(ifile, format=p_ftype%default)
      if (mol%ftype < 0) call env%set_error(mol%name)
      close(ifile)
   else
      call env%set_error("input file not found")
   endif
   call env%checkpoint("reading molecule from '"//filename//"'")

   call main_run(env, opt, mol, res)
   call env%checkpoint("running calculation for '"//filename//"'")

contains

subroutine read_command_line_arguments(env, filename, opt)
   use d3def_argparser
   use d3mod_help
   use d3par_damping_parameters
   implicit none
   class(d3_environment), intent(inout) :: env
   class(d3_options), intent(inout) :: opt
   character(len=:), allocatable, intent(out) :: filename
   character(len=:), allocatable :: string

   type(d3_argparser) :: args
   integer :: iarg
   logical :: found

   call args%new

   if (size(args) == 0) then
      call helpmsg(env%unit)
      call env%set_error("no arguments given, so there is nothing to do")
      return
   endif

   ! check for quick help exit first
   if (args%has_option('help')) then
      call helpmsg(env%unit)
      call env%set_finalize("help message printed")
      return
   endif

   ! check for requested version number
   if (args%has_option('version')) then
      call sdftd3_version(env%unit)
      call env%set_finalize("version number printed")
      return
   endif

   ! look up the number of possible input files, we need at least one
   if (file(args) == 0) then
      call helpmsg(env%unit)
      call env%set_error("no file names given, so there is nothing to do")
      return
   endif

   ! now read the arguments
   call args%get_option('func', opt%func, found)
   if (found) then
      allocate(opt%par, source=d3_damping_parameters())
      call get_d3_damping_parameters(opt%par, opt%func, found)
      if (.not.found) then
         call env%set_error("functional '"//opt%func//"' not known")
      endif
   else
      if (allocated(opt%func)) deallocate(opt%func)
   endif

   ! check for damping parameters, if we are in energy/gradient mode
   if (.not.allocated(opt%par)) then
      call env%set_error("No damping parameters given, so there is nothing to do")
   endif

   ! everything left should be a file now, but we only want one input file
   if (file(args) > 1) then
      call env%set_error("multiple file names given, cannot process more than one")
      return
   endif
   call args%get_file(filename)

   ! at this point all command line arguments must be consumed, otherwise...
   if (len(args) > 0) then
      do iarg = 1, len(args)
         call args%get_argument(string)
         call env%set_error("unknown argument: '"//string//"'")
      enddo
   endif

end subroutine read_command_line_arguments

end program sdftd3_prog
