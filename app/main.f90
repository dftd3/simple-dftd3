! This file is part of s-dftd3.
! SPDX-Identifier: LGLP-3.0-or-later
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

program dftd3_main
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env
   use mctc_io
   use dftd3
   use dftd3_output
   implicit none
   character(len=*), parameter :: prog_name = "s-dftd3"

   type :: d3_config
      logical :: json = .true.
      character(len=:), allocatable :: json_output
      logical :: tmer = .true.
      logical :: properties = .false.
      logical :: atm = .false.
      logical :: grad = .false.
      logical :: zero = .false.
      logical :: rational = .false.
      logical :: has_param = .false.
      integer :: verbosity = 2
   end type d3_config
   type(d3_config) :: config

   character(len=:), allocatable :: input
   integer, allocatable :: input_format
   type(structure_type) :: mol
   class(damping_param), allocatable :: param
   type(d3_model) :: d3
   type(d3_param) :: inp
   type(error_type), allocatable :: error
   real(wp) :: energy, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)
   character(len=:), allocatable :: method
   type(zero_damping_param), allocatable :: zparam
   type(rational_damping_param), allocatable :: rparam
   integer :: stat, unit
   logical :: echo, exist

   call get_arguments(input, input_format, config, method, inp, echo, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (input == "-") then
      if (.not.allocated(input_format)) input_format = filetype%xyz
      call read_structure(mol, input_unit, input_format, error)
   else
      call read_structure(mol, input, error, input_format)
   end if
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (config%zero) then
      allocate(zparam)
      call new_zero_damping(zparam, inp, mol%num)
      call move_alloc(zparam, param)
   end if
   if (config%rational) then
      allocate(rparam)
      call new_rational_damping(rparam, inp, mol%num)
      call move_alloc(rparam, param)
   end if

   if (config%grad) then
      allocate(gradient(3, mol%nat))
   end if

   call new_d3_model(d3, mol)

   if (config%properties) then
      call property_calc(output_unit, mol, d3, config%verbosity)
   end if

   if (allocated(param)) then
      call get_dispersion(mol, d3, param, realspace_cutoff(), energy, gradient, &
         & sigma)
      if (config%verbosity > 0) then
         call ascii_results(output_unit, mol, energy, gradient, sigma)
      end if
      if (config%tmer) then
         if (config%verbosity > 0) then
            write(output_unit, '(a)') "[Info] Dispersion energy written to .EDISP"
         end if
         open(file=".EDISP", newunit=unit)
         write(unit, '(f24.14)') energy
         close(unit)
      end if
      if (config%grad) then
         inquire(file="gradient", exist=exist)
         if (exist) then
            call turbomole_gradient(mol, "gradient", energy, gradient, stat)
            if (stat == 0) then
               write(output_unit, '(a)') &
                  & "[Info] Dispersion gradient added to Turbomole gradient file"
            else
               write(output_unit, '(a)') &
                  & "[Warn] Could not add to Turbomole gradient file"
            end if
         end if
         inquire(file="gradlatt", exist=exist)
         if (exist) then
            call turbomole_gradlatt(mol, "gradlatt", energy, sigma, stat)
            if (stat == 0) then
               write(output_unit, '(a)') &
                  & "[Info] Dispersion virial added to Turbomole gradlatt file"
            else
               write(output_unit, '(a)') &
                  & "[Warn] Could not add to Turbomole gradlatt file"
            end if
         end if
      end if

      if (config%json) then
         open(file=config%json_output, newunit=unit)
         if (config%grad) then
            call json_results(unit, "  ", energy, gradient, sigma)
         else
            call json_results(unit, "  ", energy)
         end if
         close(unit)
         write(output_unit, '(a)') &
            & "[Info] JSON dump of results written to '"//config%json_output//"'"
      end if
   end if



contains


subroutine property_calc(unit, mol, disp, verbosity)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Printout verbosity
   integer, intent(in) :: verbosity

   integer :: mref
   real(wp), allocatable :: cn(:), gwvec(:, :), c6(:, :)
   real(wp), parameter :: lattr(3, 1) = 0.0_wp

   if (verbosity > 1) then
      call ascii_atomic_radii(unit, mol, disp)
      write(unit, '(a)')
      call ascii_atomic_references(unit, mol, disp)
      write(unit, '(a)')
   end if

   mref = maxval(disp%ref)
   allocate(cn(mol%nat), gwvec(mref, mol%nat), c6(mol%nat, mol%nat))
   call get_coordination_number(mol, lattr, 30.0_wp, disp%rcov, cn)
   call disp%weight_references(mol, cn, gwvec)
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   if (verbosity > 0) then
      call ascii_system_properties(unit, mol, disp, cn, c6)
      write(unit, '(a)')
   end if

end subroutine property_calc


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      "Takes an geometry input to calculate the D3 dispersion correction.", &
      "Periodic calculations are performed automatically for periodic input formats", &
      "Specify the functional to select the correct parameters.", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "-i, --input <format>", "Hint for the format of the input file", &
      "--func <method>", "Specify functional to use", &
      "--bj", "Use rational (Becke-Johnson) damping function", &
      "--bj-param <list>", "Specify parameters for rational damping,", &
      "", "expected order is s6, s8, a1, a2 (requires four arguments)", &
      "--zero", "Use zero (Chai-Head-Gordon) damping function", &
      "--zero-param <list>", "Specify functional to use,", &
      "", "expected order is s6, s8, rs6 (requires three arguments)", &
      "--atm", "Use ATM three-body dispersion,", &
      "--atm-scale <s9>", "Use scaled ATM three-body dispersion,", &
      "--noedisp", "Disable writing of dispersion energy to .EDISP file", &
      "--json [file]", "Dump results to JSON output (default: dftd3.json)", &
      "--grad", "Request gradient evaluation", &
      "--property", "Evaluate dispersion related properties", &
      "-v, --verbose", "Show more, can be used multiple times", &
      "-s, --silent", "Show less, use twice to supress all output", &
      "--version", "Print program version and exit", &
      "--help", "Show this help message"

   write(unit, '(a)')

end subroutine help


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_dftd3_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


!> Obtain the command line argument at a given index
subroutine get_argument(idx, arg)

   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: idx

   !> Command line argument
   character(len=:), allocatable, intent(out) :: arg

   integer :: length, stat

   call get_command_argument(idx, length=length, status=stat)
   if (stat /= 0) then
      return
   endif

   allocate(character(len=length) :: arg, stat=stat)
   if (stat /= 0) then
      return
   endif

   if (length > 0) then
      call get_command_argument(idx, arg, status=stat)
      if (stat /= 0) then
         deallocate(arg)
         return
      end if
   end if

end subroutine get_argument


subroutine get_argument_as_real(iarg, val, error)

   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: iarg

   !> Real value
   real(wp), intent(out) :: val

   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat
   character(len=:), allocatable :: arg

   call get_argument(iarg, arg)
   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read real value, argument missing")
      return
   end if
   read(arg, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read real value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_real


subroutine get_arguments(input, input_format, config, method, inp, echo, error)

   character(len=:), allocatable, intent(out) :: input
   integer, allocatable, intent(out) :: input_format

   !> Configuation data
   type(d3_config), intent(out) :: config

   !> Method name
   character(len=:), allocatable, intent(out) :: method

   !> Damping paramaters
   type(d3_param), intent(out) :: inp

   !> Print information
   logical, intent(out) :: echo

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   character(len=:), allocatable :: arg

   echo = .true.
   iarg = 0
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      select case(arg)
      case("--help")
         call help(output_unit)
         stop
      case("--version")
         call version(output_unit)
         stop
      case("-v", "--verbose")
         config%verbosity = config%verbosity + 1
      case("-s", "--silent")
         config%verbosity = config%verbosity - 1
      case default
         if (.not.allocated(input)) then
            call move_alloc(arg, input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case("-i", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         input_format = get_filetype("."//arg)
      case("--json")
         config%json = .true.
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%json_output)
         else
            config%json_output = "dftd3.json"
         end if
      case("--property")
         config%properties = .true.
      case("--noedisp")
         config%tmer = .false.
      case("--grad")
         config%grad = .true.
      case("--atm")
         inp%s9 = 1.0_wp
         config%atm = .true.
      case("--atm-scale")
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%s9, error)
         if (allocated(error)) exit
         config%atm = .true.
      case("--zero-param")
         config%zero = .true.
         config%has_param = .true.
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%s6, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%s8, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%rs6, error)
         if (allocated(error)) exit
      case("--bj-param")
         config%rational = .true.
         config%has_param = .true.
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%s6, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%s8, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%a1, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%a2, error)
         if (allocated(error)) exit
      case("--func")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, method)
      end select
   end do

   if (.not.config%has_param .and. .not.allocated(method)) then
      config%properties = .true.
   end if

   if (config%zero .and. config%rational) then
      call fatal_error(error, "Can only select zero or rational damping function")
      return
   end if

   if (.not.allocated(input)) then
      if (.not.allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_arguments


end program dftd3_main
