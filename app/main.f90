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

program dftd3_main
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type, read_structure, filetype, get_filetype
   use dftd3
   use dftd3_output
   use dftd3_utils
   implicit none
   character(len=*), parameter :: prog_name = "s-dftd3"

   type :: d3_config
      logical :: json = .false.
      character(len=:), allocatable :: json_output
      logical :: wrap = .true.
      logical :: tmer = .true.
      logical :: properties = .false.
      logical :: atm = .false.
      logical :: grad = .false.
      character(len=:), allocatable :: grad_output
      logical :: zero = .false.
      logical :: rational = .false.
      logical :: mzero = .false.
      logical :: mrational = .false.
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
   type(mzero_damping_param), allocatable :: mparam
   type(rational_damping_param), allocatable :: rparam
   real(wp), allocatable :: s9
   integer :: stat, unit
   logical :: echo, exist

   call get_arguments(input, input_format, config, method, inp, echo, error)
   if (allocated(error)) then
      write(error_unit, '("[Error]", 1x, a)') error%message
      error stop
   end if

   if (config%verbosity > 1) then
      call header(output_unit)
   end if

   if (input == "-") then
      if (.not.allocated(input_format)) input_format = filetype%xyz
      call read_structure(mol, input_unit, input_format, error)
   else
      call read_structure(mol, input, error, input_format)
   end if
   if (allocated(error)) then
      write(error_unit, '("[Error]", 1x, a)') error%message
      error stop
   end if
   if (config%wrap) then
      call wrap_to_central_cell(mol%xyz, mol%lattice, mol%periodic)
   end if

   if (config%atm) s9 = inp%s9
   if (config%zero) then
      if (.not.config%has_param) then
         call get_zero_damping(inp, method, error, s9)
         if (allocated(error)) then
            write(error_unit, '("[Error]", 1x, a)') error%message
            error stop
         end if
      end if
      allocate(zparam)
      call new_zero_damping(zparam, inp, mol%num)
      call move_alloc(zparam, param)
   end if
   if (config%mzero) then
      if (.not.config%has_param) then
         call get_mzero_damping(inp, method, error, s9)
         if (allocated(error)) then
            write(error_unit, '("[Error]", 1x, a)') error%message
            error stop
         end if
      end if
      allocate(mparam)
      call new_mzero_damping(mparam, inp, mol%num)
      call move_alloc(mparam, param)
   end if
   if (config%rational .or. config%mrational) then
      if (.not.config%has_param) then
         if (config%mrational) then
            call get_mrational_damping(inp, method, error, s9)
         else
            call get_rational_damping(inp, method, error, s9)
         end if
         if (allocated(error)) then
            write(error_unit, '("[Error]", 1x, a)') error%message
            error stop
         end if
      end if
      allocate(rparam)
      call new_rational_damping(rparam, inp, mol%num)
      call move_alloc(rparam, param)
   end if

   if (allocated(param) .and. config%verbosity > 0) then
      call ascii_damping_param(output_unit, param, method)
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
         open(file=config%grad_output, newunit=unit)
         call tagged_result(unit, energy, gradient, sigma)
         close(unit)
         if (config%verbosity > 0) then
            write(output_unit, '(a)') &
               & "[Info] Dispersion results written to '"//config%grad_output//"'"
         end if

         inquire(file="gradient", exist=exist)
         if (exist) then
            call turbomole_gradient(mol, "gradient", energy, gradient, stat)
            if (config%verbosity > 0) then
               if (stat == 0) then
                  write(output_unit, '(a)') &
                     & "[Info] Dispersion gradient added to Turbomole gradient file"
               else
                  write(output_unit, '(a)') &
                     & "[Warn] Could not add to Turbomole gradient file"
               end if
            end if
         end if
         inquire(file="gradlatt", exist=exist)
         if (exist) then
            call turbomole_gradlatt(mol, "gradlatt", energy, sigma, stat)
            if (config%verbosity > 0) then
               if (stat == 0) then
                  write(output_unit, '(a)') &
                     & "[Info] Dispersion virial added to Turbomole gradlatt file"
               else
                  write(output_unit, '(a)') &
                     & "[Warn] Could not add to Turbomole gradlatt file"
               end if
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
         if (config%verbosity > 0) then
            write(output_unit, '(a)') &
               & "[Info] JSON dump of results written to '"//config%json_output//"'"
         end if
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
   real(wp), allocatable :: cn(:), gwvec(:, :), c6(:, :), lattr(:, :)

   if (verbosity > 1) then
      call ascii_atomic_radii(unit, mol, disp)
      write(unit, '(a)')
      call ascii_atomic_references(unit, mol, disp)
      write(unit, '(a)')
   end if

   mref = maxval(disp%ref)
   allocate(cn(mol%nat), gwvec(mref, mol%nat), c6(mol%nat, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, 30.0_wp, lattr)
   call get_coordination_number(mol, lattr, 30.0_wp, disp%rcov, cn)
   call disp%weight_references(mol, cn, gwvec)
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   if (verbosity > 0) then
      call ascii_system_properties(unit, mol, disp, cn, c6)
      write(unit, '(a)')
   end if

end subroutine property_calc


subroutine header(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_dftd3_version(string=version_string)
   write(unit, '(a)') &
      "-----------------------------------", &
      " s i m p l e   D F T - D 3  v"// version_string, &
      "-----------------------------------", ""

end subroutine header


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      "Takes an geometry input to calculate the D3 dispersion correction.", &
      "Periodic calculations are performed automatically for periodic input formats.", &
      "Specify the functional to select the correct parameters.", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "-i, --input <format>", "Hint for the format of the input file", &
      "--bj <method>", "Use rational (Becke-Johnson) damping function", &
      "--bj-param <list>", "Specify parameters for rational damping,", &
      "", "expected order is s6, s8, a1, a2 (requires four arguments)", &
      "--zero <method>", "Use zero (Chai-Head-Gordon) damping function", &
      "--zero-param <list>", "Specify parameters for zero damping,", &
      "", "expected order is s6, s8, rs6 (requires three arguments)", &
      "--bjm <method>", "Use modified rational damping function", &
      "--bjm-param <list>", "Specify parameters for rational damping,", &
      "", "expected order is s6, s8, a1, a2 (requires four arguments)", &
      "--zerom <method>", "Use modified zero damping function", &
      "--zerom-param <list>", "Specify parameters for modified zero damping,", &
      "", "expected order is s6, s8, rs6, bet (requires four arguments)", &
      "--atm", "Use ATM three-body dispersion", &
      "--atm-scale <s9>", "Use scaled ATM three-body dispersion", &
      "--noedisp", "Disable writing of dispersion energy to .EDISP file", &
      "--json [file]", "Dump results to JSON output (default: dftd3.json)", &
      "--grad [file]", "Request gradient evaluation,", &
      "", "write results to file (default: dftd3.txt),", &
      "", "attempts to add to Turbomole gradient and gradlatt files", &
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
         if (arg(1:1) == "-") then
            call fatal_error(error, "Unknown argument encountered: '"//arg//"'")
         else
            call fatal_error(error, "Too many positional arguments present")
         end if
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
         config%json_output = "dftd3.json"
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%json_output)
         end if
      case("--property")
         config%properties = .true.
      case("--noedisp")
         config%tmer = .false.
      case("--nowrap")
         config%wrap = .false.
      case("--grad")
         config%grad = .true.
         config%grad_output = "dftd3.txt"
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%grad_output)
         end if
      case("--atm")
         inp%s9 = 1.0_wp
         config%atm = .true.
      case("--atm-scale")
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%s9, error)
         if (allocated(error)) exit
         config%atm = .true.
      case("--zero")
         config%zero = .true.
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, method)
      case("--zerom")
         config%mzero = .true.
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, method)
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
      case("--zerom-param")
         config%mzero = .true.
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
         iarg = iarg + 1
         call get_argument_as_real(iarg, inp%bet, error)
         if (allocated(error)) exit
      case("--bj")
         config%rational = .true.
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, method)
      case("--bjm")
         config%mrational = .true.
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, method)
      case("--bj-param", "--bjm-param")
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
      end select
   end do
   if (allocated(error)) return

   if (.not.config%has_param .and. .not.allocated(method)) then
      config%properties = .true.
   end if

   if (count([config%zero, config%rational, config%mzero, config%mrational]) > 1) then
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
