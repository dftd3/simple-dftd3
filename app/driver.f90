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

module dftd3_app_driver
   use, intrinsic :: iso_fortran_env, only : output_unit, input_unit
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type, read_structure, filetype, get_filetype
   use dftd3
   use dftd3_output
   use dftd3_utils
   use dftd3_app_help, only : header
   use dftd3_app_cli, only : app_config, run_config, param_config, get_arguments
   use dftd3_app_toml, only : param_database
   implicit none
   private

   public :: app_driver

contains

subroutine app_driver(config, error)
   class(app_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   select type(config)
   type is(run_config)
      call run_driver(config, error)
   type is(param_config)
      call param_driver(config, error)
   end select
end subroutine app_driver

subroutine run_driver(config, error)
   type(run_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(damping_param), allocatable :: param
   type(d3_param) :: inp
   type(d3_model) :: d3
   real(wp), allocatable :: energies(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: pair_disp2(:, :), pair_disp3(:, :)
   real(wp), allocatable :: s9
   real(wp) :: energy
   integer :: stat, unit
   logical :: exist

   if (config%verbosity > 1) then
      call header(output_unit)
   end if

   if (config%input == "-") then
      if (.not.allocated(config%input_format)) then
         call read_structure(mol, input_unit, filetype%xyz, error)
      else
         call read_structure(mol, input_unit, config%input_format, error)
      end if
   else
      call read_structure(mol, config%input, error, config%input_format)
   end if
   if (allocated(error)) return
   if (config%wrap) then
      call wrap_to_central_cell(mol%xyz, mol%lattice, mol%periodic)
   end if

   if (config%has_param) inp = config%inp
   if (config%atm) s9 = config%inp%s9
   if (config%zero) then
      if (.not.config%has_param) then
         if (allocated(config%db)) then
            call from_db(param, config%db, config%method, "zero", error)
         else
            call get_zero_damping(inp, config%method, error, s9)
         end if
         if (allocated(error)) return
      end if
      if (.not.allocated(param)) then
         block
            type(zero_damping_param), allocatable :: zparam
            allocate(zparam)
            call new_zero_damping(zparam, inp)
            call move_alloc(zparam, param)
         end block
      end if
   end if
   if (config%mzero) then
      if (.not.config%has_param) then
         if (allocated(config%db)) then
            call from_db(param, config%db, config%method, "zerom", error)
         else
            call get_mzero_damping(inp, config%method, error, s9)
         end if
         if (allocated(error)) return
      end if
      if (.not.allocated(param)) then
         block
            type(mzero_damping_param), allocatable :: mparam
            allocate(mparam)
            call new_mzero_damping(mparam, inp)
            call move_alloc(mparam, param)
         end block
      end if
   end if
   if (config%rational .or. config%mrational) then
      if (.not.config%has_param) then
         if (config%mrational) then
            if (allocated(config%db)) then
               call from_db(param, config%db, config%method, "bjm", error)
            else
               call get_mrational_damping(inp, config%method, error, s9)
            end if
         else
            if (allocated(config%db)) then
               call from_db(param, config%db, config%method, "bj", error)
            else
               call get_rational_damping(inp, config%method, error, s9)
            end if
         end if
         if (allocated(error)) return
      end if
      if (.not.allocated(param)) then
         block
            type(rational_damping_param), allocatable :: rparam
            allocate(rparam)
            call new_rational_damping(rparam, inp)
            call move_alloc(rparam, param)
         end block
      end if
   end if
   if (config%optimizedpower) then
      if (.not.config%has_param) then
         if (allocated(config%db)) then
            call from_db(param, config%db, config%method, "op", error)
         else
            call get_optimizedpower_damping(inp, config%method, error, s9)
         end if
         if (allocated(error)) return
      end if
      if (.not.allocated(param)) then
         block
            type(optimizedpower_damping_param), allocatable :: oparam
            allocate(oparam)
            call new_optimizedpower_damping(oparam, inp)
            call move_alloc(oparam, param)
         end block
      end if
   end if

   if (allocated(param) .and. config%verbosity > 0) then
      call ascii_damping_param(output_unit, param, config%method)
   end if

   if (allocated(param)) then
      allocate(energies(mol%nat))
      if (config%grad) then
         allocate(gradient(3, mol%nat), sigma(3, 3))
      end if
   end if

   call new_d3_model(d3, mol)

   if (config%properties) then
      call property_calc(output_unit, mol, d3, config%verbosity)
   end if

   if (allocated(param)) then
      call get_dispersion(mol, d3, param, realspace_cutoff(), energies, &
         & gradient, sigma)
      energy = sum(energies)

      if (config%pair_resolved) then
         allocate(pair_disp2(mol%nat, mol%nat), pair_disp3(mol%nat, mol%nat))
         call get_pairwise_dispersion(mol, d3, param, realspace_cutoff(), pair_disp2, &
            & pair_disp3)
      end if
      if (config%verbosity > 0) then
         if (config%verbosity > 2) then
            call ascii_energy_atom(output_unit, mol, energies)
         end if
         call ascii_results(output_unit, mol, energy, gradient, sigma)
         if (config%pair_resolved) then
            call ascii_pairwise(output_unit, mol, pair_disp2, pair_disp3)
         end if
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
         if (allocated(config%grad_output)) then
            open(file=config%grad_output, newunit=unit)
            call tagged_result(unit, energy, gradient, sigma)
            close(unit)
            if (config%verbosity > 0) then
               write(output_unit, '(a)') &
                  & "[Info] Dispersion results written to '"//config%grad_output//"'"
            end if
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
         call json_results(unit, "  ", energy=energy, gradient=gradient, sigma=sigma, &
            & pairwise_energy2=pair_disp2, pairwise_energy3=pair_disp3, param=param)
         close(unit)
         if (config%verbosity > 0) then
            write(output_unit, '(a)') &
               & "[Info] JSON dump of results written to '"//config%json_output//"'"
         end if
      end if

   end if

end subroutine run_driver

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

subroutine param_driver(config, error)
   type(param_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(param_database) :: db
   class(damping_param), allocatable :: param

   call db%load(config%input, error)
   if (allocated(error)) return

   if (allocated(config%method)) then
      call db%get(param, config%method, config%damping)
      if (.not.allocated(param)) then
         call fatal_error(error, "No entry for '"//config%method//"' found in '"//config%input//"'")
         return
      end if
      call ascii_damping_param(output_unit, param, config%method)
   else
      write(output_unit, '(a, *(1x, g0))') "[Info] Found", size(db%records), &
         "damping parameters in '"//config%input//"'"
   end if

end subroutine param_driver

subroutine from_db(param, input, method, damping, error)
   class(damping_param), allocatable, intent(out) :: param
   character(len=*), intent(in) :: input
   character(len=*), intent(in) :: method
   character(len=*), intent(in) :: damping
   type(error_type), allocatable, intent(out) :: error

   type(param_database) :: db

   call db%load(input, error)
   if (allocated(error)) return

   call db%get(param, method, damping)
   if (.not.allocated(param)) then
      call fatal_error(error, "No entry for '"//method//"' found in '"//input//"'")
      return
   end if
end subroutine from_db

end module dftd3_app_driver
