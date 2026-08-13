program test_ghost_atoms
   use, intrinsic :: iso_fortran_env, only : r8 => real64
   use dftd3, only: d3_model, d3_param, rational_damping_param, get_rational_damping, &
      & new_rational_damping, new_d3_model, get_dispersion, realspace_cutoff
   use mctc_env, only: error_type
   use mctc_io, only: structure_type, new
   implicit none

   type(structure_type) :: mol
   type(error_type), allocatable :: error
   type(d3_model) :: disp
   type(d3_param) :: inp
   type(rational_damping_param) :: param
   real(r8) :: energy

   call new(mol, [8, 1, 1, 8, 1, 1], reshape([ &
      & -2.5142928288858868_r8,  1.8879764891638844_r8,  0.0000000000000000_r8, &
      & -2.0098096791094657_r8,  3.6303285583904779_r8,  0.0000000000000000_r8, &
      & -0.9545622395768059_r8, 0.92963588958819932_r8,  0.0000000000000000_r8, &
      &  1.8995798491583649_r8, -1.4196911493711153_r8,  0.0000000000000000_r8, &
      &  1.7895424492078418_r8, -2.5141248888855077_r8, -1.4467274293585810_r8, &
      &  1.7895424492078418_r8, -2.5141248888855077_r8,  1.4467274293585810_r8],&
      & [3, 6]), &
      & charge=0.0_r8, uhf=0)

   call get_rational_damping(inp, "PBE0", error, s9=1.0_r8)
   if (allocated(error)) then
      print "(2a)", "Error: ", error%message
      stop 1
   end if

   call new_rational_damping(param, inp)
   call new_d3_model(disp, mol, ghost=[4, 5, 6])
   call get_dispersion(mol, disp, param, realspace_cutoff(), energy)

   print "(a, f13.10, a)", "Dispersion energy with ghost atoms is ", energy, " Hartree"
end program test_ghost_atoms
