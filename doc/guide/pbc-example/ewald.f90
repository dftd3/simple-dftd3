program pbc_ewald
   use, intrinsic :: iso_fortran_env, only : r8 => real64
   use dftd3, only: d3_model, d3_param, rational_damping_param, get_rational_damping, &
      & new_rational_damping, new_d3_model, get_dispersion, realspace_cutoff, &
      & d3_lowrank_config
   use mctc_env, only: error_type
   use mctc_io, only: structure_type, new
   implicit none

   character(len=:), allocatable :: method
   type(structure_type) :: mol
   type(error_type), allocatable :: error
   integer, allocatable :: num(:)
   real(r8), allocatable :: xyz(:, :)
   real(r8) :: lattice(3, 3), energy

   type(d3_model) :: disp
   type(d3_param) :: inp
   type(rational_damping_param) :: param

   method = "r2SCAN"
   num = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 7, 7, 7]
   xyz = reshape([ &  ! coordinates in Bohr
     & 3.5141339317522400_r8, 2.5620901186651399_r8, 0.9887044918760900_r8,  &
     & 2.5620901186651399_r8, 0.9887044918760900_r8, 3.5141339317522400_r8,  &
     & 0.9887044918760900_r8, 3.5141339317522400_r8, 2.5620901186651399_r8,  &
     & 8.2027323927494997_r8, 2.1265083423321101_r8, 8.3884924301184203_r8,  &
     & 7.2504996070913199_r8, 3.6998939691211699_r8, 5.8628740176711904_r8,  &
     & 5.6773029528733403_r8, 1.1742755566739300_r8, 6.8151068033293702_r8,  &
     & 6.8151068033293702_r8, 5.6773029528733403_r8, 1.1742755566739300_r8,  &
     & 8.3884924301184203_r8, 8.2027323927494997_r8, 2.1265083423321101_r8,  &
     & 3.6998939691211699_r8, 5.8628740176711904_r8, 7.2504996070913199_r8,  &
     & 1.1742755566739300_r8, 6.8151068033293702_r8, 5.6773029528733403_r8,  &
     & 2.1265083423321101_r8, 8.3884924301184203_r8, 8.2027323927494997_r8,  &
     & 5.8628740176711904_r8, 7.2504996070913199_r8, 3.6998939691211699_r8,  &
     & 1.9543543300807500_r8, 1.9543543300807500_r8, 1.9543543300807500_r8,  &
     & 6.6427638185069302_r8, 2.7342441309165002_r8, 7.4228425919137502_r8,  &
     & 7.4228425919137502_r8, 6.6427638185069302_r8, 2.7342441309165002_r8,  &
     & 2.7342441309165002_r8, 7.4228425919137502_r8, 6.6427638185069302_r8], &
     & [3, size(num)])
   lattice = reshape([ &  ! lattice vectors in Bohr
     & 9.3771283249512098_r8, 0.0000000000000000_r8, 0.0000000000000000_r8,  &
     & 0.0000000000000000_r8, 9.3771283249512098_r8, 0.0000000000000000_r8,  &
     & 0.0000000000000000_r8, 0.0000000000000000_r8, 9.3771283249512098_r8], &
     & [3, 3])
   call new(mol, num, xyz, lattice=lattice, charge=0.0_r8, uhf=0)

   call get_rational_damping(inp, method, error, s9=1.0_r8)
   if (allocated(error)) then
      print "(2a)", "Error: ", error%message
      return
   end if
   call new_rational_damping(param, inp)

   ! the low-rank expansion of the reference C6 coefficients enables the
   ! reciprocal space summation, an empty config uses the defaults
   call new_d3_model(disp, mol, lowrank=d3_lowrank_config())

   ! the two-body cutoff is unused here, the remaining cutoffs still apply
   call get_dispersion(error, mol, disp, param, realspace_cutoff(), energy)
   if (allocated(error)) then
      print "(2a)", "Error: ", error%message
      return
   end if

   print "(3a, f13.10, a)", "Dispersion energy for ", method, "-D3(BJ) is ", energy, " Hartree"

end program pbc_ewald
