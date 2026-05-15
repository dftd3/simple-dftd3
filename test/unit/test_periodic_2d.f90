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

module test_periodic_2d
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use dftd3
   implicit none
   private

   public :: collect_periodic_2d

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr3 = 1000*sqrt(epsilon(1.0_wp))
   type(realspace_cutoff), parameter :: cutoff = &
      & realspace_cutoff(cn=30_wp, disp2=60.0_wp, disp3=15.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_periodic_2d(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("BP86-D3(0)", test_bp86d3zero_2d), &
      & new_unittest("gh185", test_gh185) &
      & ]

end subroutine collect_periodic_2d


subroutine test_dftd3_gen(error, mol, param, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Expected dispersion energy
   real(wp), intent(in) :: ref

   type(d3_model) :: d3
   real(wp) :: energy

   call new_d3_model(d3, mol)
   call get_dispersion(mol, d3, param, cutoff, energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
      print*,energy
   end if

end subroutine test_dftd3_gen


subroutine test_bp86d3zero_2d(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 4
   integer, parameter :: num(nat) = 6
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & -0.12918412100093_wp,  0.06210659750976_wp, -2.13384498734326_wp, &
      &  0.12856915667443_wp, -0.07403227791901_wp,  4.02358027265954_wp, &
      & -0.12317720857511_wp,  2.75170732207802_wp, -2.13345350602279_wp, &
      &  2.44816466162280_wp,  1.28612566399214_wp,  4.02317048854901_wp],&
      & shape(xyz))
   logical, parameter :: periodic(3) = [.true., .true., .false.]
   real(wp), parameter :: lattice(3, 3) = reshape([&
      &  4.68837849314507_wp, 0.00000000000000_wp, 0.00000000000000_wp, &
      & -2.36282788044783_wp, 4.04978545156612_wp, 0.00000000000000_wp, &
      &  0.00000000000000_wp, 0.00000000000000_wp, 1.00000000000000_wp],&
      & shape(lattice))
   real(wp), parameter :: lattice_vac(3, 3) = reshape([&
      &  4.68837849314507_wp, 0.00000000000000_wp, 0.00000000000000_wp, &
      & -2.36282788044783_wp, 4.04978545156612_wp, 0.00000000000000_wp, &
      &  0.00000000000000_wp, 0.00000000000000_wp, 100.000000000000_wp],&
      & shape(lattice))
   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & rs6 = 1.139_wp, s8 = 1.683_wp)

   call new(mol, num, xyz, periodic=periodic, lattice=lattice)
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -1.5983949078995030E-2_wp)
   if (allocated(error)) return

   call new(mol, num, xyz, periodic=[.true.], lattice=lattice_vac)
   call test_dftd3_gen(error, mol, param, -1.5983949078995030E-2_wp)
   if (allocated(error)) return
end subroutine test_bp86d3zero_2d


subroutine test_gh185(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 48
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "Na", "Na", "Na", "Na", "Na", "Na", "Na", "Na", "Na", "Na", "Na", "Na", "Nb", "Nb", &
      & "Nb", "Nb", "Nb", "Nb", "Nb", "Nb", "Nb", "Nb", "Nb", "Nb", "S", "S", "S", "S", &
      & "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", &
      & "S", "S", "S", "S", "S", "S"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      &  8.69858394342454_wp, 12.55856972495257_wp, 28.31851304478158_wp, &
      & 10.90531650094346_wp, 16.26501382266333_wp, 33.01544953180765_wp, &
      & 13.11204907735964_wp,-12.94575892316046_wp, 37.71238601883372_wp, &
      & 15.31878163487857_wp, -9.23931482544970_wp, 42.40932250585978_wp, &
      & 17.52551419239749_wp, -5.53287072773894_wp, 47.10625899288585_wp, &
      &  6.49185138590562_wp,  8.85212562724181_wp, 23.62157655775551_wp, &
      & 12.52689814187107_wp,-21.99147939517061_wp, 49.45472724584751_wp, &
      &  7.59521766466508_wp, -5.75326075511872_wp, 25.97004481071717_wp, &
      &  9.80195022218400_wp, -2.04681665740796_wp, 30.66698129774324_wp, &
      & 12.00868279860018_wp,  1.65962744030280_wp, 35.36391778476931_wp, &
      & 14.21541535611910_wp,  5.36607153801357_wp, 40.06085427179538_wp, &
      & 16.42214791363802_wp,  9.07251563572433_wp, 44.75779075882144_wp, &
      &  7.04353451583672_wp, 18.00804084838019_wp, 24.79581068423634_wp, &
      &  9.25026709225290_wp,-11.20273189744361_wp, 29.49274717126240_wp, &
      & 11.45699964977182_wp, -7.49628779973284_wp, 34.18968365828847_wp, &
      & 13.66373220729074_wp, -3.78984370202208_wp, 38.88662014531455_wp, &
      & 15.87046478370693_wp, -0.08339960431132_wp, 43.58355663234062_wp, &
      & 18.07719734122585_wp,  3.62304449339945_wp, 48.28049311936669_wp, &
      & 10.35363337101236_wp,  7.10909858262769_wp, 31.84121540532681_wp, &
      & 12.56036592853128_wp, 10.81554268033845_wp, 36.53815189235288_wp, &
      & 14.76709850494746_wp, 14.52198677804921_wp, 41.23508837937895_wp, &
      & 16.97383106246638_wp,-14.68878596777458_wp, 45.93202486640502_wp, &
      & 13.07858127180217_wp,-12.83556415513497_wp, 50.62896135343109_wp, &
      &  8.14690081349344_wp,  3.40265448491692_wp, 27.14427891830074_wp, &
      & 10.53425671992156_wp, 17.48677851070597_wp, 27.89639448903690_wp, &
      & 12.74098927744048_wp,-11.72399423511783_wp, 32.59333097606297_wp, &
      & 14.94772185385667_wp, -8.01755013740706_wp, 37.29026746308904_wp, &
      & 11.05247206319245_wp, -6.16432832476745_wp, 41.98720395011511_wp, &
      & 13.25920462071137_wp, -2.45788422705668_wp, 46.68414043714118_wp, &
      & 15.46593719712755_wp,  1.24855987065408_wp, 51.38107692416724_wp, &
      &  9.06964374334370_wp, 11.33680501801267_wp, 33.43756808755231_wp, &
      & 11.27637630086262_wp, 15.04324911572344_wp, 38.13450457457839_wp, &
      & 13.48310885838154_wp,-14.16752363010036_wp, 42.83144106160446_wp, &
      & 15.68984143479772_wp,-10.46107953238960_wp, 47.52837754863053_wp, &
      & 10.75816095759174_wp,  5.77713910766229_wp, 24.04369511350018_wp, &
      & 12.96489351511066_wp,  9.48358320537305_wp, 28.74063160052625_wp, &
      & 10.51194485584900_wp,  6.58783607487812_wp, 36.50744471210182_wp, &
      & 12.71867741336792_wp, 10.29428017258888_wp, 41.20438119912788_wp, &
      & 14.92540998978410_wp,-18.91649257323492_wp, 45.90131768615395_wp, &
      & 17.13214254730302_wp,-15.21004847552415_wp, 50.59825417318002_wp, &
      & 12.20046207009704_wp,  1.02817016452774_wp, 27.11357173804968_wp, &
      &  8.30521229833008_wp,  2.88139197716735_wp, 31.81050822507574_wp, &
      &  9.09195560741627_wp, 22.23574745384053_wp, 24.82651786448740_wp, &
      & 11.29868816493519_wp, -6.97502529198327_wp, 29.52345435151347_wp, &
      & 13.50542072245411_wp, -3.26858119427251_wp, 34.22039083853954_wp, &
      & 15.71215329887029_wp,  0.43786290343826_wp, 38.91732732556561_wp, &
      & 11.81690350820607_wp,  2.29108471607787_wp, 43.61426381259168_wp, &
      & 14.02363606572499_wp,  5.99752881378864_wp, 48.31120029961775_wp],&
      & shape(xyz))
   real(wp), parameter :: lattice(3, 3) = reshape([&
      &  6.10198234818314_wp, 0.00000000000000_wp, 18.99819597778184_wp, &
      &  1.85322228507115_wp,32.91721684353456_wp,-33.35210606988404_wp, &
      &  0.00000000000000_wp, 0.00000000000000_wp, 75.00265348192275_wp],&
      & shape(lattice))
   real(wp), parameter :: lattice_vac(3, 3) = reshape([&
      &  6.10198234818314_wp, 0.00000000000000_wp, 18.99819597778184_wp, &
      &  1.85322228507115_wp,32.91721684353456_wp,-33.35210606988404_wp, &
      &  0.00000000000000_wp, 0.00000000000000_wp,575.00265348192275_wp],&
      & shape(lattice_vac))
   logical, parameter :: periodic(3) = [.true., .true., .false.]

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s8 = 0.7875_wp, a1 = 0.4289_wp, a2 = 4.4407_wp, s9 = 0.0_wp)

   call new(mol, sym, xyz, periodic=[.true.], lattice=lattice_vac)
   print *, mol%lattice(:, 1)
   print *, mol%lattice(:, 2)
   print *, mol%lattice(:, 3)
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -0.34811410398591064_wp)
   ! if (allocated(error)) return

   call new(mol, sym, xyz, periodic=periodic, lattice=lattice)
   call test_dftd3_gen(error, mol, param, -0.34811410398591064_wp)
   if (allocated(error)) return
end subroutine test_gh185


end module test_periodic_2d