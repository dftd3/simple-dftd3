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

module test_model
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use dftd3_cutoff, only : get_lattice_points
   use dftd3_data, only : get_covalent_rad, get_vdw_rad, get_r4r2_val
   use dftd3_ncoord, only : get_coordination_number
   use dftd3_model
   implicit none
   private

   public :: collect_model

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_model(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("r4r2-val", test_r4r2_val), &
      & new_unittest("cov-rad", test_cov_rad), &
      & new_unittest("vdw-rad", test_vdw_rad), &
      & new_unittest("gw-mb01", test_gw_mb01), &
      & new_unittest("gw-mb02", test_gw_mb02), &
      & new_unittest("gw-mb03", test_gw_mb03), &
      & new_unittest("dgw-mb04", test_dgw_mb04), &
      & new_unittest("dgw-mb05", test_dgw_mb05) &
      & ]

end subroutine collect_model


subroutine test_gw_gen(error, mol, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Reference Gaussian weights
   real(wp), intent(in) :: ref(:, :)

   type(d3_model) :: d3
   real(wp), allocatable :: cn(:), rcov(:), gwvec(:, :)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)

   call new_d3_model(d3, mol)

   allocate(rcov(mol%nid), cn(mol%nat), gwvec(maxval(d3%ref), mol%nat))
   rcov(:) = get_covalent_rad(mol%num)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call get_coordination_number(mol, lattr, cutoff, rcov, cn)

   call d3%weight_references(mol, cn, gwvec)

   if (any(abs(gwvec - ref) > thr)) then
      call test_failed(error, "Gaussian weights do not match")
      print'(3es21.14)', gwvec
   end if

end subroutine test_gw_gen


subroutine test_dgw_gen(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   integer :: iat, mref
   type(d3_model) :: d3
   real(wp), allocatable :: cn(:), rcov(:), gwvec(:, :), gwdcn(:, :)
   real(wp), allocatable :: gwr(:, :), gwl(:, :), numdcn(:, :)
   real(wp), parameter :: cutoff = 30.0_wp, lattr(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp

   call new_d3_model(d3, mol)

   mref = maxval(d3%ref)
   allocate(rcov(mol%nid), cn(mol%nat), gwvec(mref, mol%nat), &
      & gwdcn(mref, mol%nat), gwr(mref, mol%nat), gwl(mref, mol%nat), &
      & numdcn(mref, mol%nat))
   rcov(:) = get_covalent_rad(mol%num)

   call get_coordination_number(mol, lattr, cutoff, rcov, cn)

   do iat = 1, mol%nat
      cn(iat) = cn(iat) + step
      call d3%weight_references(mol, cn, gwr)
      cn(iat) = cn(iat) - 2*step
      call d3%weight_references(mol, cn, gwl)
      cn(iat) = cn(iat) + step
      gwdcn(:, :) = 0.5_wp*(gwr - gwl)/step
      numdcn(:, iat) = gwdcn(:, iat)
      gwdcn(:, iat) = 0.0_wp
      if (any(abs(gwdcn) > thr)) then
         call test_failed(error, "Unexpected non-zero gradient element found")
         exit
      end if
   end do
   if (allocated(error)) return

   call d3%weight_references(mol, cn, gwvec, gwdcn)

   if (any(abs(gwdcn - numdcn) > thr2)) then
      call test_failed(error, "Gaussian weights derivatives do not match")
      print'(3es21.14)', gwdcn
      print'("---")'
      print'(3es21.14)', numdcn
      print'("---")'
      print'(3es21.14)', gwdcn - numdcn
   end if

end subroutine test_dgw_gen


subroutine test_gw_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(5, 16) = reshape([&
      & 4.6126800368181E-13_wp, 9.9999999999954E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.7843192515461E-01_wp, &
      & 2.1568074845392E-02_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.3325393706891E-08_wp, 1.5583083623077E-02_wp, &
      & 9.8441682305153E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 9.9942489973013E-01_wp, 5.7510026987275E-04_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.3577151993005E-02_wp, &
      & 9.8642284800700E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.8299214154930E-01_wp, 1.7007858450705E-02_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 9.9951946457202E-01_wp, 4.8053542797572E-04_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 1.1318179472266E-07_wp, &
      & 1.7150403606747E-02_wp, 9.8284948321146E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 1.2592840292219E-25_wp, 6.7327055773391E-14_wp, &
      & 1.9416632078902E-05_wp, 9.9998058336785E-01_wp, 0.0000000000000E+00_wp, &
      & 9.8640340476985E-01_wp, 1.3596595230150E-02_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.8337753046044E-01_wp, &
      & 1.6622469539564E-02_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 6.6363952683228E-06_wp, 9.9999336360473E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 4.7847441915845E-38_wp, 4.7250374771321E-24_wp, 2.6485900646248E-13_wp, &
      & 7.0840894527733E-06_wp, 9.9999291591028E-01_wp, 5.5794439953362E-26_wp, &
      & 1.4826440948761E-14_wp, 2.1971858851349E-06_wp, 1.5998035377261E-01_wp, &
      & 8.4001744904149E-01_wp, 1.1147440825485E-26_wp, 1.3347260063236E-14_wp, &
      & 8.8004839288989E-06_wp, 9.9999119951606E-01_wp, 0.0000000000000E+00_wp, &
      & 3.6442997835825E-41_wp, 1.6427702362541E-24_wp, 4.5062968981411E-11_wp, &
      & 9.9999999995494E-01_wp, 0.0000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "01")
   call test_gw_gen(error, mol, ref)

end subroutine test_gw_mb01


subroutine test_gw_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(5, 16) = reshape([&
      & 9.7676601040638E-01_wp, 2.3233989593621E-02_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 4.0048206631646E-21_wp, &
      & 3.2966477348527E-09_wp, 9.9999999670335E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 2.6642225826265E-25_wp, 4.9492317750973E-14_wp, &
      & 5.1218214268045E-06_wp, 2.6088549926377E-01_wp, 7.3910937891475E-01_wp, &
      & 1.5849430860130E-14_wp, 6.3791859302862E-06_wp, 9.9999362081405E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 3.6120178039298E-31_wp, &
      & 1.7457287541473E-14_wp, 9.9999999999998E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.8832485027058E-01_wp, 1.1675149729416E-02_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 9.7824019785714E-01_wp, 2.1759802142856E-02_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.7886593695194E-01_wp, &
      & 2.1134063048057E-02_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 7.8904759084419E-43_wp, 1.1580653027860E-27_wp, &
      & 3.5488902604170E-15_wp, 4.4104813362192E-06_wp, 9.9999558951866E-01_wp, &
      & 9.8980967099600E-01_wp, 1.0190329004003E-02_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 5.7626635072029E-28_wp, &
      & 4.1530345466038E-16_wp, 1.6741190805300E-07_wp, 3.2993948467453E-02_wp, &
      & 9.6700588412064E-01_wp, 7.5102223103533E-13_wp, 9.9999999999925E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 6.1105055092469E-06_wp, 9.9999388949449E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 9.7757382703570E-01_wp, &
      & 2.2426172964298E-02_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.9985614148140E-01_wp, 1.4385851860412E-04_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 6.5950662513562E-08_wp, 1.3359654289818E-02_wp, 9.8664027975952E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "02")
   call test_gw_gen(error, mol, ref)

end subroutine test_gw_mb02


subroutine test_gw_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(5, 16) = reshape([&
      & 3.2300021077159E-29_wp, 5.0433663928213E-17_wp, 4.9152842358743E-08_wp, &
      & 1.2091311559755E-02_wp, 9.8790863928740E-01_wp, 1.7643160161180E-14_wp, &
      & 6.7311590385752E-06_wp, 9.9999326884094E-01_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.8529818293400E-01_wp, 1.4701817066005E-02_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 1.0147310964179E-15_wp, 1.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 8.4311443142203E-45_wp, &
      & 2.1985658830340E-21_wp, 1.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 3.1341073999364E-34_wp, 7.1547887726564E-20_wp, &
      & 9.8580532921937E-09_wp, 9.9999999014195E-01_wp, 0.0000000000000E+00_wp, &
      & 2.5700026696663E-41_wp, 3.9933262735174E-26_wp, 4.6099078390827E-14_wp, &
      & 1.2393943802001E-05_wp, 9.9998760605615E-01_wp, 9.9555001924632E-01_wp, &
      & 4.4499807536766E-03_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 9.8371097706795E-01_wp, 1.6289022932048E-02_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 9.8279063160659E-01_wp, 1.7209368393411E-02_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 6.6515987191840E-06_wp, &
      & 9.9999334840128E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 8.2817609746890E-15_wp, 4.7528261312921E-06_wp, &
      & 9.9999524717386E-01_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 6.1348469416842E-41_wp, 7.6845613363214E-26_wp, 7.1126133041414E-14_wp, &
      & 1.5370563166116E-05_wp, 9.9998462943676E-01_wp, 9.9199504792296E-01_wp, &
      & 8.0049520770422E-03_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 5.6944988933592E-11_wp, 9.9999999994305E-01_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, 0.0000000000000E+00_wp, &
      & 9.9836938556521E-01_wp, 1.6306144347935E-03_wp, 0.0000000000000E+00_wp, &
      & 0.0000000000000E+00_wp, 0.0000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "03")
   call test_gw_gen(error, mol, ref)

end subroutine test_gw_mb03


subroutine test_dgw_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_dgw_gen(error, mol)

end subroutine test_dgw_mb04


subroutine test_dgw_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_dgw_gen(error, mol)

end subroutine test_dgw_mb05


subroutine test_r4r2_val(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check(error, get_r4r2_val("H"), get_r4r2_val(1))
   if (allocated(error)) return
   call check(error, get_r4r2_val("Og"), get_r4r2_val(118))
   if (allocated(error)) return
   call check(error, get_r4r2_val("X"), get_r4r2_val(-1))
end subroutine test_r4r2_val


subroutine test_cov_rad(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check(error, get_covalent_rad("C"), get_covalent_rad(6))
   if (allocated(error)) return
   call check(error, get_covalent_rad("Og"), get_covalent_rad(118))
   if (allocated(error)) return
   call check(error, get_covalent_rad("X"), get_covalent_rad(-1))
end subroutine test_cov_rad


subroutine test_vdw_rad(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check(error, get_vdw_rad("He", "Rn"), get_vdw_rad(2, 86))
   if (allocated(error)) return
   call check(error, get_vdw_rad("S", "Fl"), get_vdw_rad(16, 114))
   if (allocated(error)) return
   call check(error, get_vdw_rad("Am", "U"), get_vdw_rad(95, 92))
   if (allocated(error)) return
   call check(error, get_vdw_rad("Og", "Cn"), get_vdw_rad(118, 112))
   if (allocated(error)) return
   call check(error, get_vdw_rad("X", "X"), get_vdw_rad(-1, -1))
end subroutine test_vdw_rad


end module test_model
