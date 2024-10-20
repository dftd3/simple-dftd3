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

module test_gcp
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd3_gcp
   use dftd3_output, only : ascii_gcp_param
   implicit none
   private

   public :: collect_gcp

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp)) * 10


contains


!> Collect all exported unit tests
subroutine collect_gcp(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("hf3c", test_hf3c), &
      & new_unittest("pbeh3c", test_pbeh3c), &
      & new_unittest("hse3c", test_hse3c), &
      & new_unittest("b973c", test_b973c), &
      & new_unittest("r2scan3c", test_r2scan3c), &
      & new_unittest("hf/sv", test_hf_sv), &
      & new_unittest("hybrid/sv", test_hyb_sv), &
      & new_unittest("gga/sv", test_gga_sv), &
      & new_unittest("hf/sv(p)", test_hf_sv_p), &
      & new_unittest("dft/sv(p)", test_dft_sv_p), &
      & new_unittest("dft/svx", test_dft_svx), &
      & new_unittest("hf/svp", test_hf_svp), &
      & new_unittest("tpss/svp", test_tpss_svp), &
      & new_unittest("pw6b95/svp", test_pw6b95_svp), &
      & new_unittest("hybrid/svp", test_hyb_svp), &
      & new_unittest("gga/svp", test_gga_svp), &
      & new_unittest("hf/minis", test_hf_minis), &
      & new_unittest("tpss/minis", test_tpss_minis), &
      & new_unittest("pw6b95/minis", test_pw6b95_minis), &
      & new_unittest("hybrid/minis", test_hyb_minis), &
      & new_unittest("gga/minis", test_gga_minis), &
      & new_unittest("hf/631gd", test_hf_631gd), &
      & new_unittest("dft/631gd", test_dft_631gd), &
      & new_unittest("hf/tz", test_hf_tz), &
      & new_unittest("hybrid/tz", test_hyb_tz), &
      & new_unittest("gga/tz", test_gga_tz), &
      & new_unittest("hf/deftzvp", test_hf_deftzvp), &
      & new_unittest("dft/deftzvp", test_dft_deftzvp), &
      & new_unittest("hf/ccdz", test_hf_ccdz), &
      & new_unittest("dft/ccdz", test_dft_ccdz), &
      & new_unittest("hf/accdz", test_hf_accdz), &
      & new_unittest("dft/accdz", test_dft_accdz), &
      & new_unittest("dft/pobtz", test_dft_pobtz), &
      & new_unittest("hf/minix", test_hf_minix), &
      & new_unittest("dft/minix", test_dft_minix), &
      & new_unittest("hf/2g", test_hf_2g), &
      & new_unittest("hf/dzp", test_hf_dzp), &
      & new_unittest("dft/dzp", test_dft_dzp), &
      & new_unittest("hf/dz", test_hf_dz), &
      & new_unittest("dft/dz", test_dft_dz), &
      & new_unittest("b3pbe3c/def2mtzvp", test_b3pbe3c_def2mtzvp), &
      & new_unittest("dft/lanl2", test_dft_lanl2), &
      & new_unittest("without-args", test_without_args) &
      & ]

end subroutine collect_gcp


subroutine test_without_args(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(gcp_param) :: param

   call get_structure(mol, "MB16-43", "01")
   call get_gcp_param(param, mol)

   call check(error, .not.allocated(param%xv))
   if (allocated(error)) return
   call check(error, .not.allocated(param%emiss))
   if (allocated(error)) return
   call check(error, .not.allocated(param%slater))
   if (allocated(error)) return
   call check(error, .not.allocated(param%zeff))
   if (allocated(error)) return
end subroutine test_without_args


subroutine test_hf3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(gcp_param) :: param
   real(wp) :: energy

   call get_structure(mol, "MB16-43", "02")
   call get_gcp_param(param, mol, method="hf3c")

   call check(error, allocated(param%xv))
   if (allocated(error)) return
   call check(error, allocated(param%emiss))
   if (allocated(error)) return
   call check(error, allocated(param%slater))
   if (allocated(error)) return
   call check(error, allocated(param%zeff))
   if (allocated(error)) return
   call check(error, param%base)
   if (allocated(error)) return

   param%base = .false.
   energy = 0.0_wp
   call get_geometric_counterpoise(mol, param, energy)

   call check(error, energy, 0.0848896136_wp, thr=thr2)
   if (allocated(error)) return

   param%base = .true.
   energy = 0.0_wp
   call get_geometric_counterpoise(mol, param, energy)

   call check(error, energy, 0.0292452232_wp, thr=thr2)
   if (allocated(error)) return
end subroutine test_hf3c


subroutine test_pbeh3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(gcp_param) :: param
   real(wp) :: energy

   call get_structure(mol, "MB16-43", "03")
   call get_gcp_param(param, mol, method="pbeh3c")

   call check(error, allocated(param%xv))
   if (allocated(error)) return
   call check(error, allocated(param%emiss))
   if (allocated(error)) return
   call check(error, allocated(param%slater))
   if (allocated(error)) return
   call check(error, allocated(param%zeff))
   if (allocated(error)) return
   call check(error, .not.param%base)
   if (allocated(error)) return
   call check(error, .not.param%srb)
   if (allocated(error)) return

   ! call ascii_gcp_param(6, mol, param)

   energy = 0.0_wp
   call get_geometric_counterpoise(mol, param, energy)

   call check(error, energy, 0.0206352378_wp, thr=thr2)
end subroutine test_pbeh3c


subroutine test_hse3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(gcp_param) :: param
   real(wp) :: energy

   call get_structure(mol, "MB16-43", "04")
   call get_gcp_param(param, mol, method="hse3c")

   call check(error, allocated(param%xv), message="xv")
   if (allocated(error)) return
   call check(error, allocated(param%emiss), message="emiss")
   if (allocated(error)) return
   call check(error, allocated(param%slater), message="slater")
   if (allocated(error)) return
   call check(error, allocated(param%zeff), message="zeff")
   if (allocated(error)) return
   call check(error, .not.param%base)
   if (allocated(error)) return
   call check(error, .not.param%srb)
   if (allocated(error)) return

   energy = 0.0_wp
   call get_geometric_counterpoise(mol, param, energy)

   call check(error, energy, 0.0219213037_wp, thr=thr2)
end subroutine test_hse3c


subroutine test_b973c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(gcp_param) :: param
   real(wp) :: energy

   call get_structure(mol, "MB16-43", "05")
   call get_gcp_param(param, mol, method="b973c")

   call check(error, .not.allocated(param%xv), message="xv")
   if (allocated(error)) return
   call check(error, .not.allocated(param%emiss), message="emiss")
   if (allocated(error)) return
   call check(error, .not.allocated(param%slater), message="slater")
   if (allocated(error)) return
   call check(error, allocated(param%zeff), message="zeff")
   if (allocated(error)) return
   call check(error, .not.param%base)
   if (allocated(error)) return
   call check(error, param%srb)
   if (allocated(error)) return

   ! call ascii_gcp_param(6, mol, param)

   energy = 0.0_wp
   call get_geometric_counterpoise(mol, param, energy)

   call check(error, energy, -0.0369737816_wp, thr=thr2)
end subroutine test_b973c


subroutine test_r2scan3c(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(gcp_param) :: param
   real(wp) :: energy

   call get_structure(mol, "MB16-43", "06")
   call get_gcp_param(param, mol, method="r2scan3c")

   call check(error, allocated(param%xv))
   if (allocated(error)) return
   call check(error, allocated(param%emiss))
   if (allocated(error)) return
   call check(error, allocated(param%slater))
   if (allocated(error)) return
   call check(error, allocated(param%zeff))
   if (allocated(error)) return

   ! call ascii_gcp_param(6, mol, param)

   energy = 0.0_wp
   call get_geometric_counterpoise(mol, param, energy)

   call check(error, energy, 0.0113040952_wp, thr=thr2)
   if (allocated(error)) return
end subroutine test_r2scan3c


subroutine test_hf_sv(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "07")
   call test_generic_energy(error, mol, "hf", "sv", 0.0459359335_wp)
end subroutine test_hf_sv


subroutine test_hyb_sv(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "08")
   call test_generic_energy(error, mol, "b3lyp", "sv", 0.0817595304_wp)
end subroutine test_hyb_sv


subroutine test_gga_sv(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "09")
   call test_generic_energy(error, mol, "blyp", "sv", 0.0527050316_wp)
end subroutine test_gga_sv


subroutine test_hf_sv_p(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "10")
   call test_generic_energy(error, mol, "hf", "sv(p)", 0.0320180286_wp)
end subroutine test_hf_sv_p


subroutine test_dft_sv_p(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "11")
   call test_generic_energy(error, mol, "dft", "sv(p)", 0.0481858368_wp)
end subroutine test_dft_sv_p


subroutine test_dft_svx(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "12")
   call test_generic_energy(error, mol, "dft", "svx", 0.0442103186_wp)
end subroutine test_dft_svx


subroutine test_hf_svp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "13")
   call test_generic_energy(error, mol, "hf", "svp", 0.0389762270_wp)
end subroutine test_hf_svp


subroutine test_tpss_svp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "14")
   call test_generic_energy(error, mol, "tpss", "svp", 0.0361827169_wp)
end subroutine test_tpss_svp


subroutine test_pw6b95_svp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "15")
   call test_generic_energy(error, mol, "pw6b95", "svp", 0.0308635906_wp)
end subroutine test_pw6b95_svp


subroutine test_hyb_svp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "16")
   call test_generic_energy(error, mol, "b3lyp", "svp", 0.0310927236_wp)
end subroutine test_hyb_svp


subroutine test_gga_svp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "17")
   call test_generic_energy(error, mol, "blyp", "svp", 0.0956240694_wp)
end subroutine test_gga_svp


subroutine test_hf_minis(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "18")
   call test_generic_energy(error, mol, "hf", "minis", 0.0997731613_wp)
end subroutine test_hf_minis


subroutine test_tpss_minis(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "19")
   call test_generic_energy(error, mol, "tpss", "minis", 0.0888263651_wp)
end subroutine test_tpss_minis


subroutine test_pw6b95_minis(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "20")
   call test_generic_energy(error, mol, "pw6b95", "minis", 0.1107806553_wp)
end subroutine test_pw6b95_minis


subroutine test_hyb_minis(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "21")
   call test_generic_energy(error, mol, "b3lyp", "minis", 0.1019823129_wp)
end subroutine test_hyb_minis


subroutine test_gga_minis(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "22")
   call test_generic_energy(error, mol, "blyp", "minis", 0.1185302176_wp)
end subroutine test_gga_minis


subroutine test_hf_631gd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "23")
   call test_generic_energy(error, mol, "hf", "631gd", 0.0164851915_wp)
end subroutine test_hf_631gd


subroutine test_dft_631gd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "24")
   call test_generic_energy(error, mol, "dft", "631gd", 0.0691171349_wp)
end subroutine test_dft_631gd


subroutine test_hf_tz(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "25")
   call test_generic_energy(error, mol, "hf", "tz", 0.0079158812_wp)
end subroutine test_hf_tz


subroutine test_hyb_tz(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "26")
   call test_generic_energy(error, mol, "b3lyp", "tz", 0.0094799906_wp)
end subroutine test_hyb_tz


subroutine test_gga_tz(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "27")
   call test_generic_energy(error, mol, "blyp", "tz", 0.0010732743_wp)
end subroutine test_gga_tz


subroutine test_hf_deftzvp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "28")
   call test_generic_energy(error, mol, "hf", "deftzvp", 0.0173709952_wp)
end subroutine test_hf_deftzvp


subroutine test_dft_deftzvp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "29")
   call test_generic_energy(error, mol, "b3lyp", "deftzvp", 0.0162351039_wp)
end subroutine test_dft_deftzvp


subroutine test_hf_ccdz(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "30")
   call test_generic_energy(error, mol, "hf", "ccdz", 0.0531963593_wp)
end subroutine test_hf_ccdz


subroutine test_dft_ccdz(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "31")
   call test_generic_energy(error, mol, "b3lyp", "ccdz", 0.0369893530_wp)
end subroutine test_dft_ccdz


subroutine test_hf_accdz(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "32")
   call test_generic_energy(error, mol, "hf", "accdz", 0.0058924811_wp)
end subroutine test_hf_accdz


subroutine test_dft_accdz(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "33")
   call test_generic_energy(error, mol, "b3lyp", "accdz", 0.0065641768_wp)
end subroutine test_dft_accdz


subroutine test_dft_pobtz(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "34")
   call test_generic_energy(error, mol, "dft", "pobtz", 0.0630483570_wp)
end subroutine test_dft_pobtz


subroutine test_hf_minix(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "35")
   call test_generic_energy(error, mol, "hf", "minix", 0.0629778826_wp)
end subroutine test_hf_minix


subroutine test_dft_minix(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "36")
   call test_generic_energy(error, mol, "b3lyp", "minix", 0.1777559799_wp)
end subroutine test_dft_minix


subroutine test_hf_2g(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "37")
   call test_generic_energy(error, mol, "hf", "2g", 0.3699081462_wp, threshold=thr2*10)
end subroutine test_hf_2g


subroutine test_hf_dzp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "38")
   call test_generic_energy(error, mol, "hf", "dzp", 0.0237741643_wp)
end subroutine test_hf_dzp


subroutine test_dft_dzp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "39")
   call test_generic_energy(error, mol, "b3lyp", "dzp", 0.0621712827_wp)
end subroutine test_dft_dzp


subroutine test_hf_dz(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "40")
   call test_generic_energy(error, mol, "hf", "dz", 0.0246769955_wp)
end subroutine test_hf_dz


subroutine test_dft_dz(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "41")
   call test_generic_energy(error, mol, "b3lyp", "dz", 0.0721913136_wp)
end subroutine test_dft_dz


subroutine test_b3pbe3c_def2mtzvp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "42")
   call test_generic_energy(error, mol, "b3pbe3c", "def2mtzvp", 0.0712990259_wp)
end subroutine test_b3pbe3c_def2mtzvp


subroutine test_dft_lanl2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   call get_structure(mol, "MB16-43", "43")
   call test_generic_energy(error, mol, "b3lyp", "lanl", 0.0471216257_wp)
end subroutine test_dft_lanl2


subroutine test_generic_energy(error, mol, method, basis, reference_energy, threshold)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Method name
   character(len=*), intent(in) :: method

   !> Basis name
   character(len=*), intent(in) :: basis

   !> Reference energy
   real(wp), intent(in) :: reference_energy

   !> Threshold
   real(wp), intent(in), optional :: threshold

   type(gcp_param) :: param
   real(wp) :: energy, thr_

   thr_ = thr2
   if (present(threshold)) thr_ = threshold

   call get_gcp_param(param, mol, method=method, basis=basis)

   call check(error, allocated(param%xv), "Missing number of virtual orbitals")
   if (allocated(error)) return
   call check(error, allocated(param%emiss), "Missing BSSE parameters")
   if (allocated(error)) return
   call check(error, allocated(param%slater), "Missing Slater parameters")
   if (allocated(error)) return
   call check(error, allocated(param%zeff), "Missing effective nuclear charges")
   if (allocated(error)) return

   ! call ascii_gcp_param(6, mol, param)

   energy = 0.0_wp
   call get_geometric_counterpoise(mol, param, energy)

   call check(error, energy, reference_energy, thr=thr_)
   if (allocated(error)) then
      print *, energy, reference_energy, energy - reference_energy
   end if
end subroutine test_generic_energy


end module test_gcp