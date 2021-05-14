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

module test_param
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd3
   implicit none
   private

   public :: collect_param

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_param(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("DFT-D3(BJ)", test_d3bj_mb01), &
      & new_unittest("DFT-D3(0)", test_d3zero_mb09), &
      & new_unittest("DFT-D3(BJM)", test_d3bjm_mb02), &
      & new_unittest("DFT-D3(0M)", test_d3zerom_mb03), &
      & new_unittest("DFT-D3(BJ)-ATM", test_d3bjatm_mb17), &
      & new_unittest("DFT-D3(0)-ATM", test_d3zeroatm_mb25), &
      & new_unittest("DFT-D3(BJM)-ATM", test_d3bjmatm_mb04), &
      & new_unittest("DFT-D3(0M)-ATM", test_d3zeromatm_mb05), &
      & new_unittest("unknown-D3(BJ)", test_d3bj_unknown, should_fail=.true.), &
      & new_unittest("unknown-D3(0)", test_d3zero_unknown, should_fail=.true.), &
      & new_unittest("unknown-D3(BJM)", test_d3bjm_unknown, should_fail=.true.), &
      & new_unittest("unknown-D3(0M)", test_d3zerom_unknown, should_fail=.true.) &
      & ]

end subroutine collect_param


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
   call get_dispersion(mol, d3, param, realspace_cutoff(), energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
      print*,energy
   end if

end subroutine test_dftd3_gen


subroutine test_d3bj_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp

   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "b1b95", "b2gpplyp", "b2plyp", "b3lyp", "b3lyp/631gd", "b3pw91", "b97d", &
      & "bhlyp", "blyp", "bmk", "bop", "bp", "bpbe", "camb3lyp", "dftb3", &
      & "dsdblyp", "dsdblypfc", "hcth120", "hf", "hf3c", "hf3cv", "hf/minis", &
      & "hf/mixed", "hf/sv", "hse06", "hsesol", "lcwpbe", "mpw1b95", "mpwlyp", &
      & "olyp", "opbe", "otpss", "pbe", "pbe0", "pbeh3c", "pbesol", "ptpss", &
      & "pw1pw", "pw6b95", "pwb6k", "pwgga", "pwpb95", "revpbe", "revpbe0", &
      & "revpbe38", "revssb", "rpbe", "rpw86pbe", "ssb", "tpss", "tpss0", "tpssh", &
      & "scan", "rscan", "r2scan"]
   real(wp), parameter :: ref(*) = [&
      & -2.9551695097427895E-2_wp,-1.6638703327271943E-2_wp,-1.6725877952934207E-2_wp,&
      & -3.3014430058391525E-2_wp,-2.2051435529913663E-2_wp,-3.3481566296209062E-2_wp,&
      & -4.4319254716021633E-2_wp,-2.7768780587384141E-2_wp,-3.9702051630321589E-2_wp,&
      & -2.9455263912394638E-2_wp,-6.4253096882441671E-2_wp,-3.2374556196686173E-2_wp,&
      & -4.0184787904852644E-2_wp,-1.6566618994919072E-2_wp,-1.6283175780266980E-2_wp,&
      & -1.8049467081971568E-2_wp,-1.8968229532220542E-2_wp,-3.1793648577108383E-2_wp,&
      & -1.2496497910384512E-1_wp,-7.3593576345313688E-2_wp,-4.0603489309375830E-2_wp,&
      & -1.4163816710611299E-1_wp,-2.0131394878556811E-2_wp,-3.4278101637812045E-2_wp,&
      & -1.4158701187982837E-2_wp,-7.6707359228393463E-3_wp,-1.8974040257745569E-2_wp,&
      & -1.3757653102220277E-2_wp,-1.9765756641049737E-2_wp,-8.2952459757362679E-2_wp,&
      & -7.7457839141723550E-2_wp,-3.1481182175535194E-2_wp,-1.7882220406061061E-2_wp,&
      & -1.6523309281833876E-2_wp,-8.8508562330016100E-3_wp,-8.5090063235056107E-3_wp,&
      & -1.7255800051822631E-2_wp,-1.2406518256367386E-2_wp,-1.1629316087970519E-2_wp,&
      & -5.5085429466682810E-3_wp,-1.5678627159702693E-2_wp,-1.0788648599009749E-2_wp,&
      & -4.1947445140749925E-2_wp,-3.7430729018863761E-2_wp,-3.5154038334136578E-2_wp,&
      & -1.6039213182661993E-2_wp,-1.0199409897552987E-1_wp,-1.8473410189294232E-2_wp,&
      & -4.5113899837342458E-2_wp,-2.3465720037227366E-2_wp,-2.5011519452203857E-2_wp,&
      & -2.2089903668121399E-2_wp,-4.1432014419532683E-3_wp,-6.6123762443193899E-3_wp,&
      & -5.4507720332652972E-3_wp]

   call get_structure(mol, "MB16-43", "01")
   do ii = 1, size(func)
      call get_rational_damping(inp, trim(func(ii)), error)
      if (allocated(error)) exit
      call new_rational_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3bj_mb01


subroutine test_d3zero_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "slaterdiracexchange", "blyp", "bp", "b97d", "revpbe", "pbe", "pbesol", &
      & "rpw86pbe", "rpbe", "tpss", "b3lyp", "pbe0", "hse06", "revpbe38", &
      & "pw6b95", "tpss0", "b2plyp", "pwpb95", "b2gpplyp", "ptpss", "hf", &
      & "mpwlyp", "bpbe", "bhlyp", "tpssh", "pwb6k", "b1b95", "bop", "olyp", &
      & "opbe", "ssb", "revssb", "otpss", "b3pw91", "revpbe0", "pbe38", &
      & "mpw1b95", "mpwb1k", "bmk", "camb3lyp", "lcwpbe", "m05", "m052x", &
      & "m06l", "m06", "m062x", "m06hf", "hcth120", "scan"]
   real(wp), parameter :: ref(*) = [&
      & 1.4617020226660418E-1_wp,-1.4741291128345406E-2_wp,-1.3716392502956634E-2_wp,&
      &-2.0673076853447107E-2_wp,-1.8741322042994321E-2_wp,-6.7002109171808389E-3_wp,&
      &-4.9291044834241385E-3_wp,-7.6895623042570272E-3_wp,-2.0178785221185692E-2_wp,&
      &-9.7472918725363556E-3_wp,-1.2109385377469271E-2_wp,-7.2249765503104692E-3_wp,&
      &-4.4683705746465154E-3_wp,-1.2242157115222481E-2_wp,-5.7528731198206173E-3_wp,&
      &-9.4290371886078892E-3_wp,-6.6998082917956317E-3_wp,-4.6620816985011009E-3_wp,&
      &-4.8252542559612107E-3_wp,-5.6960455996363157E-3_wp,-1.3736390302808913E-2_wp,&
      &-8.7008860939233046E-3_wp,-1.7036708878390743E-2_wp,-9.7808060815089623E-3_wp,&
      &-9.6115287023438994E-3_wp,-3.6502663387597968E-3_wp,-1.1637689033466837E-2_wp,&
      &-2.4093411956739790E-2_wp,-3.5647772297611624E-2_wp,-3.3269300394152192E-2_wp,&
      &-6.3710244120544031E-3_wp,-5.6777824727715573E-3_wp,-1.2808061255275812E-2_wp,&
      &-1.3606462449138001E-2_wp,-1.5611940977561824E-2_wp,-7.3187679359965549E-3_wp,&
      &-7.1465707444989066E-3_wp,-6.7044219556221807E-3_wp,-1.3136115461328510E-2_wp,&
      &-8.3906504527315601E-3_wp,-8.8794919064170072E-3_wp,-4.6788894211193350E-3_wp,&
      &-9.1485044846090549E-4_wp,-4.7465342510505533E-4_wp,-1.3751400571877287E-3_wp,&
      &-4.0846462251094566E-4_wp,-8.1141664867494625E-4_wp,-9.5577230094779884E-3_wp,&
      &-1.3816199665340021E-3_wp]

   call get_structure(mol, "MB16-43", "09")
   do ii = 1, size(func)
      call get_zero_damping(inp, trim(func(ii)), error)
      if (allocated(error)) return
      call new_zero_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3zero_mb09


subroutine test_d3bjatm_mb17(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "b1b95", "b2gpplyp", "b2plyp", "b3lyp", "b3lyp/631gd", "b3pw91", "b97d", &
      & "bhlyp", "blyp", "bmk", "bop", "bp", "bpbe", "camb3lyp", "dftb3", &
      & "dsdblyp", "dsdblypfc", "hcth120", "hf", "hf3c", "hf3cv", "hf/minis", &
      & "hf/mixed", "hf/sv", "hse06", "hsesol", "lcwpbe", "mpw1b95", "mpwlyp", &
      & "olyp", "opbe", "otpss", "pbe", "pbe0", "pbeh3c", "pbesol", "ptpss", &
      & "pw1pw", "pw6b95", "pwb6k", "pwgga", "pwpb95", "revpbe", "revpbe0", &
      & "revpbe38", "revssb", "rpbe", "rpw86pbe", "ssb", "tpss", "tpss0", "tpssh", &
      & "scan", "rscan", "r2scan"]
   real(wp), parameter :: ref(*) = [&
      & -2.3886025810397844E-2_wp,-1.2511386988136818E-2_wp,-1.4044660769381583E-2_wp,&
      & -2.8422911196001738E-2_wp,-1.9411342505549978E-2_wp,-2.9014209453188066E-2_wp,&
      & -3.9594832657528958E-2_wp,-2.3209275734790969E-2_wp,-3.4356535835526364E-2_wp,&
      & -2.3548744390066521E-2_wp,-5.6222759698911239E-2_wp,-2.7737519130797358E-2_wp,&
      & -3.4944173006073054E-2_wp,-1.4207748710356847E-2_wp,-1.4807550615050694E-2_wp,&
      & -1.3592756421169463E-2_wp,-1.4297170752621349E-2_wp,-2.7202042216270681E-2_wp,&
      & -1.0723401451570813E-1_wp,-6.4226912061064456E-2_wp,-3.4476705337669271E-2_wp,&
      & -1.1444511684955419E-1_wp,-1.7961032237982111E-2_wp,-2.9706521389021402E-2_wp,&
      & -1.2161446851946924E-2_wp,-6.6471543537650248E-3_wp,-1.6387020697492129E-2_wp,&
      & -1.1272763435361845E-2_wp,-1.7452571852892073E-2_wp,-7.3575772751263660E-2_wp,&
      & -6.8870151906901350E-2_wp,-2.7534796959745162E-2_wp,-1.5749518311380703E-2_wp,&
      & -1.4436297073644914E-2_wp,-8.0176630195073831E-3_wp,-7.3547357715266467E-3_wp,&
      & -1.3145372076798436E-2_wp,-1.0640116316673890E-2_wp,-9.6504451818992445E-3_wp,&
      & -4.4996101691524041E-3_wp,-1.2740993313964484E-2_wp,-8.3393305736308033E-3_wp,&
      & -3.7212355838534522E-2_wp,-3.2852411771293132E-2_wp,-3.0614547789893518E-2_wp,&
      & -1.4346950270945317E-2_wp,-8.2244404039440305E-2_wp,-1.6294716345316058E-2_wp,&
      & -3.9038612599902538E-2_wp,-2.0567214691452097E-2_wp,-2.1558425276809661E-2_wp,&
      & -1.9336394869980178E-2_wp,-3.7332351763594033E-3_wp,-5.8233053769766027E-3_wp,&
      & -4.8268740049936159E-3_wp]

   call get_structure(mol, "MB16-43", "17")
   do ii = 1, size(func)
      call get_rational_damping(inp, trim(func(ii)), error, s9=1.0_wp)
      if (allocated(error)) return
      call new_rational_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3bjatm_mb17


subroutine test_d3zeroatm_mb25(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "slaterdiracexchange", "blyp", "bp", "b97d", "revpbe", "pbe", "pbesol", &
      & "rpw86pbe", "rpbe", "tpss", "b3lyp", "pbe0", "hse06", "revpbe38", &
      & "pw6b95", "tpss0", "b2plyp", "pwpb95", "b2gpplyp", "ptpss", "hf", &
      & "mpwlyp", "bpbe", "bhlyp", "tpssh", "pwb6k", "b1b95", "bop", "olyp", &
      & "opbe", "ssb", "revssb", "otpss", "b3pw91", "revpbe0", "pbe38", &
      & "mpw1b95", "mpwb1k", "bmk", "camb3lyp", "lcwpbe", "m05", "m052x", &
      & "m06l", "m06", "m062x", "m06hf", "hcth120", "scan"]
   real(wp), parameter :: ref(*) = [&
      & 1.0613172079515042E-1_wp,-1.8876287896469810E-2_wp,-1.7576404430843656E-2_wp,&
      &-2.3748729523009868E-2_wp,-2.2303397068410617E-2_wp,-8.6007354850641034E-3_wp,&
      &-6.1997814485238887E-3_wp,-9.8451848444439990E-3_wp,-2.2252314071025864E-2_wp,&
      &-1.2523774184016354E-2_wp,-1.5397294209036311E-2_wp,-9.1754201756565033E-3_wp,&
      &-5.8634276790166225E-3_wp,-1.5486804197012613E-2_wp,-7.1334203824096306E-3_wp,&
      &-1.2012621977313601E-2_wp,-8.3954644341124127E-3_wp,-5.7681880699224936E-3_wp,&
      &-6.0023918572656284E-3_wp,-7.0862806490480613E-3_wp,-1.7586419074267258E-2_wp,&
      &-1.1106947230183427E-2_wp,-2.1778378918668472E-2_wp,-1.2322554742451662E-2_wp,&
      &-1.2281846471432885E-2_wp,-4.4791942270594031E-3_wp,-1.4591821425616309E-2_wp,&
      &-2.9238167957521571E-2_wp,-3.7369719345347982E-2_wp,-3.6796498695353966E-2_wp,&
      &-8.1862302329976724E-3_wp,-7.2999669584442638E-3_wp,-1.6431165239300282E-2_wp,&
      &-1.7402831723160782E-2_wp,-1.8922205388850082E-2_wp,-9.2405961703656662E-3_wp,&
      &-8.9000499264262189E-3_wp,-8.3515334090636072E-3_wp,-1.6554343644100331E-2_wp,&
      &-1.0553688835175507E-2_wp,-1.1195100434896285E-2_wp,-5.8550454517077369E-3_wp,&
      &-1.0495107671637961E-3_wp,-4.4388399279473369E-4_wp,-1.7181987666819525E-3_wp,&
      &-3.6397398055820356E-4_wp,-9.0020655669676070E-4_wp,-1.2216291416317950E-2_wp,&
      &-1.7275318040978318E-3_wp]

   call get_structure(mol, "MB16-43", "25")
   do ii = 1, size(func)
      call get_zero_damping(inp, trim(func(ii)), error, s9=1.0_wp)
      if (allocated(error)) return
      call new_zero_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3zeroatm_mb25


subroutine test_d3bjm_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp

   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      "b2plyp", "b3lyp", "b97d", "blyp", "bp", "pbe", "pbe0", "lcwpbe"]
   real(wp), parameter :: ref(*) = [&
      &-2.4496063247118574E-2_wp,-4.7306509081943296E-2_wp,-5.9503372903920583E-2_wp,&
      &-5.6673388608130842E-2_wp,-4.7342640119970941E-2_wp,-2.4952685029454729E-2_wp,&
      &-2.3521899651908002E-2_wp,-2.7484708967011155E-2_wp]


   call get_structure(mol, "MB16-43", "02")
   do ii = 1, size(func)
      call get_rational_damping(inp, trim(func(ii)), error)
      if (allocated(error)) exit
      call new_rational_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3bjm_mb02


subroutine test_d3zerom_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      "b2plyp", "b3lyp", "b97d", "blyp", "bp", "pbe", "pbe0", "lcwpbe"]
   real(wp), parameter :: ref(*) = [&
      &-1.8298670704452319E-2_wp,-3.3517560882274484E-2_wp,-8.3081172538397016E-2_wp,&
      &-4.1380900243251549E-2_wp,-2.2212570818644292E-2_wp,-7.4972967973246998E-2_wp,&
      &-5.7366054727364238E-2_wp,-1.6369520315554695E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   do ii = 1, size(func)
      call get_mzero_damping(inp, trim(func(ii)), error)
      if (allocated(error)) return
      call new_mzero_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3zerom_mb03


subroutine test_d3bjmatm_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp

   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      "b2plyp", "b3lyp", "b97d", "blyp", "bp", "pbe", "pbe0", "lcwpbe"]
   real(wp), parameter :: ref(*) = [&
      &-1.8210584699554565E-2_wp,-4.9517673578593678E-2_wp,-1.1799782484083678E-1_wp,&
      &-5.0482003148319445E-2_wp,-2.1296304828924539E-2_wp,-4.2996504333752357E-2_wp,&
      &-4.2021636938441666E-2_wp,-1.8994833098554556E-2_wp]

   call get_structure(mol, "MB16-43", "04")
   do ii = 1, size(func)
      call get_mrational_damping(inp, trim(func(ii)), error, s9=1.0_wp)
      if (allocated(error)) exit
      call new_rational_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3bjmatm_mb04


subroutine test_d3zeromatm_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "b2plyp", "b3lyp", "b97d", "blyp", "bp", "pbe", "pbe0", "lcwpbe"]
   real(wp), parameter :: ref(*) = [&
      &-2.0123656701354201E-2_wp,-3.7721635943941914E-2_wp,-7.8081677058345220E-2_wp, &
      &-4.6610013189846915E-2_wp,-2.6471073414636268E-2_wp,-1.2081296874585111E-1_wp, &
      &-8.5996088468268228E-2_wp,-1.9517709619188587E-2_wp]

   call get_structure(mol, "MB16-43", "05")
   do ii = 1, size(func)
      call get_mzero_damping(inp, trim(func(ii)), error, s9=1.0_wp)
      if (allocated(error)) return
      call new_mzero_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3zeromatm_mb05


subroutine test_d3bj_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d3_param) :: inp

   call get_rational_damping(inp, "unknown", error)

end subroutine test_d3bj_unknown


subroutine test_d3zero_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d3_param) :: inp

   call get_zero_damping(inp, "unknown", error)

end subroutine test_d3zero_unknown


subroutine test_d3bjm_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d3_param) :: inp

   call get_mrational_damping(inp, "unknown", error)

end subroutine test_d3bjm_unknown


subroutine test_d3zerom_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d3_param) :: inp

   call get_mzero_damping(inp, "unknown", error)

end subroutine test_d3zerom_unknown


end module test_param
