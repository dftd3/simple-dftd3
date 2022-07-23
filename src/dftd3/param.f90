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

module dftd3_param
   use mctc_env, only : wp, error_type, fatal_error
   implicit none

   public :: d3_param
   public :: get_rational_damping, get_zero_damping
   public :: get_mrational_damping, get_mzero_damping
   public :: get_optimizedpower_damping


   type :: d3_param
      real(wp) :: s6 = 1.0_wp
      real(wp) :: s8 = 1.0_wp
      real(wp) :: s9 = 0.0_wp
      real(wp) :: rs6 = 1.0_wp
      real(wp) :: rs8 = 1.0_wp
      real(wp) :: a1 = 0.4_wp
      real(wp) :: a2 = 5.0_wp
      real(wp) :: alp = 14.0_wp
      real(wp) :: bet = 0.0_wp
   end type d3_param


   enum, bind(C)
      enumerator :: p_invalid, &
         & p_bp_df, p_blyp_df, p_revpbe_df, p_rpbe_df, p_b97d_df, p_pbe_df, &
         & p_rpw86pbe_df, p_b3lyp_df, p_tpss_df, p_hf_df, p_tpss0_df, &
         & p_pbe0_df, p_hse06_df, p_revpbe38_df, p_pw6b95_df, p_b2plyp_df, &
         & p_dsdblyp_df, p_dsdblypfc_df, p_bop_df, p_mpwlyp_df, p_olyp_df, &
         & p_pbesol_df, p_bpbe_df, p_opbe_df, p_ssb_df, p_revssb_df, p_otpss_df, &
         & p_b3pw91_df, p_bhlyp_df, p_revpbe0_df, p_tpssh_df, p_mpw1b95_df, &
         & p_pwb6k_df, p_b1b95_df, p_bmk_df, p_camb3lyp_df, p_lcwpbe_df, &
         & p_b2gpplyp_df, p_ptpss_df, p_pwpb95_df, p_hf_mixed_df, p_hf_sv_df, &
         & p_hf_minis_df, p_b3lyp_631gd_df, p_hcth120_df, p_dftb3_df, p_pw1pw_df, &
         & p_pwgga_df, p_hsesol_df, p_hf3c_df, p_hf3cv_df, p_pbeh3c_df, &
         & p_slaterdiracexchange_df, p_m05_df, p_m052x_df, p_m06l_df, p_m06_df, &
         & p_m062x_df, p_m06hf_df, p_pbe38_df, p_mpwb1k_df, p_scan_df, &
         & p_rscan_df, p_r2scan_df, p_wb97x_df, p_b97m_df, p_wb97m_df, &
         & p_pbehpbe_df, p_xlyp_df, p_mpwpw_df, p_hcth407_df, p_revtpss_df, &
         & p_tauhcth_df, p_b3p_df, p_b1p_df, p_b1lyp_df, p_mpw1pw_df, &
         & p_mpw1kcis_df, p_pbeh1pbe_df, p_pbe1kcis_df, p_x3lyp_df, p_o3lyp_df, &
         & p_b97_1_df, p_b97_2_df, p_b98_df, p_hiss_df, p_hse03_df, p_revtpssh_df, &
         & p_revtpss0_df, p_tpss1kcis_df, p_tauhcthhyb_df, p_mn15_df, p_lc_whpbe_df, &
         & p_mpw2plyp_df, p_m11_df, p_sogga11x_df, p_n12sx_df, p_mn12sx_df, &
         & p_ms2_df, p_ms2h_df, p_mpw1lyp_df, p_mpwkcis1k_df, p_pkzb_df, p_n12_df, &
         & p_m08hx_df, p_m11l_df, p_mn15l_df, p_pwp_df, p_r2scanh_df, p_r2scan0_df, &
         & p_r2scan50_df, p_b973c_df
   end enum

contains


function get_method_id(method) result(id)

   !> Name of the method to look up
   character(len=*), intent(in) :: method

   character(len=len(method)) :: lc_method

   integer :: id
   integer :: i, j

   j = 0
   do i = 1, len(method)
      if (method(i:i) /= "-") then
         j = j + 1
         lc_method(j:j) = method(i:i)
      end if
   end do
   select case(trim(lowercase(lc_method)))
   case default; id = p_invalid
   case("b1b95"); id = p_b1b95_df
   case("b1lyp"); id = p_b1lyp_df
   case("b1p", "b1p86"); id = p_b1p_df
   case("b2gpplyp"); id = p_b2gpplyp_df
   case("b2plyp"); id = p_b2plyp_df
   case("b3lyp"); id = p_b3lyp_df
   case("b3lyp/631gd"); id = p_b3lyp_631gd_df
   case("b3p", "b3p86"); id = p_b3p_df
   case("b3pw91"); id = p_b3pw91_df
   case("b971"); id = p_b97_1_df
   case("b972"); id = p_b97_2_df
   case("b97d"); id = p_b97d_df
   case("b973c"); id = p_b973c_df
   case("b97m"); id = p_b97m_df
   case("b98"); id = p_b98_df
   case("bhlyp"); id = p_bhlyp_df
   case("blyp"); id = p_blyp_df
   case("bmk"); id = p_bmk_df
   case("bop"); id = p_bop_df
   case("bp", "bp86"); id = p_bp_df
   case("bpbe"); id = p_bpbe_df
   case("camb3lyp"); id = p_camb3lyp_df
   case("dftb3"); id = p_dftb3_df
   case("dsdblyp"); id = p_dsdblyp_df
   case("dsdblypfc"); id = p_dsdblypfc_df
   case("hcth120"); id = p_hcth120_df
   case("hcth407", "hcth/407"); id = p_hcth407_df
   case("hf"); id = p_hf_df
   case("hf/minis"); id = p_hf_minis_df
   case("hf/mixed"); id = p_hf_mixed_df
   case("hf/sv"); id = p_hf_sv_df
   case("hf3c"); id = p_hf3c_df
   case("hf3cv"); id = p_hf3cv_df
   case("hiss"); id = p_hiss_df
   case("hse03"); id = p_hse03_df
   case("hse06"); id = p_hse06_df
   case("hsesol"); id = p_hsesol_df
   case("lcwhpbe", "lcomegahpbe", "lcωhpbe"); id = p_lc_whpbe_df
   case("lcwpbe"); id = p_lcwpbe_df
   case("m05"); id = p_m05_df
   case("m052x"); id = p_m052x_df
   case("m06"); id = p_m06_df
   case("m062x"); id = p_m062x_df
   case("m06hf"); id = p_m06hf_df
   case("m06l"); id = p_m06l_df
   case("m08hx"); id = p_m08hx_df
   case("m11"); id = p_m11_df
   case("m11l"); id = p_m11l_df
   case("mn12sx"); id = p_mn12sx_df
   case("mn15"); id = p_mn15_df
   case("mn15l"); id = p_mn15l_df
   case("mpw1b95"); id = p_mpw1b95_df
   case("mpw1kcis"); id = p_mpw1kcis_df
   case("mpw1pw", "mpw1pw91"); id = p_mpw1pw_df
   case("mpw2plyp"); id = p_mpw2plyp_df
   case("mpwb1k"); id = p_mpwb1k_df
   case("mpwlyp"); id = p_mpwlyp_df
   case("mpwpw", "mpwpw91"); id = p_mpwpw_df
   case("mpw1lyp"); id = p_mpw1lyp_df
   case("mpwkcis1k"); id = p_mpwkcis1k_df
   case("ms2"); id = p_ms2_df
   case("ms2h"); id = p_ms2h_df
   case("n12"); id = p_n12_df
   case("n12sx"); id = p_n12sx_df
   case("o3lyp"); id = p_o3lyp_df
   case("olyp"); id = p_olyp_df
   case("opbe"); id = p_opbe_df
   case("otpss"); id = p_otpss_df
   case("pbe"); id = p_pbe_df
   case("pbe0"); id = p_pbe0_df
   case("pbe1kcis"); id = p_pbe1kcis_df
   case("pbe38"); id = p_pbe38_df
   case("pbeh1pbe"); id = p_pbeh1pbe_df
   case("pbeh3c"); id = p_pbeh3c_df
   case("pbehpbe"); id = p_pbehpbe_df
   case("pbesol"); id = p_pbesol_df
   case("pkzb"); id = p_pkzb_df
   case("ptpss"); id = p_ptpss_df
   case("pwp", "pw91p86"); id = p_pwp_df
   case("pw1pw"); id = p_pw1pw_df
   case("pw6b95"); id = p_pw6b95_df
   case("pwb6k"); id = p_pwb6k_df
   case("pwgga"); id = p_pwgga_df
   case("pwpb95"); id = p_pwpb95_df
   case("r2scan"); id = p_r2scan_df
   case("r2scanh"); id = p_r2scanh_df
   case("r2scan0"); id = p_r2scan0_df
   case("r2scan50"); id = p_r2scan50_df
   case("revpbe"); id = p_revpbe_df
   case("revpbe0"); id = p_revpbe0_df
   case("revpbe38"); id = p_revpbe38_df
   case("revssb"); id = p_revssb_df
   case("revtpss"); id = p_revtpss_df
   case("revtpss0"); id = p_revtpss0_df
   case("revtpssh"); id = p_revtpssh_df
   case("rpbe"); id = p_rpbe_df
   case("rpw86pbe"); id = p_rpw86pbe_df
   case("rscan"); id = p_rscan_df
   case("scan"); id = p_scan_df
   case("slaterdiracexchange"); id = p_slaterdiracexchange_df
   case("sogga11x"); id = p_sogga11x_df
   case("ssb"); id = p_ssb_df
   case("tauhcth", "τhcth"); id = p_tauhcth_df
   case("tauhcthhyb", "τhcthhyb"); id = p_tauhcthhyb_df
   case("tpss"); id = p_tpss_df
   case("tpss0"); id = p_tpss0_df
   case("tpss1kcis"); id = p_tpss1kcis_df
   case("tpssh"); id = p_tpssh_df
   case("wb97m"); id = p_wb97m_df
   case("wb97x"); id = p_wb97x_df
   case("x3lyp"); id = p_x3lyp_df
   case("xlyp"); id = p_xlyp_df
   end select

end function get_method_id


subroutine get_rational_damping(param, method, error, s9)

   !> Loaded parameter record
   type(d3_param), intent(out) :: param

   !> Name of the method to look up
   character(len=*), intent(in) :: method

   !> Overwrite s9
   real(wp), intent(in), optional :: s9

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(get_method_id(method))
   case default
      call fatal_error(error, "No entry for '"//method//"' present")
      return
   case(p_bp_df)
      param = d3_param(a1=0.3946_wp, s8=3.2822_wp, a2=4.8516_wp)
   case(p_blyp_df)
      param = d3_param(a1=0.4298_wp, s8=2.6996_wp, a2=4.2359_wp)
   case(p_revpbe_df)
      param = d3_param(a1=0.5238_wp, s8=2.3550_wp, a2=3.5016_wp)
   case(p_rpbe_df)
      param = d3_param(a1=0.1820_wp, s8=0.8318_wp, a2=4.0094_wp)
   case(p_b97d_df)
      param = d3_param(a1=0.5545_wp, s8=2.2609_wp, a2=3.2297_wp)
   case(p_b973c_df)
      param = d3_param(a1=0.37_wp, s8=1.50_wp, a2=4.10_wp)
   case(p_pbe_df)
      param = d3_param(a1=0.4289_wp, s8=0.7875_wp, a2=4.4407_wp)
   case(p_rpw86pbe_df)
      param = d3_param(a1=0.4613_wp, s8=1.3845_wp, a2=4.5062_wp)
   case(p_b3lyp_df)
      param = d3_param(a1=0.3981_wp, s8=1.9889_wp, a2=4.4211_wp)
   case(p_tpss_df)
      param = d3_param(a1=0.4535_wp, s8=1.9435_wp, a2=4.4752_wp)
   case(p_hf_df)
      param = d3_param(a1=0.3385_wp, s8=0.9171_wp, a2=2.8830_wp)
   case(p_tpss0_df)
      param = d3_param(a1=0.3768_wp, s8=1.2576_wp, a2=4.5865_wp)
   case(p_pbe0_df)
      param = d3_param(a1=0.4145_wp, s8=1.2177_wp, a2=4.8593_wp)
   case(p_hse06_df)
      param = d3_param(a1=0.383_wp, s8=2.310_wp, a2=5.685_wp)
   case(p_revpbe38_df)
      param = d3_param(a1=0.4309_wp, s8=1.4760_wp, a2=3.9446_wp)
   case(p_pw6b95_df)
      param = d3_param(a1=0.2076_wp, s8=0.7257_wp, a2=6.3750_wp)
   case(p_b2plyp_df)
      param = d3_param(a1=0.3065_wp, s8=0.9147_wp, a2=5.0570_wp, s6=0.64_wp)
   case(p_dsdblyp_df)
      param = d3_param(a1=0.0000_wp, s8=0.2130_wp, a2=6.0519_wp, s6=0.50_wp)
   case(p_dsdblypfc_df)
      param = d3_param(a1=0.0009_wp, s8=0.2112_wp, a2=5.9807_wp, s6=0.50_wp)
   case(p_bop_df)
      param = d3_param(a1=0.4870_wp, s8=3.2950_wp, a2=3.5043_wp)
   case(p_mpwlyp_df)
      param = d3_param(a1=0.4831_wp, s8=2.0077_wp, a2=4.5323_wp)
   case(p_olyp_df)
      param = d3_param(a1=0.5299_wp, s8=2.6205_wp, a2=2.8065_wp)
   case(p_pbesol_df)
      param = d3_param(a1=0.4466_wp, s8=2.9491_wp, a2=6.1742_wp)
   case(p_bpbe_df)
      param = d3_param(a1=0.4567_wp, s8=4.0728_wp, a2=4.3908_wp)
   case(p_opbe_df)
      param = d3_param(a1=0.5512_wp, s8=3.3816_wp, a2=2.9444_wp)
   case(p_ssb_df)
      param = d3_param(a1=-0.0952_wp, s8=-0.1744_wp, a2=5.2170_wp)
   case(p_revssb_df)
      param = d3_param(a1=0.4720_wp, s8=0.4389_wp, a2=4.0986_wp)
   case(p_otpss_df)
      param = d3_param(a1=0.4634_wp, s8=2.7495_wp, a2=4.3153_wp)
   case(p_b3pw91_df)
      param = d3_param(a1=0.4312_wp, s8=2.8524_wp, a2=4.4693_wp)
   case(p_bhlyp_df)
      param = d3_param(a1=0.2793_wp, s8=1.0354_wp, a2=4.9615_wp)
   case(p_revpbe0_df)
      param = d3_param(a1=0.4679_wp, s8=1.7588_wp, a2=3.7619_wp)
   case(p_tpssh_df)
      param = d3_param(a1=0.4529_wp, s8=2.2382_wp, a2=4.6550_wp)
   case(p_mpw1b95_df)
      param = d3_param(a1=0.1955_wp, s8=1.0508_wp, a2=6.4177_wp)
   case(p_pwb6k_df)
      param = d3_param(a1=0.1805_wp, s8=0.9383_wp, a2=7.7627_wp)
   case(p_b1b95_df)
      param = d3_param(a1=0.2092_wp, s8=1.4507_wp, a2=5.5545_wp)
   case(p_bmk_df)
      param = d3_param(a1=0.1940_wp, s8=2.0860_wp, a2=5.9197_wp)
   case(p_camb3lyp_df)
      param = d3_param(a1=0.3708_wp, s8=2.0674_wp, a2=5.4743_wp)
   case(p_lcwpbe_df)
      param = d3_param(a1=0.3919_wp, s8=1.8541_wp, a2=5.0897_wp)
   case(p_b2gpplyp_df)
      param = d3_param(a1=0.0000_wp, s8=0.2597_wp, a2=6.3332_wp, s6=0.560_wp)
   case(p_ptpss_df)
      param = d3_param(a1=0.0000_wp, s8=0.2804_wp, a2=6.5745_wp, s6=0.750_wp)
   case(p_pwpb95_df)
      param = d3_param(a1=0.0000_wp, s8=0.2904_wp, a2=7.3141_wp, s6=0.820_wp)
   case(p_hf_mixed_df)
      param = d3_param(a1=0.5607_wp, s8=3.9027_wp, a2=4.5622_wp)
   case(p_hf_sv_df)
      param = d3_param(a1=0.4249_wp, s8=2.1849_wp, a2=4.2783_wp)
   case(p_hf_minis_df)
      param = d3_param(a1=0.1702_wp, s8=0.9841_wp, a2=3.8506_wp)
   case(p_b3lyp_631gd_df)
      param = d3_param(a1=0.5014_wp, s8=4.0672_wp, a2=4.8409_wp)
   case(p_hcth120_df)
      param = d3_param(a1=0.3563_wp, s8=1.0821_wp, a2=4.3359_wp)
   case(p_dftb3_df)
      param = d3_param(a1=0.5719_wp, s8=0.5883_wp, a2=3.6017_wp)
   case(p_pw1pw_df)
      param = d3_param(a1=0.3807_wp, s8=2.3363_wp, a2=5.8844_wp)
   case(p_pwgga_df)
      param = d3_param(a1=0.2211_wp, s8=2.6910_wp, a2=6.7278_wp)
   case(p_hsesol_df)
      param = d3_param(a1=0.4650_wp, s8=2.9215_wp, a2=6.2003_wp)
   case(p_hf3c_df)
      param = d3_param(a1=0.4171_wp, s8=0.8777_wp, a2=2.9149_wp)
   case(p_hf3cv_df)
      param = d3_param(a1=0.3063_wp, s8=0.5022_wp, a2=3.9856_wp)
   case(p_pbeh3c_df)
      param = d3_param(a1=0.4860_wp, s8=0.0000_wp, a2=4.5000_wp)
   case(p_scan_df)
      param = d3_param(a1=0.5380_wp, s8=0.0000_wp, a2=5.4200_wp)
   case(p_rscan_df)
      param = d3_param(a1=0.47023427_wp, s8=1.08859014_wp, a2=5.73408312_wp)
   case(p_r2scan_df)
      param = d3_param(a1=0.49484001_wp, s8=0.78981345_wp, a2=5.73083694_wp)
   case(p_r2scanh_df)
      param = d3_param(s8=1.1236_wp, a1=0.4709_wp, a2=5.9157_wp)
   case(p_r2scan0_df)
      param = d3_param(s8=1.1846_wp, a1=0.4534_wp, a2=5.8972_wp)
   case(p_r2scan50_df)
      param = d3_param(s8=1.3294_wp, a1=0.4311_wp, a2=5.9240_wp)
   case(p_wb97x_df)
      param = d3_param(a1=0.0000_wp, s8=0.2641_wp, a2=5.4959_wp)
   case(p_wb97m_df)
      param = d3_param(a1=0.5660_wp, s8=0.3908_wp, a2=3.1280_wp)
   case(p_b97m_df)
      param = d3_param(a1=-0.0780_wp, s8=0.1384_wp, a2=5.5946_wp)
   case(p_pbehpbe_df)
      param = d3_param(a1=0.0000_wp, s8=1.1152_wp, a2=6.7184_wp)
   case(p_xlyp_df)
      param = d3_param(a1=0.0809_wp, s8=1.5669_wp, a2=5.3166_wp)
   case(p_mpwpw_df)
      param = d3_param(a1=0.3168_wp, s8=1.7974_wp, a2=4.7732_wp)
   case(p_hcth407_df)
      param = d3_param(a1=0.0000_wp, s8=0.6490_wp, a2=4.8162_wp)
   case(p_revtpss_df)
      param = d3_param(a1=0.4326_wp, s8=1.4023_wp, a2=4.4723_wp)
   case(p_tauhcth_df)
      param = d3_param(a1=0.0000_wp, s8=1.2626_wp, a2=5.6162_wp)
   case(p_b3p_df)
      param = d3_param(a1=0.4601_wp, s8=3.3211_wp, a2=4.9858_wp)
   case(p_b1p_df)
      param = d3_param(a1=0.4724_wp, s8=3.5681_wp, a2=4.9858_wp)
   case(p_b1lyp_df)
      param = d3_param(a1=0.1986_wp, s8=2.1167_wp, a2=5.3875_wp)
   case(p_mpwb1k_df)
      param = d3_param(a1=0.1474_wp, s8=0.9499_wp, a2=6.6223_wp)
   case(p_mpw1pw_df)
      param = d3_param(a1=0.3342_wp, s8=1.8744_wp, a2=4.9819_wp)
   case(p_mpw1kcis_df)
      param = d3_param(a1=0.0576_wp, s8=1.0893_wp, a2=5.5314_wp)
   case(p_pbeh1pbe_df)
      param = d3_param(a1=0.0000_wp, s8=1.4877_wp, a2=7.0385_wp)
   case(p_pbe1kcis_df)
      param = d3_param(a1=0.0000_wp, s8=0.7688_wp, a2=6.2794_wp)
   case(p_x3lyp_df)
      param = d3_param(a1=0.2022_wp, s8=1.5744_wp, a2=5.4184_wp)
   case(p_o3lyp_df)
      param = d3_param(a1=0.0963_wp, s8=1.8171_wp, a2=5.9940_wp)
   case(p_b97_1_df)
      param = d3_param(a1=0.0000_wp, s8=0.4814_wp, a2=6.2279_wp)
   case(p_b97_2_df)
      param = d3_param(a1=0.0000_wp, s8=0.9448_wp, a2=5.9940_wp)
   case(p_b98_df)
      param = d3_param(a1=0.0000_wp, s8=0.7086_wp, a2=6.0672_wp)
   case(p_hiss_df)
      param = d3_param(a1=0.0000_wp, s8=1.6112_wp, a2=7.3539_wp)
   case(p_hse03_df)
      param = d3_param(a1=0.0000_wp, s8=1.1243_wp, a2=6.8889_wp)
   case(p_revtpssh_df)
      param = d3_param(a1=0.2660_wp, s8=1.4076_wp, a2=5.3761_wp)
   case(p_revtpss0_df)
      param = d3_param(a1=0.2218_wp, s8=1.6151_wp, a2=5.7985_wp)
   case(p_tpss1kcis_df)
      param = d3_param(a1=0.0000_wp, s8=1.0542_wp, a2=6.0201_wp)
   case(p_tauhcthhyb_df)
      param = d3_param(a1=0.0000_wp, s8=0.9585_wp, a2=10.1389_wp)
   case(p_m11_df)
      param = d3_param(a1=0.0000_wp, s8=2.8112_wp, a2=10.1389_wp)
   case(p_sogga11x_df)
      param = d3_param(a1=0.1330_wp, s8=1.1426_wp, a2=5.7381_wp)
   case(p_n12sx_df)
      param = d3_param(a1=0.3283_wp, s8=2.4900_wp, a2=5.7898_wp)
   case(p_mn12sx_df)
      param = d3_param(a1=0.0983_wp, s8=1.1674_wp, a2=8.0259_wp)
   case(p_mn15_df)
      param = d3_param(a1=2.0971_wp, s8=0.7862_wp, a2=7.5923_wp)
   case(p_lc_whpbe_df)
      param = d3_param(a1=0.2746_wp, s8=1.1908_wp, a2=5.3157_wp)
   case(p_mpw2plyp_df)
      param = d3_param(s6=0.66_wp, a1=0.4105_wp, s8=0.6223_wp, a2=5.0136_wp)
   end select

   if (present(s9)) then
      param%s9 = s9
   end if

end subroutine get_rational_damping


subroutine get_zero_damping(param, method, error, s9)

   !> Loaded parameter record
   type(d3_param), intent(out) :: param

   !> Name of the method to look up
   character(len=*), intent(in) :: method

   !> Overwrite s9
   real(wp), intent(in), optional :: s9

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(get_method_id(method))
   case default
      call fatal_error(error, "No entry for '"//method//"' present")
      return
   case(p_slaterdiracexchange_df)
      param = d3_param(rs6=0.999_wp, s8=-1.957_wp, rs8=0.697_wp)
   case(p_blyp_df)
      param = d3_param(rs6=1.094_wp, s8=1.682_wp)
   case(p_bp_df)
      param = d3_param(rs6=1.139_wp, s8=1.683_wp)
   case(p_b97d_df)
      param = d3_param(rs6=0.892_wp, s8=0.909_wp)
   case(p_b973c_df)
      param = d3_param(rs6=1.06_wp, s8=1.50_wp)
   case(p_revpbe_df)
      param = d3_param(rs6=0.923_wp, s8=1.010_wp)
   case(p_pbe_df)
      param = d3_param(rs6=1.217_wp, s8=0.722_wp)
   case(p_pbesol_df)
      param = d3_param(rs6=1.345_wp, s8=0.612_wp)
   case(p_rpw86pbe_df)
      param = d3_param(rs6=1.224_wp, s8=0.901_wp)
   case(p_rpbe_df)
      param = d3_param(rs6=0.872_wp, s8=0.514_wp)
   case(p_tpss_df)
      param = d3_param(rs6=1.166_wp, s8=1.105_wp)
   case(p_b3lyp_df)
      param = d3_param(rs6=1.261_wp, s8=1.703_wp)
   case(p_pbe0_df)
      param = d3_param(rs6=1.287_wp, s8=0.928_wp)
   case(p_hse06_df)
      param = d3_param(rs6=1.129_wp, s8=0.109_wp)
   case(p_revpbe38_df)
      param = d3_param(rs6=1.021_wp, s8=0.862_wp)
   case(p_pw6b95_df)
      param = d3_param(rs6=1.532_wp, s8=0.862_wp)
   case(p_tpss0_df)
      param = d3_param(rs6=1.252_wp, s8=1.242_wp)
   case(p_b2plyp_df)
      param = d3_param(rs6=1.427_wp, s8=1.022_wp, s6=0.64_wp)
   case(p_pwpb95_df)
      param = d3_param(rs6=1.557_wp, s8=0.705_wp, s6=0.82_wp)
   case(p_b2gpplyp_df)
      param = d3_param(rs6=1.586_wp, s8=0.760_wp, s6=0.56_wp)
   case(p_ptpss_df)
      param = d3_param(rs6=1.541_wp, s8=0.879_wp, s6=0.75_wp)
   case(p_hf_df)
      param = d3_param(rs6=1.158_wp, s8=1.746_wp)
   case(p_mpwlyp_df)
      param = d3_param(rs6=1.239_wp, s8=1.098_wp)
   case(p_bpbe_df)
      param = d3_param(rs6=1.087_wp, s8=2.033_wp)
   case(p_bhlyp_df)
      param = d3_param(rs6=1.370_wp, s8=1.442_wp)
   case(p_tpssh_df)
      param = d3_param(rs6=1.223_wp, s8=1.219_wp)
   case(p_pwb6k_df)
      param = d3_param(rs6=1.660_wp, s8=0.550_wp)
   case(p_b1b95_df)
      param = d3_param(rs6=1.613_wp, s8=1.868_wp)
   case(p_bop_df)
      param = d3_param(rs6=0.929_wp, s8=1.975_wp)
   case(p_olyp_df)
      param = d3_param(rs6=0.806_wp, s8=1.764_wp)
   case(p_opbe_df)
      param = d3_param(rs6=0.837_wp, s8=2.055_wp)
   case(p_ssb_df)
      param = d3_param(rs6=1.215_wp, s8=0.663_wp)
   case(p_revssb_df)
      param = d3_param(rs6=1.221_wp, s8=0.560_wp)
   case(p_otpss_df)
      param = d3_param(rs6=1.128_wp, s8=1.494_wp)
   case(p_b3pw91_df)
      param = d3_param(rs6=1.176_wp, s8=1.775_wp)
   case(p_revpbe0_df)
      param = d3_param(rs6=0.949_wp, s8=0.792_wp)
   case(p_pbe38_df)
      param = d3_param(rs6=1.333_wp, s8=0.998_wp)
   case(p_mpw1b95_df)
      param = d3_param(rs6=1.605_wp, s8=1.118_wp)
   case(p_mpwb1k_df)
      param = d3_param(rs6=1.671_wp, s8=1.061_wp)
   case(p_bmk_df)
      param = d3_param(rs6=1.931_wp, s8=2.168_wp)
   case(p_camb3lyp_df)
      param = d3_param(rs6=1.378_wp, s8=1.217_wp)
   case(p_lcwpbe_df)
      param = d3_param(rs6=1.355_wp, s8=1.279_wp)
   case(p_m05_df)
      param = d3_param(rs6=1.373_wp, s8=0.595_wp)
   case(p_m052x_df)
      param = d3_param(rs6=1.417_wp, s8=0.000_wp)
   case(p_m06l_df)
      param = d3_param(rs6=1.581_wp, s8=0.000_wp)
   case(p_m06_df)
      param = d3_param(rs6=1.325_wp, s8=0.000_wp)
   case(p_m062x_df)
      param = d3_param(rs6=1.619_wp, s8=0.000_wp)
   case(p_m06hf_df)
      param = d3_param(rs6=1.446_wp, s8=0.000_wp)
   case(p_hcth120_df)
      param = d3_param(rs6=1.221_wp, s8=1.206_wp)
   case(p_scan_df)
      param = d3_param(rs6=1.324_wp, s8=0.000_wp)
   case(p_wb97x_df)
      param = d3_param(rs6=1.281_wp, s8=1.0_wp, rs8=1.094_wp)
   case(p_pw1pw_df)
      param = d3_param(rs6=1.4968_wp, s8=1.1786_wp)
   case(p_pbehpbe_df)
      param = d3_param(rs6=1.5703_wp, s8=1.4010_wp)
   case(p_xlyp_df)
      param = d3_param(rs6=0.9384_wp, s8=0.7447_wp)
   case(p_mpwpw_df)
      param = d3_param(rs6=1.3725_wp, s8=1.9467_wp)
   case(p_hcth407_df)
      param = d3_param(rs6=4.0426_wp, s8=2.7694_wp)
   case(p_revtpss_df)
      param = d3_param(rs6=1.3491_wp, s8=1.3666_wp)
   case(p_tauhcth_df)
      param = d3_param(rs6=0.9320_wp, s8=0.5662_wp)
   case(p_b3p_df)
      param = d3_param(rs6=1.1897_wp, s8=1.1961_wp)
   case(p_b1p_df)
      param = d3_param(rs6=1.1815_wp, s8=1.1209_wp)
   case(p_b1lyp_df)
      param = d3_param(rs6=1.3725_wp, s8=1.9467_wp)
   case(p_mpw1lyp_df)
      param = d3_param(rs6=2.0512_wp, s8=1.9529_wp)
   case(p_mpw1pw_df)
      param = d3_param(rs6=1.2892_wp, s8=1.4758_wp)
   case(p_mpw1kcis_df)
      param = d3_param(rs6=1.7231_wp, s8=2.2917_wp)
   case(p_mpwkcis1k_df)
      param = d3_param(rs6=1.4853_wp, s8=1.7553_wp)
   case(p_pbeh1pbe_df)
      param = d3_param(rs6=1.3719_wp, s8=1.0430_wp)
   case(p_pbe1kcis_df)
      param = d3_param(rs6=3.6355_wp, s8=1.7934_wp)
   case(p_x3lyp_df)
      param = d3_param(rs6=1.0_wp, s8=0.2990_wp)
   case(p_o3lyp_df)
      param = d3_param(rs6=1.4060_wp, s8=1.8058_wp)
   case(p_b97_1_df)
      param = d3_param(rs6=3.7924_wp, s8=1.6418_wp)
   case(p_b97_2_df)
      param = d3_param(rs6=1.7066_wp, s8=1.6418_wp)
   case(p_b98_df)
      param = d3_param(rs6=2.6895_wp, s8=1.9078_wp)
   case(p_hiss_df)
      param = d3_param(rs6=1.3338_wp, s8=0.7615_wp)
   case(p_hse03_df)
      param = d3_param(rs6=1.3944_wp, s8=1.0156_wp)
   case(p_revtpssh_df)
      param = d3_param(rs6=1.3224_wp, s8=1.2504_wp)
   case(p_revtpss0_df)
      param = d3_param(rs6=1.2881_wp, s8=1.0649_wp)
   case(p_tpss1kcis_df)
      param = d3_param(rs6=1.7729_wp, s8=2.0902_wp)
   case(p_tauhcthhyb_df)
      param = d3_param(rs6=1.5001_wp, s8=1.6302_wp)
   case(p_pkzb_df)
      param = d3_param(rs6=0.6327_wp, s8=0.0_wp)
   case(p_n12_df)
      param = d3_param(rs6=1.3493_wp, s8=2.3916_wp)
   case(p_mpw2plyp_df)
      param = d3_param(s6=0.66_wp, rs6=1.5527_wp, s8=0.7529_wp)
   case(p_m08hx_df)
      param = d3_param(rs6=1.6247_wp, s8=0.0_wp)
   case(p_m11l_df)
      param = d3_param(rs6=2.3933_wp, s8=1.1129_wp)
   case(p_mn15l_df)
      param = d3_param(rs6=3.3388_wp, s8=0.0_wp)
   case(p_pwp_df)
      param = d3_param(rs6=2.1040_wp, s8=0.8747_wp)
   end select

   if (present(s9)) then
      param%s9 = s9
   end if

end subroutine get_zero_damping


subroutine get_mrational_damping(param, method, error, s9)

   !> Loaded parameter record
   type(d3_param), intent(out) :: param

   !> Name of the method to look up
   character(len=*), intent(in) :: method

   !> Overwrite s9
   real(wp), intent(in), optional :: s9

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(get_method_id(method))
   case default
      call fatal_error(error, "No entry for '"//method//"' present")
      return
   case(p_b2plyp_df)
      param = d3_param(a1=0.486434_wp, s8=0.672820_wp, a2=3.656466_wp, s6=0.640000_wp)
   case(p_b3lyp_df)
      param = d3_param(a1=0.278672_wp, s8=1.466677_wp, a2=4.606311_wp)
   case(p_b97d_df)
      param = d3_param(a1=0.240184_wp, s8=1.206988_wp, a2=3.864426_wp)
   case(p_blyp_df)
      param = d3_param(a1=0.448486_wp, s8=1.875007_wp, a2=3.610679_wp)
   case(p_bp_df)
      param = d3_param(a1=0.821850_wp, s8=3.140281_wp, a2=2.728151_wp)
   case(p_pbe_df)
      param = d3_param(a1=0.012092_wp, s8=0.358940_wp, a2=5.938951_wp)
   case(p_pbe0_df)
      param = d3_param(a1=0.007912_wp, s8=0.528823_wp, a2=6.162326_wp)
   case(p_lcwpbe_df)
      param = d3_param(a1=0.563761_wp, s8=0.906564_wp, a2=3.593680_wp)
   end select

   if (present(s9)) then
      param%s9 = s9
   end if

end subroutine get_mrational_damping


subroutine get_mzero_damping(param, method, error, s9)

   !> Loaded parameter record
   type(d3_param), intent(out) :: param

   !> Name of the method to look up
   character(len=*), intent(in) :: method

   !> Overwrite s9
   real(wp), intent(in), optional :: s9

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(get_method_id(method))
   case default
      call fatal_error(error, "No entry for '"//method//"' present")
      return
   case(p_b2plyp_df)
      param = d3_param(rs6=1.313134_wp, s8=0.717543_wp, bet=0.016035_wp, s6=0.640000_wp)
   case(p_b3lyp_df)
      param = d3_param(rs6=1.338153_wp, s8=1.532981_wp, bet=0.013988_wp)
   case(p_b97d_df)
      param = d3_param(rs6=1.151808_wp, s8=1.020078_wp, bet=0.035964_wp)
   case(p_blyp_df)
      param = d3_param(rs6=1.279637_wp, s8=1.841686_wp, bet=0.014370_wp)
   case(p_bp_df)
      param = d3_param(rs6=1.233460_wp, s8=1.945174_wp, bet=0.000000_wp)
   case(p_pbe_df)
      param = d3_param(rs6=2.340218_wp, s8=0.000000_wp, bet=0.129434_wp)
   case(p_pbe0_df)
      param = d3_param(rs6=2.077949_wp, s8=0.000081_wp, bet=0.116755_wp)
   case(p_lcwpbe_df)
      param = d3_param(rs6=1.366361_wp, s8=1.280619_wp, bet=0.003160_wp)
   end select

   if (present(s9)) then
      param%s9 = s9
   end if

end subroutine get_mzero_damping


subroutine get_optimizedpower_damping(param, method, error, s9)

   !> Loaded parameter record
   type(d3_param), intent(out) :: param

   !> Name of the method to look up
   character(len=*), intent(in) :: method

   !> Overwrite s9
   real(wp), intent(in), optional :: s9

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(get_method_id(method))
   case default
      call fatal_error(error, "No entry for '"//method//"' present")
      return
   case(p_pbe_df)
      param = d3_param(s6=0.91826_wp, s8=0.0_wp, a1=0.200_wp, a2=4.750_wp, bet=6.0_wp)
   case(p_pbe0_df)
      param = d3_param(s6=0.88290_wp, s8=0.0_wp, a1=0.150_wp, a2=4.750_wp, bet=6.0_wp)
   case(p_revtpss_df)
      param = d3_param(s6=1.0_wp, s8=0.27632_wp, a1=0.700_wp, a2=2.500_wp, bet=8.0_wp)
   case(p_revtpssh_df)
      param = d3_param(s6=1.0_wp, s8=0.12467_wp, a1=0.575_wp, a2=3.000_wp, bet=10.0_wp)
   case(p_blyp_df)
      param = d3_param(s6=1.0_wp, s8=1.31867_wp, a1=0.425_wp, a2=3.50_wp, bet=2.0_wp)
   case(p_b3lyp_df)
      param = d3_param(s6=1.0_wp, s8=0.78311_wp, a1=0.300_wp, a2=4.25_wp, bet=4.0_wp)
   case(p_b97d_df)
      param = d3_param(s6=1.0_wp, s8=1.46861_wp, a1=0.600_wp, a2=2.50_wp, bet=0.0_wp)
   case(p_b97_1_df, p_b97_2_df)
      param = d3_param(s6=0.97388_wp, s8=0.0_wp, a1=0.150_wp, a2=4.25_wp, bet=6.0_wp)
   case(p_revpbe_df)
      param = d3_param(s6=1.0_wp, s8=1.44765_wp, a1=0.600_wp, a2=2.50_wp, bet=0.0_wp)
   case(p_revpbe0_df)
      param = d3_param(s6=1.0_wp, s8=1.25684_wp, a1=0.725_wp, a2=2.25_wp, bet=0.0_wp)
   case(p_tpss_df)
      param = d3_param(s6=1.0_wp, s8=0.51581_wp, a1=0.575_wp, a2=3.00_wp, bet=8.0_wp)
   case(p_tpssh_df)
      param = d3_param(s6=1.0_wp, s8=0.43185_wp, a1=0.575_wp, a2=3.00_wp, bet=8.0_wp)
   case(p_ms2_df)
      param = d3_param(s6=1.0_wp, s8=0.90743_wp, a1=0.700_wp, a2=4.00_wp, bet=2.0_wp)
   case(p_ms2h_df)
      param = d3_param(s6=1.0_wp, s8=1.69464_wp, a1=0.650_wp, a2=4.75_wp, bet=0.0_wp)
   end select

   if (present(s9)) then
      param%s9 = s9
   end if

end subroutine get_optimizedpower_damping


!> Convert string to lower case
pure function lowercase(str) result(lcstr)
   character(len=*), intent(in)  :: str
   character(len=len_trim(str)) :: lcstr
   integer :: ilen, ioffset, iquote, i, iav, iqc

   ilen=len_trim(str)
   ioffset=iachar('A')-iachar('a')
   iquote=0
   lcstr=str
   do i=1, ilen
      iav=iachar(str(i:i))
      if(iquote==0 .and. (iav==34 .or.iav==39)) then
         iquote=1
         iqc=iav
        cycle
      endif
      if(iquote==1 .and. iav==iqc) then
         iquote=0
         cycle
      endif
      if (iquote==1) cycle
      if(iav >= iachar('A') .and. iav <= iachar('Z')) then
         lcstr(i:i)=achar(iav-ioffset)
      else
         lcstr(i:i)=str(i:i)
      endif
   enddo

end function lowercase


end module dftd3_param
