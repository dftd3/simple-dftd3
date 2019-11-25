module d3par_damping_parameters
   use iso_fortran_env, only: wp => real64
   use d3def_damping_parameters, only: d3par => d3_damping_parameters
   implicit none
   public :: get_d3_damping_parameters
   private

   enum, bind(c)
      enumerator :: &
         &  p_df_none           =   0,  &
         &  p_df_hf             =   1,  &
         &  p_df_blyp           =   2,  &
         &  p_df_bpbe           =   3,  &
         &  p_df_bp             =   4,  &
         &  p_df_bpw            =   5,  &
         &  p_df_lb94           =   6,  &
         &  p_df_mpwlyp         =   7,  &
         &  p_df_mpwpw          =   8,  &
         &  p_df_olyp           =   9,  &
         &  p_df_opbe           =  10,  &
         &  p_df_pbe            =  11,  &
         &  p_df_rpbe           =  12,  &
         &  p_df_revpbe         =  13,  &
         &  p_df_pw86pbe        =  14,  &
         &  p_df_rpw86pbe       =  15,  &
         &  p_df_pw91           =  16,  &
         &  p_df_pwp            =  17,  &
         &  p_df_xlyp           =  18,  &
         &  p_df_b97            =  19,  &
         &  p_df_tpss           =  20,  &
         &  p_df_revtpss        =  21,  &
         &  p_df_scan           =  22,  &
         &  p_df_b1lyp          =  23,  &
         &  p_df_b3lyp          =  24,  &
         &  p_df_bhlyp          =  25,  &
         &  p_df_b1p            =  26,  &
         &  p_df_b3p            =  27,  &
         &  p_df_b1pw           =  28,  &
         &  p_df_b3pw           =  29,  &
         &  p_df_o3lyp          =  30,  &
         &  p_df_revpbe0        =  31,  &
         &  p_df_revpbe38       =  32,  &
         &  p_df_pbe0           =  33,  &
         &  p_df_pwp1           =  34,  &
         &  p_df_pw1pw          =  35,  &
         &  p_df_mpw1pw         =  36,  &
         &  p_df_mpw1lyp        =  37,  &
         &  p_df_pw6b95         =  38,  &
         &  p_df_tpssh          =  39,  &
         &  p_df_tpss0          =  40,  &
         &  p_df_x3lyp          =  41,  &
         &  p_df_m06l           =  42,  &
         &  p_df_m06            =  43,  &
         &  p_df_m062x          =  44,  &
         &  p_df_wb97           =  45,  &
         &  p_df_wb97x          =  46,  &
         &  p_df_camb3lyp       =  47,  &
         &  p_df_lcblyp         =  48,  &
         &  p_df_lh07tsvwn      =  49,  &
         &  p_df_lh07ssvwn      =  50,  &
         &  p_df_lh12ctssirpw92 =  51,  &
         &  p_df_lh12ctssifpw92 =  52,  &
         &  p_df_lh14tcalpbe    =  53,  &
         &  p_df_b2plyp         =  54,  &
         &  p_df_b2gpplyp       =  55,  &
         &  p_df_mpw2plyp       =  56,  &
         &  p_df_pwpb95         =  57,  &
         &  p_df_dsdblyp        =  58,  &
         &  p_df_dsdpbe         =  59,  &
         &  p_df_dsdpbeb95      =  60,  &
         &  p_df_dsdpbep86      =  61,  &
         &  p_df_dsdsvwn        =  62,  &
         &  p_df_dodblyp        =  63,  &
         &  p_df_dodpbe         =  64,  &
         &  p_df_dodpbeb95      =  65,  &
         &  p_df_dodpbep86      =  66,  &
         &  p_df_dodsvwn        =  67,  &
         &  p_df_pbe0_2         =  68,  &
         &  p_df_pbe0_dh        =  69,  &
         &  p_df_hf3c           =  70,  &
         &  p_df_hf3cv          =  71,  &
         &  p_df_pbeh3c         =  72,  &
         &  p_df_b973c          =  73,  &
         &  p_df_hsesol         =  74,  &
         &  p_df_pwgga          =  75,  &
         &  p_df_dftb3          =  76,  &
         &  p_df_hcth120        =  77,  &
         &  p_df_ptpss          =  78,  &
         &  p_df_lcwpbe         =  79,  &
         &  p_df_bmk            =  80,  &
         &  p_df_b1b95          =  81,  &
         &  p_df_pwb6k          =  82,  &
         &  p_df_otpss          =  83,  &
         &  p_df_ssb            =  84,  &
         &  p_df_revssb         =  85,  &
         &  p_df_pbesol         =  86,  &
         &  p_df_hse06          =  87,  &
         &  p_df_pbexalpha      =  88,  &
         &  p_df_pbehpbe        =  89,  &
         &  p_df_hcth407        =  90,  &
         &  p_df_n12            =  91,  &
         &  p_df_pkzb           =  92,  &
         &  p_df_thcth          =  93,  &
         &  p_df_m11l           =  94,  &
         &  p_df_mn15l          =  95,  &
         &  p_df_mpwb1k         =  96,  &
         &  p_df_mpw1kcis       =  97,  &
         &  p_df_mpwkcis1k      =  98,  &
         &  p_df_pbeh1pbe       =  99,  &
         &  p_df_pbe1kcis       = 100,  &
         &  p_df_b97_1          = 101,  &
         &  p_df_b97_2          = 102,  &
         &  p_df_b98            = 103,  &
         &  p_df_hiss           = 104,  &
         &  p_df_hse03          = 105,  &
         &  p_df_revtpssh       = 106,  &
         &  p_df_tpss1kcis      = 107,  &
         &  p_df_m05            = 108,  &
         &  p_df_m052x          = 109,  &
         &  p_df_m08hx          = 110,  &
         &  p_df_lcwhpbe        = 111,  &
         &  p_df_mn12l          = 112,  &
         &  p_df_tauhcthhyb     = 113,  &
         &  p_df_sogga11x       = 114,  &
         &  p_df_n12sx          = 115,  &
         &  p_df_mn12sx         = 116,  &
         &  p_df_mn15           = 117,  &
         &  p_df_glyp           = 118,  &
         &  p_df_bop            = 119,  &
         &  p_df_mpw1b95        = 120,  &
         &  p_df_revpbe0dh      = 121,  &
         &  p_df_revtpss0       = 122
   end enum
   integer, parameter :: ik = kind(p_df_none)

contains

pure elemental function func_to_enum(func) result(enum)
   integer(ik) :: enum
   character(len=*),intent(in) :: func
   select case(func)
   case default;                             enum = p_df_none
   case('hf');                               enum = p_df_hf
   case('b-lyp','blyp');                     enum = p_df_blyp
   case('bpbe');                             enum = p_df_bpbe
   case('b-p','bp86','bp','b-p86');          enum = p_df_bp
   case('bpw','b-pw');                       enum = p_df_bpw
   case('lb94');                             enum = p_df_lb94
   case('mpwlyp','mpw-lyp');                 enum = p_df_mpwlyp
   case('mpwpw','mpw-pw','mpwpw91');         enum = p_df_mpwpw
   case('o-lyp','olyp');                     enum = p_df_olyp
   case('opbe');                             enum = p_df_opbe
   case('pbe');                              enum = p_df_pbe
   case('rpbe');                             enum = p_df_rpbe
   case('revpbe');                           enum = p_df_revpbe
   case('pw86pbe');                          enum = p_df_pw86pbe
   case('rpw86pbe');                         enum = p_df_rpw86pbe
   case('pw91');                             enum = p_df_pw91
   case('pwp','pw-p','pw91p86');             enum = p_df_pwp
   case('x-lyp','xlyp');                     enum = p_df_xlyp
   case('b97');                              enum = p_df_b97
   case('tpss');                             enum = p_df_tpss
   case('revtpss');                          enum = p_df_revtpss
   case('scan');                             enum = p_df_scan
   case('b1lyp','b1-lyp');                   enum = p_df_b1lyp
   case('b3-lyp','b3lyp');                   enum = p_df_b3lyp
   case('bh-lyp','bhlyp');                   enum = p_df_bhlyp
   case('b1p','b1-p','b1p86');               enum = p_df_b1p
   case('b3p','b3-p','b3p86');               enum = p_df_b3p
   case('b1pw','b1-pw','b1pw91');            enum = p_df_b1pw
   case('b3pw','b3-pw','b3pw91');            enum = p_df_b3pw
   case('o3-lyp','o3lyp');                   enum = p_df_o3lyp
   case('revpbe0');                          enum = p_df_revpbe0
   case('revpbe38');                         enum = p_df_revpbe38
   case('pbe0');                             enum = p_df_pbe0
   case('pwp1');                             enum = p_df_pwp1
   case('pw1pw','pw1-pw');                   enum = p_df_pw1pw
   case('mpw1pw','mpw1-pw','mpw1pw91');      enum = p_df_mpw1pw
   case('mpw1lyp','mpw1-lyp');               enum = p_df_mpw1lyp
   case('pw6b95');                           enum = p_df_pw6b95
   case('tpssh');                            enum = p_df_tpssh
   case('tpss0');                            enum = p_df_tpss0
   case('x3-lyp','x3lyp');                   enum = p_df_x3lyp
   case('m06l');                             enum = p_df_m06l
   case('m06');                              enum = p_df_m06
   case('m06-2x','m062x');                   enum = p_df_m062x
   case('wb97','ωb97','omegab97');           enum = p_df_wb97
   case('wb97x','ωb97x','omegab97x');        enum = p_df_wb97x
   case('cam-b3lyp');                        enum = p_df_camb3lyp
   case('lc-blyp');                          enum = p_df_lcblyp
   case('lh07tsvwn','lh07t-svwn');           enum = p_df_lh07tsvwn
   case('lh07ssvwn','lh07s-svwn');           enum = p_df_lh07ssvwn
   case('lh12ctssirpw92','lh12ct-ssirpw92'); enum = p_df_lh12ctssirpw92
   case('lh12ctssifpw92','lh12ct-ssifpw92'); enum = p_df_lh12ctssifpw92
   case('lh14tcalpbe','lh14t-calpbe');       enum = p_df_lh14tcalpbe
   case('b2plyp','b2-plyp');                 enum = p_df_b2plyp
   case('b2gpplyp','b2gp-plyp');             enum = p_df_b2gpplyp
   case('mpw2plyp');                         enum = p_df_mpw2plyp
   case('pwpb95');                           enum = p_df_pwpb95
   case('dsdblyp','dsd-blyp');               enum = p_df_dsdblyp
   case('dsdpbe','dsd-pbe');                 enum = p_df_dsdpbe
   case('dsdpbeb95','dsd-pbeb95');           enum = p_df_dsdpbeb95
   case('dsdpbep86','dsd-pbep86');           enum = p_df_dsdpbep86
   case('dsdsvwn','dsd-svwn');               enum = p_df_dsdsvwn
   case('dodblyp','dod-blyp');               enum = p_df_dodblyp
   case('dodpbe','dod-pbe');                 enum = p_df_dodpbe
   case('dodpbeb95','dod-pbeb95');           enum = p_df_dodpbeb95
   case('dodpbep86','dod-pbep86');           enum = p_df_dodpbep86
   case('dodsvwn','dod-svwn');               enum = p_df_dodsvwn
   case('pbe02','pbe0-2');                   enum = p_df_pbe0_2
   case('pbe0dh','pbe0-dh');                 enum = p_df_pbe0_dh
   case('hf-3c','hf3c');                     enum = p_df_hf3c
   case('hf-3cv','hf3cv');                   enum = p_df_hf3cv
   case('pbeh3c','pbeh-3c');                 enum = p_df_pbeh3c
   case('b973c','b97-3c');                   enum = p_df_b973c
   case('hsesol');                           enum = p_df_hsesol
   case('pwgga');                            enum = p_df_pwgga
   case('dftb3');                            enum = p_df_dftb3
   case('hcth120');                          enum = p_df_hcth120
   case('ptpss');                            enum = p_df_ptpss
   case('lc-wpbe','lcwpbe');                 enum = p_df_lcwpbe
   case('bmk');                              enum = p_df_bmk
   case('b1b95');                            enum = p_df_b1b95
   case('bwb6k');                            enum = p_df_pwb6k
   case('otpss');                            enum = p_df_otpss
   case('ssb');                              enum = p_df_ssb
   case('revssb');                           enum = p_df_revssb
   case('pbesol');                           enum = p_df_pbesol
   case('hse06');                            enum = p_df_hse06
   case('pbexalpha');                        enum = p_df_pbexalpha
   case('pbehpbe');                          enum = p_df_pbehpbe
   case('hcth407');                          enum = p_df_hcth407
   case('n12');                              enum = p_df_n12
   case('pkzb');                             enum = p_df_pkzb
   case('thcth','tauhctc');                  enum = p_df_thcth
   case('m11l');                             enum = p_df_m11l
   case('mn15l');                            enum = p_df_mn15l
   case('mpwb1k');                           enum = p_df_mpwb1k
   case('mpw1kcis');                         enum = p_df_mpw1kcis
   case('mpwkcis1k');                        enum = p_df_mpwkcis1k
   case('pbeh1pbe');                         enum = p_df_pbeh1pbe
   case('pbe1kcis');                         enum = p_df_pbe1kcis
   case('b97-1');                            enum = p_df_b97_1
   case('b97-2');                            enum = p_df_b97_2
   case('b98');                              enum = p_df_b98
   case('hiss');                             enum = p_df_hiss
   case('hse03');                            enum = p_df_hse03
   case('revtpssh');                         enum = p_df_revtpssh
   case('tpss1kcis');                        enum = p_df_tpss1kcis
   case('m05');                              enum = p_df_m05
   case('m052x','m05-2x');                   enum = p_df_m052x
   case('m08hx','m08-hx');                   enum = p_df_m08hx
   case('lcwhpbe','lc-whpbe');               enum = p_df_lcwhpbe
   case('mn12l');                            enum = p_df_mn12l
   case('tauhcthhyb');                       enum = p_df_tauhcthhyb
   case('sogga11x');                         enum = p_df_sogga11x
   case('n12sx');                            enum = p_df_n12sx
   case('mn12sx');                           enum = p_df_mn12sx
   case('mn15');                             enum = p_df_mn15
   case('glyp','g-lyp');                     enum = p_df_glyp
   case('revpbe0dh','revpbe0-dh');           enum = p_df_revpbe0dh
   case('revtpss0');                         enum = p_df_revtpss0
   end select
end function func_to_enum

subroutine get_d3_damping_parameters(param, func, found)
   type(d3par), intent(out) :: param
   character(len=*), intent(in) :: func
   logical, intent(out), optional :: found
   integer(ik) :: num

   if (present(found)) found = .false.

   num = func_to_enum(func)

   select case(num)
   case default; return
   case (p_df_bp); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.3946_wp, s8 =3.2822_wp, a2=4.8516_wp )
   case(p_df_blyp); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4298_wp, s8 =2.6996_wp, a2=4.2359_wp )
   case(p_df_revpbe); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.5238_wp, s8 =2.3550_wp, a2=3.5016_wp )
   case(p_df_rpbe); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.1820_wp, s8 =0.8318_wp, a2=4.0094_wp )
   case(p_df_b97); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.5545_wp, s8 =2.2609_wp, a2=3.2297_wp )
   case(p_df_pbe); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4289_wp, s8 =0.7875_wp, a2=4.4407_wp )
   case(p_df_rpw86pbe); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4613_wp, s8 =1.3845_wp, a2=4.5062_wp )
   case(p_df_b3lyp); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.3981_wp, s8 =1.9889_wp, a2=4.4211_wp )
   case(p_df_tpss); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4535_wp, s8 =1.9435_wp, a2=4.4752_wp )
   case(p_df_hf); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.3385_wp, s8 =0.9171_wp, a2=2.8830_wp )
   case(p_df_tpss0); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.3768_wp, s8 =1.2576_wp, a2=4.5865_wp )
   case(p_df_pbe0); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4145_wp, s8 =1.2177_wp, a2=4.8593_wp )
   case(p_df_hse06); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.383_wp, s8 =2.310_wp, a2=5.685_wp )
   case(p_df_revpbe38); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4309_wp, s8 =1.4760_wp, a2=3.9446_wp )
   case(p_df_pw6b95); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.2076_wp, s8 =0.7257_wp, a2=6.3750_wp )
   case(p_df_b2plyp); param = d3par ( &
   &  s6=0.6400_wp, a1 =0.3065_wp, s8 =0.9147_wp, a2=5.0570_wp )
   case(p_df_dsdblyp); param = d3par ( &
   &  s6=0.5000_wp, a1 =0.0000_wp, s8 =0.2130_wp, a2=6.0519_wp )
   !case(p_df_dsdblypfc); param = d3par ( &
   !&  s6=0.5000_wp, a1 =0.0009_wp, s8 =0.2112_wp, a2=5.9807_wp )
   case(p_df_bop); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4870_wp, s8 =3.2950_wp, a2=3.5043_wp )
   case(p_df_mpwlyp); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4831_wp, s8 =2.0077_wp, a2=4.5323_wp )
   case(p_df_olyp); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.5299_wp, s8 =2.6205_wp, a2=2.8065_wp )
   case(p_df_pbesol); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4466_wp, s8 =2.9491_wp, a2=6.1742_wp )
   case(p_df_bpbe); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4567_wp, s8 =4.0728_wp, a2=4.3908_wp )
   case(p_df_opbe); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.5512_wp, s8 =3.3816_wp, a2=2.9444_wp )
   case(p_df_ssb); param = d3par ( &
   &  s6=1.0000_wp, a1 =-.0952_wp, s8 =-.1744_wp, a2=5.2170_wp )
   case(p_df_revssb); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4720_wp, s8 =0.4389_wp, a2=4.0986_wp )
   case(p_df_otpss); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4634_wp, s8 =2.7495_wp, a2=4.3153_wp )
   case(p_df_b3pw); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4312_wp, s8 =2.8524_wp, a2=4.4693_wp )
   case(p_df_bhlyp); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.2793_wp, s8 =1.0354_wp, a2=4.9615_wp )
   case(p_df_revpbe0); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4679_wp, s8 =1.7588_wp, a2=3.7619_wp )
   case(p_df_tpssh); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4529_wp, s8 =2.2382_wp, a2=4.6550_wp )
   case(p_df_mpw1b95); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.1955_wp, s8 =1.0508_wp, a2=6.4177_wp )
   case(p_df_pwb6k); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.1805_wp, s8 =0.9383_wp, a2=7.7627_wp )
   case(p_df_b1b95); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.2092_wp, s8 =1.4507_wp, a2=5.5545_wp )
   case(p_df_bmk); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.1940_wp, s8 =2.0860_wp, a2=5.9197_wp )
   case(p_df_camb3lyp); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.3708_wp, s8 =2.0674_wp, a2=5.4743_wp )
   case(p_df_lcwpbe); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.3919_wp, s8 =1.8541_wp, a2=5.0897_wp )
   case(p_df_b2gpplyp); param = d3par ( &
   &  s6=0.5600_wp, a1 =0.0000_wp, s8 =0.2597_wp, a2=6.3332_wp )
   case(p_df_ptpss); param = d3par ( &
   &  s6=0.7500_wp, a1 =0.0000_wp, s8 =0.2804_wp, a2=6.5745_wp )
   case(p_df_pwpb95); param = d3par ( &
   &  s6=0.8200_wp, a1 =0.0000_wp, s8 =0.2904_wp, a2=7.3141_wp )
   case(p_df_hcth120); param = d3par ( &
   &  s6=1.0000_wp, a1=0.3563,   s8=1.0821,   a2=4.3359_wp )
   case(p_df_dftb3); param = d3par ( &
   &  s6=1.0000_wp, a1=0.5719_wp, s8=0.5883_wp, a2=3.6017_wp )
   case(p_df_pw1pw); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.3807_wp, s8 =2.3363_wp, a2=5.8844_wp )
   case(p_df_pwgga); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.2211_wp, s8 =2.6910_wp, a2=6.7278_wp )
   case(p_df_hsesol); param = d3par ( &
   &  s6=1.0000_wp, a1 =0.4650_wp, s8 =2.9215_wp, a2=6.2003_wp )
!  special HF-D3-gCP-SRB/MINIX parametrization
   case(p_df_hf3c); param = d3par ( &
   &  s6=1.0000_wp, a1=0.4171_wp, s8=0.8777_wp, a2=2.9149_wp )
!  special HF-D3-gCP-SRB2/ECP-2G parametrization
   case(p_df_hf3cv); param = d3par ( &
   &  s6=1.0000_wp, a1=0.3063_wp, s8=0.5022_wp, a2=3.9856_wp )
!  special PBEh-D3-gCP/def2-mSVP parametrization
   case(p_df_pbeh3c); param = d3par ( &
   &  s6=1.0000_wp, a1=0.4860_wp, s8=0.0000_wp, a2=4.5000_wp )
!  parameters from GMTKN55
   case(p_df_pbehpbe); param = d3par (&
   &  s6=1.0000_wp, a1=0.0000_wp, s8=1.1152_wp, a2=6.7184_wp )
   case(p_df_xlyp); param = d3par (&
   &  s6=1.0000_wp, a1=0.0809_wp, s8=2.0077_wp, a2=4.5323_wp )
   case(p_df_mpwpw); param = d3par (&
   &  s6=1.0000_wp, a1=0.3168_wp, s8=1.7074_wp, a2=4.7732_wp )
   case(p_df_hcth407); param = d3par (&
   &  s6=1.0000_wp, a1=0.0000_wp, s8=0.6490_wp, a2=4.8162_wp )
   case(p_df_revtpss); param = d3par (&
   &  s6=1.0000_wp, a1=0.4426_wp, s8=1.4023_wp, a2=4.4723_wp )
   case(p_df_thcth); param = d3par (&
   &  s6=1.0000_wp, a1=0.0000_wp, s8=1.2626_wp, a2=5.6162_wp )
   case(p_df_b3p); param = d3par (&
   &  s6=1.0000_wp, a1=0.4601_wp, s8=3.3211_wp, a2=4.9294_wp )
   case(p_df_b1p); param = d3par (&
   &  s6=1.0000_wp, a1=0.4724_wp, s8=3.5681_wp, a2=4.9858_wp )
   case(p_df_b1lyp); param = d3par (&
   &  s6=1.0000_wp, a1=0.1986_wp, s8=2.1167_wp, a2=5.3875_wp )
   case(p_df_mpw1pw); param = d3par (&
   &  s6=1.0000_wp, a1=0.3342_wp, s8=1.8744_wp, a2=4.9819_wp )
   case(p_df_mpw1kcis); param = d3par (&
   &  s6=1.0000_wp, a1=0.0576_wp, s8=1.0893_wp, a2=5.5314_wp )
   case(p_df_mpwkcis1k); param = d3par (&
   &  s6=1.0000_wp, a1=0.0855_wp, s8=1.2875_wp, a2=5.8961_wp )
   case(p_df_pbeh1pbe); param = d3par (&
   &  s6=1.0000_wp, a1=0.0000_wp, s8=1.4877_wp, a2=7.0385_wp )
   case(p_df_pbe1kcis); param = d3par (&
   &  s6=1.0000_wp, a1=0.0000_wp, s8=0.7688_wp, a2=6.2794_wp )
   case(p_df_x3lyp); param = d3par (&
   &  s6=1.0000_wp, a1=0.2022_wp, s8=1.5744_wp, a2=5.4184_wp )
   case(p_df_o3lyp); param = d3par (&
   &  s6=1.0000_wp, a1=0.0963_wp, s8=1.8171_wp, a2=5.9940_wp)
   end select

   if (present(found)) found = .true.

end subroutine get_d3_damping_parameters

end module d3par_damping_parameters
