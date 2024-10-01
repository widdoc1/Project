module qcd
  use types; use consts_dp

  integer,  parameter :: nf_def = 5
  real(dp), parameter :: ca_def = 3, cf_def = four/three, tr_def = half
  real(dp), parameter :: tf_def = tr_def * nf_def

  !-- the following are all modifiable, but have default values
  real(dp), public :: ca = ca_def
  real(dp), public :: cf = cf_def
  real(dp), public :: tr = tr_def

  ! beta2 is from Tarasov JINR P2-82-900 (Dubna 1982)
  ! Larin & Vermaseren NIKHEF-H/92-17 hep-ph/9302208
  real(dp), public :: beta0 = (11*ca_def - four*tf_def)/(12*pi)
  real(dp), public :: twopi_beta0 = (11*ca_def - four*tf_def)/6.0_dp
  real(dp), public :: beta1 = (17*ca_def**2 - tf_def*(10*ca_def + 6*cf_def))&
       &                        / (24*pisq)
  real(dp), public :: beta2 = (2857*ca_def**3  &
       & + (54*cf_def**2-615*cf_def*ca_def-1415*ca_def**2)*two*tf_def&
       & + (66*cf_def + 79*ca_def)*(two*tf_def)**2)/(3456*pi**3)


  real(dp), public :: cmw_K = ca_def*(67.0_dp/18.0_dp - pisq/6.0_dp) &
       &                      - tf_def * 10.0_dp/9.0_dp
  
  !!!! Taken from Moch, Vermaseren & Vogt: NB TR dependence not in...
  real(dp), parameter :: cmw_K2_def = &
       &  ca_def**2 * ( 245._dp/24._dp - 67._dp/9._dp*zeta2 &
       &             + 11.0_dp/6._dp * zeta3 + 11.0_dp/5._dp * zeta2**2)&
       &+ two*cf_def*tf_def * (-55._dp/24._dp + 2*zeta3)&
       &+ two*ca_def*tf_def * (-209._dp/108._dp + 10._dp/9._dp*zeta2 &
       &                   - 7._dp/3._dp * zeta3)&
       &- four*tf_def**2 /27._dp
  real(dp), public :: cmw_K2  = cmw_K2_def
  real(dp), public :: mvv_A3  = 16*cf_def*cmw_K2_def
  real(dp), public :: mvv_A3G = 16*ca_def*cmw_K2_def
  
end module qcd
