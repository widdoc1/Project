!     scales used in the JetVHeto resummation
      logical :: resum          ! whether resummation is to be performed

!     Boson invariant mass
      real(dp) :: M_B, M_B2

! scale variation of ptjet
      real(dp) :: Q_scale       ! resummation scale Q
      real(dp) :: Q_scalestart
! small R resummation
      real(dp) :: R_scale       ! resummation scale R
      real(dp) :: R_scalestart

! observable
      character(len=10) :: kobs
      character(len=2) :: Bconf ! Born level configuration qq or gg initiated
      character(len=1) :: mscheme  ! multiplicative or additive?

      logical :: pure_lumi

      common/resminputs/Q_scalestart,R_scalestart,resum,pure_lumi,kobs
      common/resminputs/Bconf,mscheme
      common/resmscale/Q_scale,R_scale
!$omp threadprivate(/resmscale/)
