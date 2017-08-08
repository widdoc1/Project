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
      character(len=10) :: robs
      character(len=1) :: match_scheme  ! multiplicative or additive?

      common/resminputs/Q_scalestart,R_scalestart,resum,robs,
     &     match_scheme
      common/resmscale/Q_scale,R_scale
!$omp threadprivate(/resmscale/)
