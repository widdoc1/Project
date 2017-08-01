!     scales used in the JetVHeto resummation
      logical :: resm          ! whether resummation is to be performed

!     Boson invariant mass
      real(dp) :: M_B, M_B2

! scale variation of ptjet
      real(dp) :: Q_scale       ! resummation scale Q
      real(dp) :: Q_scalestart
! small R resummation
      real(dp) :: R_scale       ! resummation scale R
      real(dp) :: R_scalestart

! container for resummation ingredients
      type(process_and_parameters) :: resm_opts

      common/resminputs/Q_scalestart,R_scalestart
      common/resmscale/Q_scale,R_scale
      common/resmparams/resm_opts
!$omp threadprivate(/resmscale/)
