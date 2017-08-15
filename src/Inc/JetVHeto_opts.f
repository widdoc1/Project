!     container for resummation ingredients
      type(process_and_parameters) :: resm_opts
      real(dp) :: L_tilde, facscaleLtilde
      real(dp) :: sudakov

      common/resmparams/resm_opts,L_tilde,facscaleLtilde,sudakov
!$omp threadprivate(/resmparams/)
