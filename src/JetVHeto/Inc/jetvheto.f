!--- include file for controlling JetVHeto resummation variables

!--- the following are set by the input file
      logical  :: jetvheto ! whether resummation process should be used
      logical  :: pure_lumi ! whether sudakov reweighting should be used
      real(dp) :: ptj_veto ! pt(jet) < pt(jet_veto) scale
      real(dp) :: q_scalestart ! JetVHeto scale in input file

!--- the follow are parameters that should not generally be changed
      real(dp), parameter :: p_pow = 5._dp ! exponent for power supression of resummed terms

!--- the following is used when computing JetVHeto processes
      real(dp) :: q_scale ! JetVHeto scale
      real(dp) :: L_tilde ! L = log(pt/Q)

! Should I include R scales for small R resummation?

      common/jetvhetoinputs/ptj_veto,q_scalestart,jetvheto,pure_lumi
      common/jetvhetovars/q_scale,L_tilde
!$omp threadprivate(/jetvhetovars/)
