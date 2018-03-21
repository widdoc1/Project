!--- include file for controlling JetVHeto resummation variables

!--- the following are set by the input file
      logical          :: jetvheto ! whether a resummation process is being performed
      logical          :: do_lumi ! do luminosity reweight
      logical          :: do_suda ! do sudakov reweight
      character(len=3) :: observable ! observable to resum
      logical          :: small_r ! should small R resummation be used
      real(dp)         :: q_scalestart ! JetVHeto scale in input file
      real(dp)         :: r_scale ! R scale in input file, unchanged by dynamics
      real(dp)         :: ptj_veto ! pt(jet_veto) scale

!--- the following are parameters that should not generally be changed
      real(dp), parameter :: p_pow = 5._dp ! exponent for power supression of resummed terms

!--- intermediate values
      character*24 :: obs_string 

!--- the following is set by MCFM and used when computing JetVHeto
!--- processes
      character(len=2) :: born_config ! the Born configuration: H or DY
      real(dp) :: q_scale ! JetVHeto scale
      real(dp) :: L_tilde ! L = log(pt/Q)

      common/jetvhetoinputs/q_scalestart,r_scale,ptj_veto
      common/jetvhetoinputs/jetvheto,do_lumi,do_suda,small_r
      common/jetvhetoinputs/observable,born_config
      common/jetvhetointerim/obs_string
      common/jetvhetovars/q_scale,L_tilde
!$omp threadprivate(/jetvhetovars/)
