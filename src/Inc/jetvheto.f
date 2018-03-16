!--- include file for controlling JetVHeto resummation variables

!--- the following are set by the input file
      logical  :: jetvheto ! whether a resummation process is being performed
      real(dp) :: ptj_veto ! pt(jet) < pt(jet_veto) scale
      real(dp) :: q_scalestart ! JetVHeto scale in input file
!      real(dp) :: r_scalestart ! R scale in input file
      logical  :: do_lumi ! do luminosity reweight
      logical  :: do_suda ! do sudakov reweight

!--- the follow are parameters that should not generally be changed
      character(len=*), parameter :: observable = 'pt_jet' ! maybe add experimental LL_R resummation later?
      real(dp), parameter :: p_pow = 5._dp ! exponent for power supression of resummed terms

!--- the following is set by MCFM and used when computing JetVHeto
!--- processes
      real(dp) :: q_scale ! JetVHeto scale
!      real(dp) :: r_scale ! R scale
      real(dp) :: L_tilde ! L = log(pt/Q)
!---  The configuration of the Born level hard process for JetVHeto 
!---  resummation 'DY' for qqb initiated or 'H' for gg initiated
      character(len=2) :: born_config 


      
! Should I include R scales for small R resummation?

      common/jetvhetoinputs/ptj_veto,q_scalestart,jetvheto
      common/jetvhetoinputs/do_lumi,do_suda
      common/jetvhetovars/q_scale,L_tilde,born_config
!$omp threadprivate(/jetvhetovars/)
