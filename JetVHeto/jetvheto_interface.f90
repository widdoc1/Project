!========================================================
!--------------------------------------------------------
! Interface to access the Sudakov form factor
! -------------------------------------------------------
!========================================================

function sudakov(proc, M, muR, muF, Q, as, p, jet_radius,&
     &observable, small_r, small_r_R0, ptj_veto, order) result(res)
  use types; use consts_dp
  use rad_tools
  use resummation
  implicit none
  character(len=*),            intent(in)  :: proc, observable
  integer,                     intent(in)  :: order
  real(dp),                    intent(in)  :: M, muR, muF, Q, as, p,&
       &small_r_R0, jet_radius, ptj_veto(:)
  logical :: small_r
  type(process_and_parameters)             :: cs
  real(dp) :: res(size(ptj_veto))

  call set_process_and_parameters(cs, proc, M, muR, muF, Q, as, p,&
       &jet_radius, observable)

  if (small_r) then
     cs%small_r = .true.
     cs%small_r_R0 = small_r_R0
  else
     cs%small_r = .false.
  endif

  res = resummed_sigma_nolumi(ptj_veto, cs, order)
end function sudakov
