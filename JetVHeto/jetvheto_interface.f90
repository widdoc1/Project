!========================================================
!--------------------------------------------------------
! Interface to access the Sudakov form factor
! -------------------------------------------------------
!========================================================

!  real(kind(1d0))
function sudakov(proc, M, muR, muF, as, Q, p, jet_radius,&
     &observable, ptj_veto, order) result(res)
  use types;  use consts_dp
  use rad_tools
  use resummation
  implicit none
  character(len=*),            intent(in)  :: proc, observable
  integer,                     intent(in)  :: order
  real(dp),                    intent(in)  :: M, muR, muF, Q, p, as,&
                                             &jet_radius, ptj_veto(:)
  type(process_and_parameters)             :: cs
  real(dp) :: res(size(ptj_veto))
  
  call set_process_and_parameters(cs, proc, M, muR, muF, as, Q, p, jet_radius,&
       &observable)
    res = resummed_sigma_nolumi(ptj_veto, cs, order)
end function sudakov
