!========================================================
!--------------------------------------------------------
! Interface to access the Sudakov form factor
! -------------------------------------------------------
!========================================================

!  real(kind(1d0))
function sudakov(proc, M, muR, muF, as, Q, p, jet_radius,&
     &observable, ptj_veto, order) result(res)
  use types
  use consts
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
  !  sudakov = resummed_sigma(ptj_veto, cs, order)
    res = resummed_sigma(ptj_veto(1), cs, order)
end function sudakov
