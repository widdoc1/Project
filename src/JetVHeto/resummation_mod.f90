module resummation_mod
  use types_mod; use consts_mod
  use rad_tools_mod
  use emsn_tools_mod
  implicit none

  private
  public :: resummed_sigma

contains
  !======================================================================
  function resummed_sigma(pt, cs, order) result(sigma)
    real(dp),                  intent(in) :: pt
    type(process_and_parameters), intent(in) :: cs
    integer,                   intent(in) :: order
    real(dp)  :: sigma
    !------------------------------
    real(dp) :: L_tilde, lambda
    real(dp) :: normalisation
    integer :: i 
    

    L_tilde = Ltilde(pt/cs%Q, cs%p)
    lambda  = get_lambda(L_tilde, cs)

    sigma = zero
    if (lambda < half)  sigma = exp(Rad(L_tilde, cs, order))

! this is redundant, tidy up
    if (order == order_LL) then
       sigma = sigma
    else if (order == order_NLL) then
       sigma = sigma
    else
!       if (order /= order_NNLL) call wae_error("expected order_NNLL, found", intval=order)  
             ! - Original code -
        sigma = sigma * (1 + &
             & non_incl(cs%jet_radius,'all') * Rad_p(lambda)*two*cs%as2pi/(1-two*lambda))
!       endif
    endif

  end function resummed_sigma

end module resummation_mod
