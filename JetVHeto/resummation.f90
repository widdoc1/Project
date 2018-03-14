module resummation
  use types; use consts
  use rad_tools
  use emsn_tools
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

    L_tilde = Ltilde(pt/cs%Q, cs%p)
    lambda  = get_lambda(L_tilde, cs)


    sigma = zero
    if (lambda < half)  sigma = exp(Rad(L_tilde, cs, order))
    if (order == order_NNLL) then
       sigma = sigma * (1 + &
            & non_incl(cs%jet_radius,'all') * Rad_p(lambda)*two*cs%as2pi/(1-two*lambda))
    endif

  end function resummed_sigma

end module resummation
