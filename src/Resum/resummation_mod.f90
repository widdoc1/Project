module resummation_mod
  use types_mod; use consts_dp
  use rad_tools
  use pdfs_tools
  use warnings_and_errors
  use emsn_tools
  use ew_parameters; use mass_corr
  use special_functions_mod
  implicit none

  private
  public :: resummed_sigma

contains
  !======================================================================
  function resummed_sigma(pt, cs, order, dlumi_lumi) result(sigma)
    real(dp),                  intent(in) :: pt(:)
    type(process_and_parameters), intent(in) :: cs
    integer,                   intent(in) :: order
    real(dp)  :: sigma(size(pt))
    !------------------------------
    real(dp) :: L_tilde(size(pt)), lambda(size(pt)), lNLL(size(pt)), &
         &lNNLL(size(pt)), dlNLL(size(pt))
    real(dp) :: normalisation
    real(dp) :: rp(size(pt)),rs(size(pt)),resum_fact(size(pt)),drp(size(pt))
    real(dp) :: av_lnz(size(pt)), av_ln2z(size(pt)), as2pi_pt(size(pt)), non_incl_largeR(size(pt))
    integer :: i 
    real(dp) :: tmp(size(pt))
    

    L_tilde = Ltilde(pt/cs%Q, cs%p)
    lambda  = get_lambda(L_tilde, cs)

    sigma = zero
    where (lambda < half)  sigma = exp(Rad(L_tilde, cs, order))

    if (present(dlumi_lumi))  dlumi_lumi = zero    

    if (order == order_LL) then
       sigma = sigma
    else if (order == order_NLL) then
       sigma = sigma
    else
       if (order /= order_NNLL) call wae_error("expected order_NNLL, found", intval=order)  
             ! - Original code -
        sigma = sigma * (1 + &
             & non_incl(cs%jet_radius,'all') * Rad_p(lambda)*two*cs%as2pi/(1-two*lambda))
       endif
    endif

  end function resummed_sigma


  !! extend dgamma to negative arguments
  function gammafull(x) result (res)
    use types; use consts_dp
    use special_functions
    real(dp)             :: res
    real(dp), intent(in) :: x
    if (x < zero) then
       res = pi/dgamma(one-x)/sin(pi*x)
    else
       res = dgamma(x)
    end if
  end function gammafull


end module resummation_mod
