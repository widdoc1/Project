module resummation
  use types; use consts_dp
  use rad_tools
  use warnings_and_errors
  use emsn_tools
  implicit none

  private
  public :: resummed_sigma_nolumi

contains
  !======================================================================
  function resummed_sigma_nolumi(pt, cs, order) result(sigma)
    real(dp),                  intent(in) :: pt(:)
    type(process_and_parameters), intent(in) :: cs
    integer,                   intent(in) :: order
    real(dp)  :: sigma(size(pt))
    !------------------------------
    real(dp) :: L_tilde(size(pt)), lambda(size(pt))
    real(dp) :: av_lnz(size(pt)), av_ln2z(size(pt)), as2pi_pt(size(pt)), non_incl_largeR(size(pt))
    integer :: i 
    real(dp) :: tmp(size(pt))
    

    if (cs%observable /= 'ptj') then
       call wae_error("observable not implemented")
    end if
    
    L_tilde = Ltilde(pt/cs%Q, cs%p)
    lambda  = get_lambda(L_tilde, cs)

    sigma = zero
    where (lambda < half)  sigma = exp(Rad(L_tilde, cs, order))

    select case(order)
    case(order_LL, order_NLL)
       ! do nothing
    case(order_NNLL)
          
       if (cs%small_r) then
          ! - Added by FD -
          ! get the all-order result
          as2pi_pt = cs%as2pi/(one-two*lambda)
          ! for ln R resummation, calculate b0 and t
          ! b0 = (11.*ca_def-2.*nf_def)/6.
          ! t = log(1/(1-as2pi_pt*b0*log(1/(cs%jet_radius**2))/(2.*pi)))/b0
          av_lnz = av_lnz_smallR(as2pi_pt,cs%jet_radius,cs%small_r_R0)
          ! subtract lnR term at O(as) from non-inclusive correction
          non_incl_largeR = as2pi_pt*(non_incl(cs%jet_radius,'all') - &
               & non_incl_lnR(cs%jet_radius,cs%small_r_R0))
          ! include the non-inclusive correction
          ! this should be checked carefully
          sigma = sigma * (exp(-Rad_p(lambda)*av_lnz) + &
               & Rad_p(lambda)*two*non_incl_largeR)
          ! sigma = sigma * (1 - Rad_p(lambda)*av_lnz + &
          !      & Rad_p(lambda)*two*non_incl_largeR)
          
          if (cs%small_r_ln2z) then
             ! extra subleading terms: MD & GPS temporary investigations (2015-02-16)
             av_ln2z = av_ln2z_smallR(as2pi_pt,cs%jet_radius,cs%small_r_R0)
             !
             ! fix a normalisation so as not to change results at large pt
             tmp(1:1) = exp(-A(1)*cs%as2pi*two &
                  &          * av_ln2z_smallR((/cs%as2pi/),cs%jet_radius,cs%small_r_R0))
             ! put it together
             sigma = sigma * exp(-A(1)*as2pi_pt*two * av_ln2z) / tmp(1)
             
          end if
       else
          ! - Original code -
          sigma = sigma * (1 + &
               & non_incl(cs%jet_radius,'all') * Rad_p(lambda)*two*cs%as2pi/(1-two*lambda))
       end if
    case default
       call wae_error("unexpected order", intval=order)
    end select

  end function resummed_sigma_nolumi


end module resummation
