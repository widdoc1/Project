!========================================================
!--------------------------------------------------------
! Includes a module providing the resummation ingredients
! -------------------------------------------------------
!========================================================
module rad_tools_mod
  use types_mod; use consts_mod
  use qcd_mod
  implicit none

  private

  ! define the orders that we will support
  integer, parameter, public :: order_LL = 0
  integer, parameter, public :: order_NLL = 1
  integer, parameter, public :: order_NNLL = 2

  ! a type to conveniently hold the coupling and information on 
  ! scale ratios
  type process_and_parameters
     character(len=5) :: collider, proc
     real(dp) :: rts
     real(dp) :: alphas_muR, as2pi
     real(dp) :: M, muR, muF, Q, p
     real(dp) :: ln_muR2_M2
     real(dp) :: ln_Q2_M2
     real(dp) :: ln_Q2_muR2
     real(dp) :: ln_muF2_M2
     real(dp) :: M2_rts2
     real(dp) :: jet_radius
  end type process_and_parameters


  public :: order_string
  public :: process_and_parameters, set_process_and_parameters, get_lambda
  public  :: init_proc, Ltilde, Rad, Rad_p, Rad_pNNLL, Rad_s, g1, g2, g3
  real(dp), public :: A(3), B(2), H(2), Hm(1)  
  real(dp), public :: CC, BB, H1, Hm1 ! colour factor & our B

contains

  !======================================================================
  function order_string(order) result(res)
    integer,    intent(in) :: order
    character(len=4) :: res
    !------------------------------------
    select case(order)
    case(order_LL)
       res = "LL"
    case(order_NLL)
       res = "NLL"
    case(order_NNLL)
       res = "NNLL"
!    case default
!       call wae_error("value of order was not recognized in order_string",intval=order)
    end select
  end function order_string

  !======================================================================
  subroutine init_proc(cs)
    implicit none
    type(process_and_parameters), intent(in) :: cs
  !======================================================================
    real(dp) :: as_pow

    select case(trim(cs%proc))
       !! the hard part of the coefficient function (form factor) is not part of B(2)
       !! -----------> the form factor (H function) MUST be evaluated at the hard scale.
       !!
       !! H(1) and B(2) determined from Grazzini, De Florian
       !! hep-ph/0108273, with some additional manipulations
       !! reflecting that in our convention H(1) multiplies alphas(muR).
       !!
       !! A(3) from Becher, Neubert 1007.4005 [hep-ph]  
       !!
       !! The explicit formulae are in section 1 of the supplementary
       !! material of BMSZ.
    case('H')
       A(1) = two*ca_def
       CC   = ca_def
       A(2) = cmw_K*A(1) 
       A(3) = cmw_K2*A(1)+pi*beta0*ca_def*(ca_def*(808._dp/27._dp-28._dp*zeta3)-224._dp/27._dp*tf_def)
       BB   = -beta0 * pi/ca_def
       B(1) = -two*twopi*beta0     
       B(2) = -two*(ca_def**2*(8._dp/3._dp+three*zeta3)-cf_def*tf_def-four/three*ca_def*tf_def) &
            & +twopi_beta0*zeta2*ca_def !! Becher & Neubert arxiv:1205.3806v1 had additional: +8._dp*zeta3*ca_def**2
    case('DY')
       A(1) = two*cf_def
       CC   = cf_def
       A(2) = cmw_K*A(1)
       A(3) = cmw_K2*A(1)+pi*beta0*cf_def*(ca_def*(808._dp/27._dp-28._dp*zeta3)-224._dp/27._dp*tf_def) 
       BB   = -three/four
       B(1) = -three*cf_def
       B(2) =   -two*(cf_def**2*(-half*pisq+3._dp/8._dp+6._dp*zeta3)&
            & + cf_def*ca_def*(11._dp/18._dp*pisq+17._dp/24._dp-three*zeta3) &
            & + cf_def*tf_def*(-one/6._dp-two/9._dp*pisq)) &
            & + twopi_beta0*zeta2*cf_def 
   end select

  end subroutine init_proc

  !======================================================================
  ! 
  subroutine set_process_and_parameters(cs, collider, proc, rts, M, muR, muF, as, Q, p, &
       &jet_radius,loop_mass,observable)
    type(process_and_parameters), intent(out) :: cs
    character(len=*),          intent(in)  :: collider, proc, loop_mass, observable
    real(dp),                  intent(in)  :: rts, M, muR, muF, Q, p, as, jet_radius
    !-------------------------------------------------------
    cs%collider   = collider
    cs%proc       = proc
    cs%rts        = rts
    cs%M          = M
    cs%Q          = Q
    cs%muR        = muR
    cs%muF        = muF
    cs%alphas_muR = as
    cs%as2pi      = cs%alphas_muR/twopi
    cs%p          = p
    cs%ln_muR2_M2 = 2*log(muR/M)
    cs%ln_Q2_M2   = 2*log(Q/M)
    cs%ln_Q2_muR2 = 2*log(Q/muR)
    cs%ln_muF2_M2 = 2*log(muF/M)
    cs%M2_rts2    = M**2/rts**2
    cs%jet_radius = jet_radius
    ! maybe put a check here to make sure only pass 'H' or 'DY'
  end subroutine set_process_and_parameters

  !=========================================================
  ! Computes the modified logarithm
  function Ltilde(v,p) result(res)
    real(dp), intent(in) :: v(:), p
    real(dp) :: Z(size(v))
    real(dp) :: res(size(v))
    
    if (p >= 0) then
       ! Ltilde is defined setting ptmax->infty
       Z = one
       res = 1/p*log(1/v**p+1)*Z 
    else
       write(0,*) 'Warning: using unmodified logarithm'
       res = log(1/v)
    end if
    
  end function Ltilde

  !======================================================================
  function get_lambda(L, cs) result(res)
    real(dp),                  intent(in) :: L(:)
    type(process_and_parameters), intent(in) :: cs
    !-----------------------------------------
    real(dp):: res(size(L))

    res = cs%alphas_muR * L * beta0
  end function get_lambda
  
    
  !=============================================================================
  ! Compute the sudakov exponent (minus the radiator)
  ! for an array values of the logarithm L =
  ! ln(Q/ptveto), or (at the user's choice) a modified logarithm that
  ! is equivalent for small ptveto.
  function Rad(L, cs, order) result(res)
    real(dp),                  intent(in) :: L(:)
    type(process_and_parameters), intent(in) :: cs
    integer,                   intent(in) :: order !(0=LL,1=NLL,2=NNLL)
    real(dp)                              :: res(size(L))
    !----------------------------------
    real(dp):: lambda(size(L))

    lambda = get_lambda(L, cs)

    res = L*g1(lambda) 
    if (order >= order_NLL)  res = res + g2(lambda,cs)
    if (order >= order_NNLL) then
       res = res + cs%as2pi*two*g3(lambda,cs) 
    end if
    

!    if (order > order_NNLL .or. order < order_LL) call wae_error("Illegal value for order", intval=order)
  end function Rad

  !=========================================================
  ! The first derivative of the radiator @ NLL accuracy (sufficient
  ! for jet veto!)
  ! 
  ! It is the derivative of Lg1(as L) with respect to L.
  function Rad_p(lambda) result(res)
    real(dp), intent(in) :: lambda(:)
    !----------------------------------
    real(dp):: res(size(lambda))

    res = A(1)/pi/beta0*two*lambda/(one-two*lambda)
  end function Rad_p


  !=========================================================
  ! The NNLL contribution to the first derivative of the
  ! radiator necessary for ptB
  function Rad_pNNLL(lambda,cs) result(res)
    real(dp), intent(in) :: lambda(:)
    type(process_and_parameters), intent(in) :: cs
    !----------------------------------
    real(dp):: res(size(lambda))

    res =  twopi*cs%as2pi*((-two*lambda*pi*beta1*A(1)*log(one-two*lambda) + &
         &  beta0*(A(2)*lambda + beta0*(one - two*lambda)*pi*B(1) + &
         &  beta0*pi*A(1)*((one - two*lambda)*(-cs%ln_Q2_M2) + &
         &  two*lambda*(-cs%ln_Q2_muR2))))/(beta0**2*(one - two*lambda)**2*pi**2))
  end function Rad_pNNLL

  !=========================================================
  ! The second derivative of the radiator @ NLL accuracy 
  function Rad_s(lambda,cs) result(res)
    real(dp), intent(in) :: lambda(:)
    type(process_and_parameters), intent(in) :: cs
    !----------------------------------
    real(dp):: res(size(lambda))

    res = 4*A(1)*cs%as2pi/((one - two*lambda)**2)
  end function Rad_s


  !=========================================================
  ! Eq.10a from BMSZ
  function g1(lambda) result(res)
    real(dp), intent(in)  :: lambda(:)
    !----------------------------------
    real(dp) :: res(size(lambda))

    where (lambda /= 0)
       res = A(1)/pi/beta0*(one+half*log(one-two*lambda)/lambda)   
    else where
       res = zero
    end where
  end function g1

  !=========================================================
  ! Eq.10b from BMSZ
  function g2(lambda, cs) result(res)
    real(dp),                  intent(in) :: lambda(:)
    type(process_and_parameters), intent(in) :: cs
    !----------------------------------
    real(dp) :: res(size(lambda))

    where (lambda /= 0)
       res = B(1)/twopi/beta0*log(one-two*lambda)&
            &-A(2)/four/pisq/beta0**2*(two*lambda/(one-two*lambda)+log(one-two*lambda))&
            & +A(1)/twopi/beta0*((two*lambda/(one-two*lambda)&
            &+log(one-two*lambda))*cs%ln_Q2_muR2+log(one-two*lambda)*(-cs%ln_Q2_M2)) &
            & +A(1)/two*(pisq*beta1)/(pi*beta0)**3*(half*log(one-two*lambda)**2&
            &+(log(one-two*lambda)+two*lambda)/(one-two*lambda))
    else where
       res = zero
    end where

  end function g2

  !=========================================================  
  ! Eq.10c from BMSZ
  function g3(lambda, cs) result(res)
    real(dp),                  intent(in) :: lambda(:)
    type(process_and_parameters), intent(in) :: cs
    !----------------------------------
    real(dp) :: res(size(lambda))

    where (lambda /= 0)
       res = -half*A(3)/8._dp/pisq/beta0**2*(two*lambda/(one-two*lambda))**2 &
            & -(B(2)+A(2)*(-cs%ln_Q2_M2))/four/pi/beta0*two*lambda/(one-two*lambda) &
            & +A(2)/four*(pisq*beta1)/(pi*beta0)**3*(lambda*(three*two*lambda-two)/(one-two*lambda)**2 &
            & -(one-four*lambda)/(one-two*lambda)**2*log(one-two*lambda)) &
            & +(B(1)+A(1)*(-cs%ln_Q2_M2))/two*(pisq*beta1)/pisq/beta0**2&
            &        *(two*lambda/(one-two*lambda)+log(one-two*lambda)/(one-two*lambda)) &
            & -half*A(1)/two*(two*lambda/(one-two*lambda))**2*cs%ln_Q2_muR2**2 &
            & +((B(1)+A(1)*(-cs%ln_Q2_M2))/two*two*lambda/(one-two*lambda)&
            & +A(2)/four/pi/beta0*(two*lambda/(one-two*lambda))**2 &
            & +A(1)/two*pisq*beta1/pisq/beta0**2*(two*lambda/(one-two*lambda) &
            & +(one-four*lambda)/(one-two*lambda)**2*log(one-two*lambda)))*cs%ln_Q2_muR2 &
            & +A(1)/two*((pisq*beta1)**2/two/(pi*beta0)**4*(one-four*lambda)/&
            & (one-two*lambda)**2*log(one-two*lambda)**2 &
            & +log(one-two*lambda)*((pi*beta0*pi**3*beta2-(pisq*beta1)**2)/&
            &          (pi*beta0)**4+(pisq*beta1)**2/(pi*beta0)**4/(one-two*lambda))&
            & +one/(pi*beta0)**4*lambda/(one-two*lambda)**2&
            &     *(pi*beta0*pi**3*beta2*(two-three*two*lambda)+(pisq*beta1)**2*two*lambda))
    else where
       res = zero
    end where
  end function g3

  !========================================================

end module rad_tools_mod
