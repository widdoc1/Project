!========================================================
!--------------------------------------------------------
! Includes a module to correct the finite part of the
! virtual corrections to our conventions
! -------------------------------------------------------
!========================================================
module jetvheto_interface
  use types_mod
  use consts_mod
  use rad_tools_mod
  use resummation_mod
  implicit none
  private

  public :: sudakov

contains

  function sudakov(p) result(res)
    implicit none
    include 'mxpart.f'
    include 'kpart.f'
    include 'born_config.f'
    include 'scale.f'
    include 'facscale.f'
    include 'jetvheto.f'
    include 'qcdcouple.f'
    real(dp) :: Rcut
    common/Rcut/Rcut

    real(dp), intent(in) :: p(mxpart,4)
    real(dp) :: res
    integer  :: order
    real(dp) :: dot, M_B
    type(process_and_parameters) :: cs

    M_B = sqrt(two*dot(p,1,2))
    call set_process_and_parameters(cs, trim(born_config), M_B, scale, &
         facscale, as, q_scale, p_pow, Rcut, observable)

    select case(kpart)
    case(kll)
       order = order_LL
    case(knll)
       order = order_NLL
    case(knnll)
       order = order_NNLL
    end select

    res = resummed_sigma(ptj_veto, cs ,order)

  end function sudakov

end module
