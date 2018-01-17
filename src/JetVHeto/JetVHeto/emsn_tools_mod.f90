!========================================================
!--------------------------------------------------------
! Includes a module for the non-inclusive correction
! -------------------------------------------------------
!========================================================
module emsn_tools_mod
  use types_mod; use consts_mod
  use rad_tools_mod; use qcd_mod
  implicit none
  private

  public :: non_incl

  real(dp), parameter :: lntwo =&
       & 0.693147180559945309417232121458176568075_dp

contains

  ! contains the analytic expressions for the non-inclusive (R-dependent)
  ! corrections
  function non_incl(R,part) result(res)
    real(dp), intent(in) :: R
    real(dp) :: correl, indep, res
    character(len=*) :: part

    ! Taylor series of the correlated emission contribution neglecting terms of O(R**8)
    ! Good agreement with MC below R=2.6
    !
    ! Taken from 1203.5773 A.16 (and independently derived by PFM)
    correl = -(-131._dp+12._dp*pisq+132._dp*lntwo)/72._dp*ca_def*log(R/1.74_dp) &
         & -(23.0_dp-24.0_dp*lntwo)/72._dp*nf_def*log(R/0.84_dp) &
         & +((1429._dp+3600._dp*pisq+12480._dp*lntwo)*ca_def+&
         & (3071._dp-1680._dp*lntwo)*nf_def)/172800._dp*R**2           &
         & +((-9383279._dp-117600._dp*pisq+1972320._dp*lntwo)*ca_def+&
         & two*(178080._dp*lntwo-168401._dp)*nf_def)/406425600._dp*R**4 &
         &  +((74801417._dp-33384960._dp*lntwo)*ca_def+&
         & (7001023._dp-5322240._dp*lntwo)*nf_def)/97542144000._dp*R**6

    if (R>pi) then
       ! Taken from 1206.4998 Eq.(45)
       indep = -half*A(1)*((pi/6._dp*R**2-1/pi/8._dp*R**4)*&
            & atan(pi/sqrt(R**2-pisq)) &!this avoids Floating point exceptions for R<pi
            & +(R**2/8._dp-pisq/12._dp)*sqrt(R**2-pisq))
    else
       ! Taken from 1203.5773 A.18 (and independently derived by PFM)
       indep = -half*A(1)*(pisq*R**2/12._dp-R**4/16._dp)
    end if

    select case(trim(part))
    case('correl')
       res = correl
    case('indep')
       res = indep
    case('all')
       res = correl+indep
    end select
  end function non_incl

end module emsn_tools_mod
