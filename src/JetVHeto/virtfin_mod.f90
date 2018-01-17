!========================================================
!--------------------------------------------------------
! Includes a module to correct the finite part of the
! virtual corrections to our conventions
! -------------------------------------------------------
!========================================================
module virtfin_mod
  use types_mod; use consts_mod
  use qcd_mod
  use rad_tools_mod
  implicit none
  private

  public :: virtfin

contains

  subroutine virtfin(p, msq, msqv)
    implicit none
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'scale.f'
    include 'facscale.f'
    include 'jetvheto.f'
    include 'epinv.f'
    include 'epinv2.f'
    include 'qcdcouple.f'
    integer :: j,k
    real(dp), intent(in) :: p(mxpart, 4)
    real(dp), intent(inout) :: msq(-nf:nf,-nf:nf), msqv(-nf:nf,-nf:nf)
    real(dp) :: dot, virt, xl12, I

    xl12=log(two*dot(p,1,2)/musq)

    ! set up insertion operator, the
    ! universal singular structure at the 1-loop level
    I = A(1)*(epinv*epinv2-epinv*xl12+half*xl12**2) &
            -B(1)*(epinv-xl12)

    do j=-nf,nf
       do k=-nf,nf
          ! catch matrix elements that are zero at the born
          ! level, these do not have virtual corrections to
          ! modify.
          if (msq(j,k) == 0._dp) then
             msqv(j,k) = 0._dp
          else
          ! subtract divergences with insertion operator
          msqv(j,k) = msqv(j,k) + ason2pi*msq(j,k)*I

          ! additional C*pi**2/6 due to coupling mismatch
          msqv(j,k) = msqv(j,k) + ason2pi*msq(j,k)*(half*A(1))*pisqo6

          ! change into form for resummation
          msqv(j,k) = msqv(j,k) - ason2pi*msq(j,k)*(B(1) &
               + half*A(1)*log(two*dot(p,1,2)/q_scale**2))&
               * log(two*dot(p,1,2)/q_scale**2)
          endif
       enddo
    enddo

  end subroutine virtfin

end module virtfin_mod
