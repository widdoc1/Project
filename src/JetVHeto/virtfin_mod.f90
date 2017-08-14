!========================================================
!--------------------------------------------------------
! Includes a module providing the 
! -------------------------------------------------------
!========================================================
module virtfin_mod
  use types_mod; use consts_mod
  use qcd_mod
  use rad_tools_mod
  implicit none

  private
  ! real(dp), parameter :: Tq2 = 

  public :: virtfin

contains

  subroutine virtfin(p, msq, msqv, cs)
    implicit none
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'scale.f'
    include 'facscale.f'
    include 'epinv.f'
    include 'epinv2.f'
    include 'qcdcouple.f'
    include 'JetVHeto.f'
    include 'JetVHeto_opts.f'
    integer :: j,k
    real(dp), intent(in) :: p(mxpart, 4)
    real(dp), intent(inout) :: msq(-nf:nf,-nf:nf), msqv(-nf:nf,-nf:nf)
    real(dp) :: dot, virt, xl12, I
    type(process_and_parameters), intent(in) :: cs
    
    
    xl12=log(two*dot(p,1,2)/musq)

    ! set up insertion operator
    I = rad_A(1)*(epinv*epinv2-epinv*xl12+half*xl12**2) &
            +(-rad_B(1))*(epinv-xl12)

    do j=-nf,nf
       do k=-nf,nf
          ! subtract divergences with insertion operator
          msqv(j,k) = msqv(j,k) + ason2pi*msq(j,k)*I

          ! additional C*pi**2/6 due to coupling mismatch
          msqv(j,k) = msqv(j,k) + ason2pi*msq(j,k)*(half*rad_A(1))*pisqo6

          ! change into form for resummation
          msqv(j,k) = msqv(j,k)+ason2pi*msq(j,k)*(-half*rad_A(1) & 
               *cs%ln_Q2_M2+rad_B(1))*cs%ln_Q2_M2 &
               +two*as_pow*pi*beta0*cs%ln_muR2_M2
       enddo
    enddo


    ! loop over all incoming partons
    ! do j=nf,nf
    !    do k=-nf,nf

    !       if (abs(msq(j,k)) < 1E-9_dp) then
    !       else
    !          xl12=log(two*dot(p,1,2)/musq)
    !          msqv(j,k)/msq(j,k)/ason2pi/two ! get the factorised divergent piece

    !       ! subtract away universal divergent structure
    !       msqv(j,k)=msqv(j,k)-(-(epinv*epinv2-epinv*xl12+half*xl12**2) &
    !            - three/two*(epinv-xl12)) ! divergent parts

    !       ! additional pi**2/6
    !       msqv(j,k)=msqv(j,k)+pisqo6*(half*rad_A(1)) ! fix casimir

    !       ! change into form for resummation
    !       msqv(j,k)=msqv(j,k)+(-half*rad_A(1) &
    !       &        *cs%ln_Q2_M2+rad_B(1))*cs%ln_Q2_M2 &
    !       &        + two*as_pow*pi*beta0*cs%ln_muR2_M2

    !       ! restore prefactors
    !       msqv(j,k)=ason2pi*two*msq(j,k)*msqv(j,k)

    !     endif
    !   enddo
    ! enddo


  end subroutine 

end module virtfin_mod
