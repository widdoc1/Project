!========================================================
!--------------------------------------------------------
! Includes a module providing the 
! -------------------------------------------------------
!========================================================
module virtfin_mod
  use types_mod; use consts_mod
  use qcd_mod
  implicit none

  private

contains

  subroutine virtfin(p, msq, msqv, resm_opts)
    implicit none
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'epinv.f'
    include 'epinv2.f'
    include 'JetVHeto.f'
    include 'JetVHeto_opts.f'
    integer :: j,k
    real(dp), intent(in) :: p(mxpart, 4)
    real(dp), intent(in) :: msq(-nf:nf,-nf:nf), msqv(-nf:nf,-nf:nf)
    real(dp) :: dot, virt, xl12
    type(process_and_parameters), intent(in) :: resm_opts
    
    ! loop over all incoming partons
    do j=nf,nf
      do k=-nf,nf

        if (abs(msq(j,k)) < 1E-9_dp) then
        else
          xl12=log(two*dot(p,1,2)/musq)
          msqv(j,k)/msq(j,k)/ason2pi/two ! get the factorised divergent piece

          ! subtract away universal divergent structure
          msqv(j,k)=msqv(j,k)-() ! divergent parts

          ! additional pi**2/6
          msqv(j,k)=msqv(j,k)+pisqo6*C ! fix casimir

          ! change into form for resummation
          msqv(j,k)=msqv(j,k)+(-half*rad_A(1)
          &        *resm_opts%ln_Q2_M2+rad_B(1))*resm_opts%ln_Q2_M2
          &        + two*as_pow*pi*beta0*resm_opts%ln_muR2_M2

          ! restore prefactors
          msqv(j,k)=ason2pi*two*msq(j,k)*msqv(j,k)

        endif
      enddo
    enddo


  end subroutine 
