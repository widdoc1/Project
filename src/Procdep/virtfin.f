      subroutine virtfin(p, msq, msqv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'b0.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'jetvheto.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'qcdcouple.f'
      integer  :: j,k
      real(dp) :: p(mxpart, 4)
      real(dp) :: msq(-nf:nf,-nf:nf), msqv(-nf:nf,-nf:nf)
      real(dp) :: dot, virt, xl12, T2, ga, I

      xl12=log(two*dot(p,1,2)/musq)

!     set up insertion operator, the
!     universal singular structure at the 1-loop level
      select case(trim(born_config))
      case('H')
         T2 = ca
         ga = b0
      case('DY')
         T2 = cf
         ga = three/two * cf
      end select
      I = two*T2*(epinv*epinv2-epinv*xl12+half*xl12**2)
     &     +two*ga*(epinv-xl12)

      do j=-nf,nf
         do k=-nf,nf
!           catch matrix elements that are zero at the born
!           level, these do not have virtual corrections to
!           modify.
            if (msq(j,k) == 0._dp) then
               msqv(j,k) = 0._dp
            else
!              subtract divergences with insertion operator
               msqv(j,k) = msqv(j,k) + ason2pi*msq(j,k)*I
!              additional C*pi**2/6 due to coupling mismatch
     &              + ason2pi*msq(j,k)*T2*pisqo6
!              change into form for jet veto resummation
     &              - ason2pi*msq(j,k)*(-two*ga
     &              + T2*log(two*dot(p,1,2)/q_scale**2))
     &              * log(two*dot(p,1,2)/q_scale**2)
            endif
         enddo
      enddo
      end subroutine virtfin
