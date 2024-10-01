! ! common blocks to think about including
! !
! ! - plabel: particle labels
! ! - ptilde: particle (ptilde) and jet (ptildejet) momenta
! ! - npart: number of particles
! ! - constants - things like pi, mxpart, etc.
! ! - nplot.f
! subroutine usernplotter()
!
! end subroutine usernplotter

!----------------------------------------------------------------------
! Warning: if the mcfm result here is false, then overriding it may
! cause inconsistent results elsewhere (e.g. scale setting, nplotter),
! because the MCFM jets will not necessarily have been properly
! defined and stored.
!
! (this could be worked around - but is it wise to?)
function userincludedipole(nd, ppart, mcfm_result)
  implicit none
  include 'types.f'
  include 'constants.f'
  include 'npart.f'
  include 'mxpart.f'
  include 'jetlabel.f'
  include 'nproc.f'
  include 'jetvheto.f'
  integer,  intent(in) :: nd
  real(dp), intent(in) :: ppart(mxpart,4)
  logical,  intent(in) :: mcfm_result
  logical :: userincludedipole
   !------------
  logical :: bin,makecuts,ATLAS_hww2017
  common/bin/bin
  common/makecuts/makecuts

  real(dp) :: ptj

  ! take the MCFM result as the default choice
  userincludedipole = mcfm_result
  ! return

  if (makecuts) then
     if ( (nproc == 61) .or. (nproc == 66) .or. (nproc == 69) &
          & .or. (nproc == 113) .or. (nproc == 123) .or. (nproc == 124) & 
          & .or. (nproc == 125) .or. (nproc == 126) ) then
        if (ATLAS_hww2017(ppart)) then
           userincludedipole = mcfm_result
           return
        else
           userincludedipole = .false.
           return
        endif
     endif

     ! DY validation
     if (nproc==31) then
        ptj = sqrt(ppart(5,1)**2+ppart(5,2)**2)
        if (ptj > ptj_veto) then
           userincludedipole = .false.
           return
        endif
     endif

  else
     userincludedipole = mcfm_result
     return
  endif

  return
end function userincludedipole

!----------------------------------------------------------------------
! Variables passed to this routine:
!
!      p:  4-momenta of jets in the format p(i,4)
!          with the particles numbered according to the input file
!          and components labelled by (px,py,pz,E).
!
!     wt:  weight of this event
!
!    wt2:  weight^2 of this event
!
!     nd:  an integer specifying the dipole number of this contribution
!          (if applicable), otherwise equal to zero

subroutine userplotter(pjet, wt, wt2, nd)
  use rad_tools
  implicit none
  include 'types.f'
  include 'constants.f'
  include 'nf.f'
  include 'mxpart.f'
  include 'jetlabel.f'
  include 'ptilde.f'
  include 'npart.f'
  include 'nplot.f'

  ! needed for the resummation
  include 'kpart.f'
  include 'nproc.f'
  include 'kprocess.f'
  include 'noglue.f'
  include 'scale.f'
  include 'facscale.f'
  include 'jetvheto.f'
  include 'qcdcouple.f'
  real(dp) :: Rcut
  common/Rcut/Rcut
  real(dp) :: dot, M_B
  integer  :: order
  real(dp) :: sudakov_arr(1)
  interface
     function sudakov(proc, M, muR, muF, Q, as, p, jet_radius,&
          &observable, small_r, small_r_R0, ptj_veto, order) result(res)
       use types; use consts_dp
       use rad_tools
       use resummation
       implicit none
       character(len=*),            intent(in)  :: proc, observable
       integer,                     intent(in)  :: order
       real(dp),                    intent(in)  :: M, muR, muF, Q, as, p,&
            &small_r_R0, jet_radius, ptj_veto(:)
       logical :: small_r
       type(process_and_parameters)             :: cs
       real(dp) :: res(size(ptj_veto))
     end function sudakov
  end interface
  real(dp), intent(in) :: pjet(mxpart,4)
  real(dp), intent(inout) :: wt,wt2
  integer,  intent(in) :: nd
  integer, parameter :: tagbook=1, tagplot=2
  !---------------------------------------------
  integer :: i, iplot
  real(dp) :: pt, pttwo
  real(dp) :: pt3, pt4, pt5, pt6, pt7
  real(dp) :: m34, m36, m45, m56, m3456
  real(dp) :: MT1, MT2, MT3
  real(dp) :: ptll, pt45(4), ptmiss, pt36(4), MTll, MTmiss
  real(dp) :: r2, delphi
  real(dp) :: dphi
  real(dp) :: ptrel
  logical, save :: first = .true.
  integer :: tag
  
  if (first) then
     tag   = tagbook

     ! setup resummation parameters if applicable
     if (jetvheto) then
        if (      (kcase==kW_only) .or. (kcase==kZ_only)&
             .or. (kcase==kWWnpol) .or. (kcase==kWZbbar)&
             .or. (kcase==kWHbbar) .or. (kcase==kWHgaga)&
             .or. (kcase==kWH__WW) .or. (kcase==kWH__ZZ)&
             .or. (kcase==kZHbbar) .or. (kcase==kZHgaga)&
             .or. (kcase==kZH__WW) .or. (kcase==kZH__ZZ) ) then
           born_config='DY'
        elseif (  (kcase==kggfus0)&
             .or. (kcase==kHWW_4l) .or. (kcase==kHWW_tb)&
             .or. (kcase==kHWWint) .or. (kcase==kHWWHpi)&
             .or. (kcase==kggWW4l) .or. (kcase==kggWWbx)&
             .or. (kcase==kHWW2lq)&
             .or. (kcase==kHZZ_4l) .or. (kcase==kHZZ_tb)&
             .or. (kcase==kHZZint) .or. (kcase==kHZZHpi)&
             .or. (kcase==kggZZ4l) .or. (kcase==kggZZbx)&
             .or. (kcase==kHi_Zga)&
             .or. (kcase==kHVV_tb) .or. (kcase==kHVVint)&
             .or. (kcase==kHVVHpi) .or. (kcase==kggVV4l)&
             .or. (kcase==kggVVbx))  then
           born_config='H'
        elseif (  (kcase==kWWqqbr) .or. (kcase==kZZlept)) then
           if     (omitgg .eqv. .true.) then
              born_config='DY'
           elseif (ggonly .eqv. .true.) then
              born_config='H'
           else
              write(6,*) 'process=', nproc, 'contains contributions from qqb and' 
              write(6,*) 'gg initial states. To perform jet veto resummation'
              write(6,*) 'set omitgg=.true. in the input card to resum the qqb'
              write(6,*) 'contribution or ggonly=.true. to resum the gg contribution.'
              write(6,*)
              stop
           endif
        else
           write(6,*) 'process=', nproc, 'is not valid'
           write(6,*) 'for jet veto resummation'
           write(6,*)
           stop
        endif
        call init_proc(born_config)
     endif

     first = .false.
  else
     tag = tagplot
  endif
  iplot = nextnplot

  ! Sudakov
  if (do_suda) then

     M_B = sqrt(two*dot(pjet,1,2))

     select case(kpart)
     case(kll)
        order = 0
     case(knll)
        order = 1
     case(knnll)
        order = 2
     end select

     sudakov_arr=sudakov(born_config, M_B, scale, facscale, q_scale, as, p_pow,&
          &Rcut, observable, small_r, r_scale, (/ptj_veto/), order) 
     wt = wt*sudakov_arr(1)
     wt2 = wt**2
  end if

  !define quantities to plot
  pt3=pt(3,pjet)
  pt4=pt(4,pjet) 
  pt5=pt(5,pjet)
  pt6=pt(6,pjet) 
  pt7=pt(7,pjet)

  pt36 = pjet(3,:) + pjet(6,:)
  pt45 = pjet(4,:) + pjet(5,:)

  ptll = pttwo(4,5,pjet)
  ptmiss = pttwo(3,6,pjet)

  r2= (pjet(4,1)*pjet(5,1)+pjet(4,2)*pjet(5,2)) &
       /sqrt((pjet(4,1)**2+pjet(4,2)**2)*(pjet(5,1)**2+pjet(5,2)**2))
  if (r2 > +0.9999999_dp) r2=+1._dp
  if (r2 < -0.9999999_dp) r2=-1._dp
  delphi=acos(r2)

  ! define ptrel
  pt36 = pjet(3, :) + pjet(6, :)
  ptrel = ptmiss
  dphi = pi

  ! lepton 1
  dphi = min(dphi, &
       acos((pjet(4,1)*pt36(1) + pjet(4,2)*pt36(2))/pt4/ptmiss))
  ! lepton 2
  dphi = min(dphi, &
       acos((pjet(5,1)*pt36(1) + pjet(5,2)*pt36(2))/pt5/ptmiss))
  ! jets
  ! if (jets > 0) then
  !    dphi = min(dphi, &
  !         acos((pjet(7,1)*pt36(1) + pjet(7,2)*pt36(2))/pt(7,pjet)/ptmiss))
  ! endif

  if (dphi < pi/two) then
     ptrel = ptrel * sin(dphi)
  endif

  m45 = (pjet(4,4) + pjet(5,4))**2
  do i = 1, 3
     m45 = m45 - (pjet(4,i) + pjet(5,i))**2
  enddo
  m45 = sqrt(m45)

  m36 = (pjet(3,4) + pjet(6,4))**2
  do i = 1, 3
     m36 = m36 - (pjet(3,i) + pjet(6,i))**2
  enddo
  m36 = sqrt(m36)

  m34 = (pjet(3,4) + pjet(4,4))**2
  m3456 = (pjet(3,4) + pjet(4,4) + pjet(5,4) + pjet(6,4))**2
  do i = 1, 3
     m34 = m34 - (pjet(3,i) + pjet(4,i))**2
     m3456 = m3456 - (pjet(3,i) + pjet(4,i) + pjet(5,i) + pjet(6,i))**2
  enddo
  m34 = sqrt(m34)
  m3456 = sqrt(m3456)
  m56 = (pjet(5,4) + pjet(6,4))**2
  do i = 1,3
     m56 = m56 - (pjet(5,i) + pjet(6,i))**2
  enddo
  m56 = sqrt(m56)

  ! transverse mass proxy for MWW
  ! MT1
  MTll = sqrt(ptll**2 + m45**2)
  MT1 = (MTll + ptmiss)**2
  do i = 1, 2
     MT1 = MT1 - (pt45(i) + pt36(i))**2
  enddo
  MT1 = sqrt(MT1)

  ! MT2
  MT2 = sqrt(two*(ptll*ptmiss - (pt36(1)*pt45(1)+pt36(2)*pt45(2))))

  ! MT3
  MTmiss = sqrt(ptmiss**2 + MTll**2)
  MT3 = (MTll + MTmiss)**2
  do i = 1, 2
     MT3 = MT3 - (pt45(i) + pt36(i))**2
  enddo
  MT3 = sqrt(MT3)

  !!!! total cross section with the veto
  call bookplot(iplot,tag,'xsec',0.5_dp,wt,wt2,zip,one,one,'lin')
  iplot = iplot + 1

  !!!! plots for the paper

  ! m(3456) dists
  call bookplot(iplot,tag,'Mvv-1',m3456,wt,wt2,zip,2000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'Mvv-2',m3456,wt,wt2,2000._dp,4000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'Mvv-rest',m3456,wt,wt2,4000._dp,14000._dp,200._dp,'log')
  iplot = iplot + 1

  ! MT1
  call bookplot(iplot,tag,'MT1-1',MT1,wt,wt2,zip,2000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'MT1-2',MT1,wt,wt2,2000._dp,4000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'MT1-rest',MT1,wt,wt2,4000._dp,14000._dp,200._dp,'log')
  iplot = iplot + 1

  ! MT2
  call bookplot(iplot,tag,'MT2-1',MT2,wt,wt2,zip,2000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'MT2-2',MT2,wt,wt2,2000._dp,4000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'MT2-rest',MT2,wt,wt2,4000._dp,14000._dp,200._dp,'log')
  iplot = iplot + 1

  ! MT3
  call bookplot(iplot,tag,'MT3-1',MT3,wt,wt2,zip,2000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'MT3-2',MT3,wt,wt2,2000._dp,4000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'MT3-rest',MT3,wt,wt2,4000._dp,14000._dp,200._dp,'log')
  iplot = iplot + 1

  !!!! misc plots mostly for validation

    !!! transverse components of visible and invisible particle pairs

  ! m(45), m_ll
  call bookplot(iplot,tag,'mll',m45,wt,wt2,zip,1000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'mll_full',m45,wt,wt2,zip,14000._dp,200._dp,'log')
  iplot = iplot + 1

  ! pt(45), pt_ll
  call bookplot(iplot,tag,'ptll',ptll,wt,wt2,zip,1000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'ptll_full',ptll,wt,wt2,zip,14000._dp,200._dp,'log')
  iplot = iplot + 1

  ! m(36), m_nunu
  call bookplot(iplot,tag,'m_nunu',m36,wt,wt2,zip,1000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'m_nunu_full',m36,wt,wt2,zip,14000._dp,200._dp,'log')
  iplot = iplot + 1

  ! pt(36), pt_miss
  call bookplot(iplot,tag,'ptmiss',ptmiss,wt,wt2,zip,1000._dp,20._dp,'log')
  iplot = iplot + 1
  call bookplot(iplot,tag,'ptmiss_full',ptmiss,wt,wt2,zip,14000._dp,200._dp,'log')
  iplot = iplot + 1

  ! dphi(l+,l-)
  call bookplot(iplot,tag,'delphi',delphi,wt,wt2,zip,3.14_dp,0.1_dp,'lin')
  iplot = iplot + 1

    !!! pt of leading and sub-leading particles

  ! pt leading lepton
  call bookplot(iplot,tag,'pt(leading lept)',max(pt4,pt5),wt,wt2,zip,500._dp,5._dp,'log')
  iplot = iplot + 1

  ! pt subleading lepton
  call bookplot(iplot,tag,'pt(sub leading lept)',min(pt4,pt5),wt,wt2,zip,500._dp,5._dp,'log')
  iplot = iplot + 1

  ! pt leading neutrino
  call bookplot(iplot,tag,'pt(leading nu)',max(pt3,pt6),wt,wt2,zip,500._dp,5._dp,'log')
  iplot = iplot + 1

  ! pt sub leading neutrino
  call bookplot(iplot,tag,'pt(sub leading nu)',min(pt3,pt6),wt,wt2,zip,500._dp,5._dp,'log')
  iplot = iplot + 1

  ! ptmiss (pt neutrino pair)
  call bookplot(iplot,tag,'ptmiss',ptmiss,wt,wt2,zip,500._dp,5._dp,'log')
  iplot = iplot + 1

  ! ptmissrel
  call bookplot(iplot,tag,'ptmissrel',ptrel,wt,wt2,zip,500._dp,5._dp,'log')
  iplot = iplot + 1

    !!! special plots in the case of an additional jet from the real

  ! Njets
  call bookplot(iplot,tag,'N(jets)',real(jets,dp),wt,wt2,zip,5._dp,1._dp,'log')
  iplot = iplot + 1

  ! pt jets
  call bookplot(iplot,tag,'pt(leading jet)',pt7,wt,wt2,zip,500._dp,5._dp,'log')
  iplot = iplot + 1

  ! update nextnplot so we get userplots and generic plots from nplotter routines
  nextnplot = iplot

end subroutine userplotter

!----------------------------------------------------------------------
! user code to write info
subroutine userwriteinfo(unitno, comment_string, xsec, xsec_err, itno)
  implicit none
  include 'types.f'
  integer,          intent(in) :: unitno
  character*2,      intent(in) :: comment_string
  real(dp),         intent(in) :: xsec, xsec_err
  integer,          intent(in) :: itno

  !write(6,*) "have reached iteration number", itno
  !write(unitno,"(a,a)") comment_string, "any additional user comments"
  !call mcfmfwrite(unitno, comment_string//"any additional user comments")
end subroutine userwriteinfo

subroutine userhistofin(xsec,xsec_err,itno,itmx)
!	This function allows for extra user-defined operations
!	at the end of each iteration (itno>0) and at the end of
!	the run of the program (itno=0).
  implicit none
  include 'types.f'
  integer,  intent(in) :: itno,itmx
  real(dp), intent(in) :: xsec,xsec_err
end subroutine userhistofin

! subroutine userscale(event_momenta, muR, muF)
!   double precision, intent(out) :: muR, muF
! end subroutine userscale

function ATLAS_hww2017(ppart) result(res)
  implicit none
  include 'types.f'
  include 'constants.f'
  include 'npart.f'
  include 'mxpart.f'
  include 'jetlabel.f'
  include 'energy.f'
  include 'jetvheto.f'

  real(dp), intent(in) :: ppart(mxpart,4)
  logical :: res

  real(dp) :: etarap,pt
  real(dp) :: y4,y5,pt3,pt4,pt5,pt6
  real(dp) :: pt34,pttwo
  real(dp) :: pt45,pt56,m45,mt45,mtrans
  real(dp) :: et_vec(4),etmiss,r2,delphi,m34,m56,m3456
  integer :: i
  integer, parameter :: VVcut=3 ! set cuts for e mu
  logical :: passcuts, passveto
  real(dp) :: etaj,ptj,ptmiss,rjl1,rjl2,r,eta4,eta5,ptll
  real(dp) :: dphi,ptrel,pt36(4)
  real(dp) :: etajveto

  !f(p1) + f(p2) --> W^-(-->nu(p3) + e^+(p4)) + W^+(-->e^-(p5) + nu~(p6))

  ! for debug set etajveto to large
  etajveto = 99._dp

  pt3 = pt(3,ppart)
  pt4 = pt(4,ppart)
  pt5 = pt(5,ppart)
  pt6 = pt(6,ppart)
  y4 = etarap(4,ppart)
  y5 = etarap(5,ppart)
  pt45 = pttwo(4,5,ppart)
  pt34 = pttwo(3,4,ppart)
  pt56 = pttwo(5,6,ppart)
  r2 = (ppart(4,1)*ppart(5,1)+ppart(4,2)*ppart(5,2)) &
       /sqrt((ppart(4,1)**2+ppart(4,2)**2)*(ppart(5,1)**2+ppart(5,2)**2))
  if (r2 > +0.9999999_dp) r2=+1._dp
  if (r2 < -0.9999999_dp) r2=-1._dp
  delphi=acos(r2)
!--- mll cut
  m45 = (ppart(4,4) + ppart(5,4))**2
  do i = 1, 3
     m45 = m45 - (ppart(4,i) + ppart(5,i))**2
  enddo
  m45 = sqrt(m45)
  m34 = (ppart(3,4) + ppart(4,4))**2
  m3456 = (ppart(3,4) + ppart(4,4) + ppart(5,4) + ppart(6,4))**2
  do i = 1, 3
     m34 = m34 - (ppart(3,i) + ppart(4,i))**2
     m3456 = m3456 - (ppart(3,i) + ppart(4,i) + ppart(5,i) + ppart(6,i))**2
  enddo
  m34 = sqrt(m34)
  m3456 = sqrt(m3456)
  m56 = (ppart(5,4) + ppart(6,4))**2
  do i = 1, 3
     m56 = m56 - (ppart(5,i) + ppart(6,i))**2
  enddo
  m56 = sqrt(m56)

  mt45 = zero
  mt45 = (sqrt(sqrt(pttwo(4,5,ppart)**2 + m45**2) + etmiss(ppart,et_vec))**2)

!     PFM & GZ place cuts for WW
  passcuts = .true.
  passveto = .true.

! define ptmiss
  ptmiss = pttwo(3,6,ppart)

! define ptrel
  pt36 = ppart(3, :) + ppart(6, :)
  ptrel = ptmiss
  dphi = pi

! lepton 1
  dphi = min(dphi, &
       acos((ppart(4,1)*pt36(1) + ppart(4,2)*pt36(2))/pt4/ptmiss))
! lepton 2
  dphi = min(dphi, &
       acos((ppart(5,1)*pt36(1) + ppart(5,2)*pt36(2))/pt5/ptmiss))
! jets
  ! if (jets > 0) then
  !    dphi = min(dphi, &
  !         acos((ppart(7,1)*pt36(1) + ppart(7,2)*pt36(2))/pt(7,ppart)/ptmiss))
  ! endif

  if (dphi < pi/two) then
     ptrel = ptrel * sin(dphi)
  endif

! define ptj, etaj, rjl1, rjl2
  if (jets > 0) then
     ptj = sqrt(ppart(7,1)**2+ppart(7,2)**2)
     etaj = etarap(7,ppart)
     rjl1 = r(ppart,4,7)
     rjl2 = r(ppart,5,7)
  endif

! rapidities of leptons
  eta4 = etarap(4,ppart)
  eta5 = etarap(5,ppart)

! pt of leptons
  ptll = pttwo(4,5,ppart)

!==========================================
!     13 TeV cuts ATLAS CONF 016-090

  if (abs(sqrts - 13000._dp) < 1._dp) then

     ! ptminl = 25._dp

     ! etamaxmu = 2.4_dp
     ! etamaxe = 2.47_dp
     ! etaegapmin = 1.37_dp
     ! etaegapmax = 1.52_dp

     ! mllmin = 10._dp
     ! metrelmin = 15._dp

     if (VVcut == 3) then
        ! emu
        if (abs(eta5) > 2.5_dp) passcuts = .false.
        if (abs(eta4) > 2.5_dp) passcuts = .false.
        if (m45 < 10._dp) passcuts = .false.
        if (ptmiss < 20._dp) passcuts = .false.
        if (ptrel < 15._dp) passcuts = .false.

        if (pt(4,ppart) < 25._dp) passcuts=.false.
        if (pt(5,ppart) < 25._dp) passcuts=.false.

        if (jets > 0) then
           if (ptj > ptj_veto .or. abs(etaj) > etajveto) &
                passveto = .false.
        endif

     elseif (VVcut == 4) then
        ! mue
        if (abs(eta5) > 2.5_dp) passcuts = .false.
        if (abs(eta4) > 2.5_dp) passcuts = .false.
        if (m45 < 10._dp) passcuts = .false.
        if (ptmiss < 20._dp) passcuts = .false.
        if (ptrel < 15._dp) passcuts = .false.

        if (pt(4,ppart) < 25._dp) passcuts=.false.
        if (pt(5,ppart) < 25._dp) passcuts=.false.

        if (jets > 0) then
           if (ptj > ptj_veto .or. abs(etaj) > etajveto) &
                passveto = .false.
        endif

     elseif (VVcut == 5) then
        ! only jet veto cuts
        if (jets > 0) then
           if (ptj > ptj_veto) passcuts = .false.
        endif

     else
        stop 'VVcuts not set'
     endif

!==========================================
!     8 TeV cuts ATLAS CONF 014-033

  elseif (abs(sqrts-8000._dp) < 1._dp) then


     if (VVcut == 1) then
        ! mumu
        if (abs(eta4) > 2.4_dp) passcuts = .false.
        if (abs(eta5) > 2.4_dp) passcuts = .false.
        if (m45 < 15._dp) passcuts = .false.
        if (abs(m45-91.188_dp) < 15._dp) passcuts = .false.
        if (ptmiss < 45._dp) passcuts = .false.

        if (max(pt(4,ppart), pt(5,ppart)) < 25._dp) passcuts = .false.
        if (min(pt(4,ppart), pt(5,ppart)) < 20._dp) passcuts = .false.

        if (jets > 0) then
           if (ptj > ptj_veto .or. abs(etaj) > etajveto) &
                passveto = .false.
        endif
        if (ptrel < 45._dp) passcuts = .false.

     elseif (VVcut == 2) then
        ! ee
        if (abs(eta4) > 2.47_dp) passcuts = .false.
        if (abs(eta5) > 2.47_dp) passcuts = .false.
        if (abs(eta4) > 1.37_dp .and. abs(eta4) < 1.52_dp) &
             passcuts = .false.
        if (abs(eta5) > 1.37_dp .and. abs(eta5) < 1.52_dp) &
             passcuts = .false.
        if (m45 < 15._dp) passcuts = .false.
        if (abs(m45-91.188_dp) < 15._dp) passcuts = .false.
        if (ptmiss < 45._dp) passcuts = .false.

        if (max(pt(4,ppart), pt(5,ppart)) < 25._dp) passcuts = .false.
        if (min(pt(4,ppart), pt(5,ppart)) < 20._dp) passcuts = .false.

        if (jets > 0) then
           if (ptj > ptj_veto .or. abs(etaj) > etajveto .and. &
                rjl1 > 0.3_dp .and. rjl2 > 0.3_dp) passveto = .false.
        endif
        if (ptrel < 45._dp) passcuts = .false.

     elseif (VVcut == 3) then
        ! emu
        if (abs(eta5) > 2.4_dp) passcuts = .false.
        if (abs(eta4) > 2.47_dp) passcuts = .false.
        if (abs(eta4) > 1.37_dp .and. abs(eta4) < 1.52_dp) &
             passcuts = .false.
        if (m45 < 10._dp) passcuts = .false.
        if (ptmiss < 20._dp) passcuts = .false.

        if (max(pt(4,ppart), pt(5,ppart)) < 25._dp) passcuts = .false.
        if (min(pt(4,ppart), pt(5,ppart)) < 20._dp) passcuts = .false.

        if (jets > 0) then
           if (ptj > ptj_veto .or. abs(etaj) > etajveto .and. &
                rjl1 > 0.3_dp) passveto = .false.
        endif
        if (ptrel < 15._dp) passcuts = .false.

     elseif (VVcut == 4) then
        ! mue
        if (abs(eta4) > 2.4_dp) passcuts = .false.
        if (abs(eta5) > 2.47_dp) passcuts = .false.
        if (abs(eta5) > 1.37_dp .and. abs(eta5) < 1.52_dp) &
             passcuts = .false.
        if (m45 < 10._dp) passcuts = .false.
        if (ptmiss < 20._dp) passcuts = .false.

        if (max(pt(4,ppart), pt(5,ppart)) < 25._dp) passcuts = .false.
        if (min(pt(4,ppart), pt(5,ppart)) < 20._dp) passcuts = .false.

        if (jets > 0) then
           if (ptj > ptj_veto .or. abs(etaj) > etajveto .and. &
                rjl2 > 0.3_dp) passveto = .false.
        endif
        if (ptrel < 15._dp) passcuts = .false.

     elseif (VVcut == 5) then
        ! only jet veto cuts
        if (jets > 0) then
           if (ptj > ptj_veto) passcuts = .false.
        endif

     else
        stop 'VVcuts not set'
     endif
!==========================================
!     7 Tev cuts 1210.2979

  elseif (abs(sqrts-7000._dp) < 1._dp) then


     if (VVcut == 1) then
        ! mumu
        if (abs(eta4) > 2.4_dp) passcuts = .false.
        if (abs(eta5) > 2.4_dp) passcuts = .false.
        if (m45 < 15._dp) passcuts = .false.
        if (abs(m45-91.188_dp) < 15._dp) passcuts = .false.
        if (ptll < 30._dp) passcuts = .false.

        if (max(pt(4,ppart), pt(5,ppart)) < 25._dp) passcuts = .false.
        if (min(pt(4,ppart), pt(5,ppart)) < 20._dp) passcuts = .false.

        if (jets > 0) then
           if (ptj > ptj_veto .or. abs(etaj) > etajveto) &
                passveto = .false.
        endif
        if (ptrel < 45._dp) passcuts = .false.

     elseif (VVcut == 2) then
        ! ee
        if (abs(eta4) > 2.47_dp) passcuts = .false.
        if (abs(eta5) > 2.47_dp) passcuts = .false.
        if (abs(eta4) > 1.37_dp .and. abs(eta4) < 1.52_dp) &
             passcuts = .false.
        if (abs(eta5) > 1.37_dp .and. abs(eta5) < 1.52_dp) &
             passcuts = .false.
        if (m45 < 15._dp) passcuts = .false.
        if (abs(m45-91.188_dp) < 15._dp) passcuts = .false.
        if (ptll < 30._dp) passcuts = .false.

        if (max(pt(4,ppart), pt(5,ppart)) < 25._dp) passcuts = .false.
        if (min(pt(4,ppart), pt(5,ppart)) < 20._dp) passcuts = .false.

        if (jets > 0) then
           if (ptj > ptj_veto .or. abs(etaj) > etajveto .and. &
                rjl1 > 0.3_dp .and. rjl2 > 0.3_dp ) &
                passveto = .false.
        endif
        if (ptrel < 45._dp) passcuts = .false.

     elseif (VVcut == 3) then
        ! emu
        if (abs(eta5) > 2.4_dp) passcuts = .false.
        if (abs(eta4) > 2.47_dp) passcuts = .false.
        if (abs(eta4) > 1.37_dp .and. abs(eta4) < 1.52_dp) &
             passcuts = .false.
        if (m45 < 10._dp) passcuts = .false.
        if (ptll < 30._dp) passcuts = .false.

        if (max(pt(4,ppart), pt(5,ppart)) < 25._dp) passcuts = .false.
        if (min(pt(4,ppart), pt(5,ppart)) < 20._dp) passcuts = .false.

        if (jets > 0) then
           if (ptj > ptj_veto .or. abs(etaj) > etajveto .and. &
                rjl1 > 0.3_dp) passveto = .false.
        endif
        if (ptrel < 25._dp) passcuts = .false.

     elseif (VVcut == 4) then
        ! mue
        if (abs(eta4) > 2.4_dp) passcuts = .false.
        if (abs(eta5) > 2.47_dp) passcuts = .false.
        if (abs(eta5) > 1.37_dp .and. abs(eta5) < 1.52_dp) &
             passcuts = .false.
        if (m45 < 10._dp) passcuts = .false.
        if (ptll < 30._dp) passcuts = .false.

        if (max(pt(4,ppart), pt(5,ppart)) < 25._dp) passcuts=.false.
        if (min(pt(4,ppart), pt(5,ppart)) < 20._dp) passcuts=.false.

        if (jets > 0) then
           if (ptj > ptj_veto .or. abs(etaj) > etajveto .and. &
                rjl2 > 0.3_dp) passveto = .false.
        endif
        if (ptrel < 25._dp) passcuts = .false.

     elseif (VVcut == 5) then
        ! only jet veto cuts
        if (jets > 0) then
           if (ptj > ptj_veto) passcuts = .false.
        endif

     else
        stop 'VVcuts not set'
     endif
  else
     stop 'Energy not OK for this study'
  endif

  if (passveto .and. passcuts) then
     res = .true.
  else
     res = .false.
  endif

end function
