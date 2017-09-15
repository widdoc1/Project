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
  use types_mod
  implicit none
  include 'constants.f'
  include 'npart.f'
  include 'mxpart.f'
  include 'jetlabel.f'
  include 'nproc.f'
  integer,  intent(in) :: nd
  real(dp), intent(in) :: ppart(mxpart,4)
  logical,  intent(in) :: mcfm_result
  logical :: userincludedipole
   !------------
  logical :: bin,makecuts,ATLAS_hww2017
  common/bin/bin
  common/makecuts/makecuts

  include 'ptjveto.f'
  real(dp) :: ptj

  ! take the MCFM result as the default choice
  userincludedipole = mcfm_result
  ! return

  if (makecuts) then
     if ( (nproc .eq. 61) .or. (nproc .eq. 66) .or. (nproc .eq. 69) &
          & .or. (nproc .eq. 123) .or. (nproc .eq. 124) .or. (nproc .eq. 125) &
          & .or. (nproc .eq. 126) ) then
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
        if (ptj > ptjveto) then
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
  use types_mod
  implicit none
  include 'constants.f'
  include 'nf.f'
  include 'mxpart.f'
  include 'ptilde.f'
  include 'npart.f'
  include 'nplot.f'
  real(dp), intent(in) :: pjet(mxpart,4)
  real(dp), intent(inout) :: wt,wt2
  integer,  intent(in) :: nd
  integer, parameter :: tagbook=1, tagplot=2
  !---------------------------------------------
  integer :: i, iplot 
  real(dp) :: m34, m45, m56, m3456
  real(dp) :: MT1, MT2
  real(dp) :: ptll, pt45(2), ptmiss, pt36(2), MTll, &
       & dphillmiss, pttwo
  real(dp) :: ht, htjet
  logical, save :: first = .true.
  integer :: tag

  if (first) then
     tag   = tagbook
     first = .false.
  else
     tag = tagplot
  endif
  iplot = nextnplot

  !ht    = sum(sqrt(sum(ptilde   (nd,3:2+npart,1:2)**2,dim=2)))
  !htjet = sum(sqrt(sum(ptildejet(nd,3:2+npart,1:2)**2,dim=2)))

  !call bookplot(iplot,tag,'UserHT',ht,wt,wt2,0d0,500d0,20d0,'lin') 
  !iplot = iplot + 1

  !call bookplot(iplot,tag,'UserHTJet',htjet,wt,wt2,0d0,500d0,20d0,'lin')
  !iplot = iplot + 1

  !define quantities to plot
  m45 = (pjet(4,4) + pjet(5,4))**2 
  do i = 1, 3
     m45 = m45 - (pjet(4,i) + pjet(5,i))**2
  enddo
  m45 = sqrt(m45)
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
  ptll = pttwo(4, 5, pjet) 
  pt45(1) = pjet(4,1) + pjet(5,1)
  pt45(2) = pjet(4,2) + pjet(5,2)

  ptmiss = pttwo(3,6,pjet) 
  pt36(1) = pjet(3,1) + pjet(6,1)

  pt36(2) = pjet(3,2) + pjet(6,2)

  MTll = sqrt(ptll**2 + m45**2)

  MT1 = (MTll + ptmiss)**2
  do i = 1, 2
     MT1 = MT1 - (pt45(i) + pt36(i))**2
  enddo
  MT1 = sqrt(MT1)

  dphillmiss = acos((pt36(1)*pt45(1)+pt36(2)*pt45(2))/ptmiss/ptll)

  MT2 = sqrt(two*ptll*ptmiss*(1-cos(dphillmiss)))

  ! m(45) dists
  call bookplot(iplot,tag,'m(45)',m45,wt,wt2,0d0,8000d0,80d0,'log')
  iplot = iplot + 1

  call bookplot(iplot,tag,'m(45)',m45,wt,wt2,0d0,1000d0,20d0,'log')
  iplot = iplot + 1
  
  call bookplot(iplot,tag,'m(45)',m45,wt,wt2,100d0,1000d0,20d0,'log')
  iplot = iplot + 1

  ! m(3456) dists
  call bookplot(iplot,tag,'m(3456)',m3456,wt,wt2,0d0,8000d0,80d0,'log')
  iplot = iplot + 1

  call bookplot(iplot,tag,'m(3456)',m3456,wt,wt2,0d0,1000d0,20d0,'log')
  iplot = iplot + 1
  
  call bookplot(iplot,tag,'m(3456)',m3456,wt,wt2,100d0,1000d0,20d0,'log')
  iplot = iplot + 1

  ! MT1
  call bookplot(iplot,tag,'MT1',MT1,wt,wt2,0d0,8000d0,80d0,'log')
  iplot = iplot + 1

  call bookplot(iplot,tag,'MT1',MT1,wt,wt2,0d0,1000d0,20d0,'log')
  iplot = iplot + 1

  call bookplot(iplot,tag,'MT1',MT1,wt,wt2,100d0,1000d0,20d0,'log')
  iplot = iplot + 1

  ! MT2
  call bookplot(iplot,tag,'MT2',MT2,wt,wt2,0d0,8000d0,80d0,'log')
  iplot = iplot + 1

  call bookplot(iplot,tag,'MT2',MT2,wt,wt2,0d0,1000d0,20d0,'log')
  iplot = iplot + 1

  call bookplot(iplot,tag,'MT2',MT2,wt,wt2,100d0,1000d0,20d0,'log')
  iplot = iplot + 1

  ! update nextnplot so we get userplots and generic plots from nplotter routines
  nextnplot = iplot

end subroutine userplotter

!----------------------------------------------------------------------
! user code to write info 
subroutine userwriteinfo(unitno, comment_string, xsec, xsec_err, itno)
  use types_mod
  implicit none
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
  use types_mod
	implicit none
  integer,  intent(in) :: itno,itmx
	real(dp), intent(in) :: xsec,xsec_err
end subroutine userhistofin

! subroutine userscale(event_momenta, muR, muF)
!   double precision, intent(out) :: muR, muF
! end subroutine userscale

function ATLAS_hww2017(ppart) result(res)
  use types_mod
  implicit none
  include 'constants.f'
  include 'npart.f'
  include 'mxpart.f'
  include 'jetlabel.f'
  include 'energy.f'
  include 'ptjveto.f'

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
  real(dp) etaj,ptj,ptmiss,rjl1,rjl2,r,eta4,eta5,ptll
  real(dp) dphi,ptrel,pt36(2) 


  !f(p1) + f(p2) --> W^-(-->nu(p3) + e^+(p4)) + W^+(-->e^-(p5) + nu~(p6))

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
  pt36(1) = ppart(3,1) + ppart(6,1)
  pt36(2) = ppart(3,2) + ppart(6,2)
  dphi = acos((ppart(4,1) * pt36(1) + ppart(4,2)*pt36(2))/pt4/ptmiss)
  dphi = min(dphi, &
       acos((ppart(5,1)*pt36(1) + ppart(5,2)*pt36(2))/pt5/ptmiss))

  if (jets > 0 .and. ptj > ptjveto) then 
     dphi = min(dphi, &
          acos((ppart(7,1)*pt36(1)+ppart(7,2)*pt36(2))/pt(7,ppart)/ptmiss))
  endif
  if (dphi > pi/two) then 
     ptrel = ptmiss 
  else
     ptrel = ptmiss * sin(dphi) 
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

     if (VVcut .eq. 3) then
        ! emu
        if (abs(eta5) > 2.4_dp) passcuts=.false.
        if (abs(eta4) > 2.47_dp) passcuts=.false.
        if (abs(eta4) > 1.37_dp .and. abs(eta4) < 1.52_dp) &
             passcuts=.false.
        if (m45 < 10.) passcuts = .false.
        if (ptmiss < 20. ) passcuts = .false.
        if (ptrel < 15d0) passcuts = .false.

        if (pt(4,ppart) < 25d0) passcuts=.false.
        if (pt(5,ppart) < 25d0) passcuts=.false.

        if (jets > 0) then
           if (ptj > ptjveto .and. abs(etaj) < 4.5) &
                passveto = .false.
        endif

     elseif (VVcut .eq. 4) then 
        ! mue 
        if (abs(eta5) > 2.4) passcuts=.false.
        if (abs(eta4) > 2.47) passcuts=.false.
        if (abs(eta4) > 1.37 .and. abs(eta4) < 1.52) &
             passcuts=.false. 
        if (m45 < 10) passcuts = .false. 
        if (ptmiss < 20 ) passcuts = .false.
        if (ptrel < 15d0) passcuts = .false. 

        if (pt(4,ppart) < 25d0) passcuts=.false.
        if (pt(5,ppart) < 25d0) passcuts=.false.

        if (jets > 0) then 
           if (ptj > ptjveto .and. abs(etaj) < 4.5) &
                passveto = .false. 
        endif

     elseif (VVcut .eq. 5) then 
        ! only jet veto cuts
        if (jets > 0) then 
           if (ptj > ptjveto) passcuts = .false. 
        endif

     else
        stop 'VVcuts not set' 
     endif

!==========================================
!     8 TeV cuts ATLAS CONF 014-033

  elseif (abs(sqrts-8000.d0) < 1d0) then 


     if (VVcut .eq. 1) then 
        ! mumu 
        if (abs(eta4).gt.2.4) passcuts=.false.
        if (abs(eta5).gt.2.4) passcuts=.false.
        if (m45 .lt. 15) passcuts = .false. 
        if (abs(m45-91.188d0) < 15) passcuts = .false. 
        if (ptmiss < 45 ) passcuts = .false. 

        if(max(pt(4,ppart),pt(5,ppart)).lt.25) passcuts=.false.
        if(min(pt(4,ppart),pt(5,ppart)).lt.20) passcuts=.false.

        if (jets > 0) then 
           if (ptj > ptjveto .and. abs(etaj) < 4.5) &
                passveto = .false. 
        endif
        if (ptrel < 45d0) passcuts = .false.

     elseif (VVcut .eq. 2) then 
        ! ee 
        if (abs(eta4).gt.2.47) passcuts=.false.
        if (abs(eta5).gt.2.47) passcuts=.false.
        if (abs(eta4).gt.1.37 .and. abs(eta4).lt.1.52) &
             passcuts=.false. 
        if (abs(eta5).gt.1.37 .and. abs(eta5).lt.1.52) &
             passcuts=.false. 
        if (m45 .lt. 15) passcuts = .false. 
        if (abs(m45-91.188d0) < 15) passcuts = .false. 
        if (ptmiss < 45 ) passcuts = .false. 

        if(max(pt(4,ppart),pt(5,ppart)).lt.25) passcuts=.false.
        if(min(pt(4,ppart),pt(5,ppart)).lt.20) passcuts=.false.

        if (jets > 0) then 
           if (ptj > ptjveto .and. abs(etaj) < 4.5 .and. &
                rjl1 > 0.3d0 .and. rjl2 > 0.3d0 ) passveto = .false. 
        endif
        if (ptrel < 45d0) passcuts = .false.

     elseif (VVcut .eq. 3) then 
        ! emu 
        if (abs(eta5).gt.2.4) passcuts=.false.
        if (abs(eta4).gt.2.47) passcuts=.false.
        if (abs(eta4).gt.1.37 .and. abs(eta4).lt.1.52) & 
             passcuts=.false. 
        if (m45 .lt. 10) passcuts = .false. 
        if (ptmiss < 20 ) passcuts = .false. 

        if(max(pt(4,ppart),pt(5,ppart)).lt.25) passcuts=.false.
        if(min(pt(4,ppart),pt(5,ppart)).lt.20) passcuts=.false.

        if (jets > 0) then 
           if (ptj > ptjveto .and. abs(etaj) < 4.5 .and. &
                rjl1 > 0.3d0) passveto = .false. 
        endif
        if (ptrel < 15d0) passcuts = .false.

     elseif (VVcut .eq. 4) then 
        ! mue 
        if (abs(eta4).gt.2.4) passcuts=.false.
        if (abs(eta5).gt.2.47) passcuts=.false.
        if (abs(eta5).gt.1.37 .and. abs(eta5).lt.1.52) &
             passcuts=.false. 
        if (m45 .lt. 10) passcuts = .false. 
        if (ptmiss < 20 ) passcuts = .false. 

        if(max(pt(4,ppart),pt(5,ppart)).lt.25) passcuts=.false.
        if(min(pt(4,ppart),pt(5,ppart)).lt.20) passcuts=.false.

        if (jets > 0) then 
           if (ptj > ptjveto .and. abs(etaj) < 4.5 .and. &
                rjl2 > 0.3d0) passveto = .false. 
        endif
        if (ptrel < 15d0) passcuts = .false.

     elseif (VVcut .eq. 5) then 
        ! only jet veto cuts
        if (jets > 0) then 
           if (ptj > ptjveto) passcuts = .false. 
        endif

     else
        stop 'VVcuts not set' 
     endif
!==========================================
!     7 Tev cuts 1210.2979 

  elseif (abs(sqrts-7000.d0) < 1d0) then 


     if (VVcut .eq. 1) then 
        ! mumu 
        if (abs(eta4).gt.2.4) passcuts=.false.
        if (abs(eta5).gt.2.4) passcuts=.false.
        if (m45 .lt. 15) passcuts = .false. 
        if (abs(m45-91.188d0) < 15) passcuts = .false.
        if (ptll < 30 ) passcuts = .false. 

        if(max(pt(4,ppart),pt(5,ppart)).lt.25) passcuts=.false.
        if(min(pt(4,ppart),pt(5,ppart)).lt.20) passcuts=.false.

        if (jets > 0) then 
           if (ptj > ptjveto .and. abs(etaj) < 4.5) &
                passveto = .false. 
        endif
        if (ptrel < 45d0) passcuts = .false.

     elseif (VVcut .eq.2) then 
        ! ee 
        if (abs(eta4).gt.2.47) passcuts=.false.
        if (abs(eta5).gt.2.47) passcuts=.false.
        if (abs(eta4).gt.1.37 .and. abs(eta4).lt.1.52) &
             passcuts=.false. 
        if (abs(eta5).gt.1.37 .and. abs(eta5).lt.1.52) &
             passcuts=.false. 
        if (m45 .lt. 15) passcuts = .false.
        if (abs(m45-91.188d0) < 15) passcuts = .false.
        if (ptll < 30 ) passcuts = .false. 

        if(max(pt(4,ppart),pt(5,ppart)).lt.25) passcuts=.false.
        if(min(pt(4,ppart),pt(5,ppart)).lt.20) passcuts=.false.

        if (jets > 0) then 
           if (ptj > ptjveto .and. abs(etaj) < 4.5 .and. &
                rjl1 .gt. 0.3 .and. rjl2 .gt. 0.3 ) &
                passveto = .false. 
        endif
        if (ptrel < 45d0) passcuts = .false.

     elseif (VVcut .eq. 3) then 
        ! emu 
        if (abs(eta5).gt.2.4) passcuts=.false.
        if (abs(eta4).gt.2.47) passcuts=.false.
        if (abs(eta4).gt.1.37 .and. abs(eta4).lt.1.52) &
             passcuts=.false. 
        if (m45 .lt. 10) passcuts = .false.
        if (ptll < 30) passcuts = .false.  

        if(max(pt(4,ppart),pt(5,ppart)).lt.25) passcuts=.false.
        if(min(pt(4,ppart),pt(5,ppart)).lt.20) passcuts=.false.

        if (jets > 0) then 
           if (ptj > ptjveto .and. abs(etaj) < 4.5 .and. &
                rjl1 > 0.3d0) passveto = .false. 
        endif
        if (ptrel < 25d0) passcuts = .false.

     elseif (VVcut .eq. 4) then 
        ! mue 
        if (abs(eta4).gt.2.4) passcuts=.false.
        if (abs(eta5).gt.2.47) passcuts=.false.
        if (abs(eta5).gt.1.37 .and. abs(eta5).lt.1.52) &
             passcuts=.false. 
        if (m45 .lt. 10) passcuts = .false.
        if (ptll < 30) passcuts = .false.  

        if(max(pt(4,ppart),pt(5,ppart)).lt.25) passcuts=.false.
        if(min(pt(4,ppart),pt(5,ppart)).lt.20) passcuts=.false.

        if (jets > 0) then 
           if (ptj > ptjveto .and. abs(etaj) < 4.5 .and. &
                rjl2 > 0.3d0) passveto = .false. 
        endif
        if (ptrel < 25d0) passcuts = .false.

     elseif (VVcut .eq. 5) then 
        ! only jet veto cuts
        if (jets > 0) then 
           if (ptj > ptjveto) passcuts = .false. 
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
