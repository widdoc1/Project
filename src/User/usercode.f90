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
logical function userincludedipole(nd, ppart, mcfm_result)
  use types_mod
  implicit none
  include 'constants.f'
  include 'mxpart.f'
  include 'npart.f'
  integer,          intent(in) :: nd
  real(dp), intent(in) :: ppart(mxpart,4)
  logical,          intent(in) :: mcfm_result
  !------------
  logical    bin
  common/bin/bin
  double precision :: ht

  ! take the MCFM result as the default choice
  userincludedipole = mcfm_result

  ! an example of placing a cut on HT
  !ht = sum(sqrt(sum(ppart(3:2+npart,1:2)**2,dim=2)))
  !if (ht < 100d0) userincludedipole = .false.
  !if (bin) write(0,*) nd, ht, userincludedipole

  !write(0,*) nd, npart, sqrt(sum(ppart(3:2+npart,1:2)**2,dim=2))
  !if (sum(ppart(3,1:2)**2) < 400) userincludedipole = .false.
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

subroutine userplotter(pjet, wt,wt2, nd)
  use types_mod
  implicit none
  include 'constants.f'
  include 'mxpart.f'
  include 'ptilde.f'
  include 'npart.f'
  include 'nplot.f'
  double precision, intent(in)  :: pjet(mxpart,4)
  double precision, intent(in)  :: wt,wt2
  integer,          intent(in)  :: nd
  !---------------------------------------------
  integer :: iplot 
  double precision :: ht, htjet
  logical, save :: first = .true.
  character*4   :: tag
  if (first) then
    tag   = "book"
    first = .false.
  else                
    tag = "plot"
  end if
  
  iplot = nextnplot

  ht    = sum(sqrt(sum(ptilde   (nd,3:2+npart,1:2)**2,dim=2)))
  htjet = sum(sqrt(sum(ptildejet(nd,3:2+npart,1:2)**2,dim=2)))
  !write(0,*) nd, ht, htjet

  call bookplot(iplot,tag,'UserHT',ht,wt,wt2,0d0,500d0,20d0,'lin'); 
  iplot = iplot + 1
  ! 
  ! call bookplot(iplot,tag,'UserHTJet',htjet,wt,wt2,0d0,500d0,20d0,'lin')
  ! iplot = iplot + 1

end subroutine userplotter

!----------------------------------------------------------------------
! user code to write info 
subroutine userwriteinfo(unitno, comment_string, xsec, xsec_err, itno)
  use types_mod
  implicit none
  integer,          intent(in) :: unitno
  character*2,      intent(in) :: comment_string
  double precision, intent(in) :: xsec, xsec_err
  integer,          intent(in) :: itno
  
  write(6,*) "have reached iteration number", itno
  write(unitno,"(a,a)") comment_string, "any additional user comments"
  call mcfmfwrite(unitno, comment_string//"any additional user comments")
end subroutine userwriteinfo

subroutine userhistofin(xsec,xsec_err,itno,itmx)
  !	This function allows for extra user-defined operations 
  !	at the end of each iteration (itno>0) and at the end of 
  !	the run of the program (itno=0).
	implicit none
	integer itno,itmx
	double precision xsec,xsec_err
end subroutine userhistofin

! subroutine userscale(event_momenta, muR, muF)
!   double precision, intent(out) :: muR, muF
! end subroutine userscale
