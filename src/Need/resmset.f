      subroutine resmset(p)
      use types_mod
      use rad_tools_mod
      use resummation_mod
      implicit none
      include 'constants.f'
      include 'scale.f'
      include 'facscale.f'
      include 'mxpart.f'
      include 'JetVHeto.f'
      include 'JetVHeto_opts.f'
      include 'dynamicscale.f'
      include 'qcdcouple.f'
      include 'noglue.f'
      include 'kpart.f'
      include 'ptjveto.f'

      real(dp) :: Rcut
      common/Rcut/Rcut
      integer :: ilomomenta
      common/ilomomenta/ilomomenta

      character(len=5) :: collider
      character(len=2) :: process
      character(len=6) :: observable
      real(dp) :: muR, muF, muResum
      real(dp) :: inv_mass
      real(dp) :: coupling, pow_sup, jet_radius
      real(dp) :: p(mxpart,4), mu0
      real(dp) :: pBorn(4)
      integer :: j
      
      ! calculate invariant mass of bosons
      ! resummation lives in born phase space
      do j=1,4
        pBorn(j) = sum(p(3:ilomomenta,j))
      enddo
      M_B2 = abs(pBorn(4)**2 - pBorn(1)**2 - pBorn(2)**2 - pBorn(3)**2)
      M_B = sqrt(M_B2)

      process    = trim(Bconf)
      inv_mass   = M_B
      muR        = scale
      muF        = facscale
      muResum    = Q_scale
      coupling   = as
      pow_sup    = 5.0_dp
      jet_radius = Rcut
      observable = 'pt_jet' 

      ! set up parameters for the radiator
      call set_process_and_parameters(resm_opts, process, inv_mass,
     &   muR, muF, coupling, muResum, pow_sup, jet_radius, observable) 

      ! initialise parameters needed in the radiator
      call init_proc(resm_opts)

      ! compute the sudakov (1+dF)e**-R
      if (pure_lumi) then
        sudakov = one
      else
        if (kpart == knnll) then
          sudakov = resummed_sigma(ptjveto,resm_opts,order_NNLL)
        elseif (kpart == knll) then
          sudakov = resummed_sigma(ptjveto,resm_opts,order_NLL)
        elseif (kpart == kll) then
          sudakov = resummed_sigma(ptjveto,resm_opts,order_LL)
        else
          sudakov = one
        endif
      endif

!     rescale factorisation scale for PDFs
!     maybe add in scales for facscale_L/H for stops?

      L_tilde = Ltilde(ptjveto/resm_opts%Q,resm_opts%p)

      if ( (kpart==knll) .or. (kpart==knnll) .or.
     &     (kpart==klumi0) .or. (kpart==klumi1) .or.
     &     (kpart==klumi) ) then
         facscaleLtilde = facscale * exp(-L_tilde)
      else
         facscaleLtilde = facscale
      endif

      end
