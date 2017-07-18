      subroutine resumset(p)
      use types_mod
      use rad_tools_mod
      implicit none
      include 'constants.f'
      include 'scale.f'
      include 'facscale.f'
      include 'resumscale.f'
      include 'resum_params.f'
      include 'dynamicscale.f'
      include 'qcdcouple.f'
      include 'part.f'
      include 'ptveto.f'
      
!      type(process_and_parameters) cs
      
      real(dp) Rcut  
      common/Rcut/Rcut
      integer ilomomenta
      common/ilomomenta/ilomomenta

      character(len=5) collider
      character(len=2) process
      character(len=6) observable
      real(dp) muR, muF, muResum
      real(dp) inv_mass
      real(dp) coupling, pow_sup, jet_radius
      real(dp) p(mxpart,4), mu0
      real(dp) pb(4)
      integer j
      
      ! calculate invariant mass of bosons
      ! resummation lives in born phase space
      do j=1,4
        pb(j) = sum(p(3:ilomomenta,j))
      enddo
      mu0 = sqrt(abs(pb(4)**2 - pb(1)**2 - pb(2)**2 - pb(3)**2))
      
      ! parameters for process
      process    = 'DY'
      inv_mass   = mu0
      muR        = scale
      muF        = facscale
      muResum    = resumscale
      coupling   = as
      pow_sup    = 5.0_dp
      jet_radius = Rcut
      observable = 'pt_jet' 
      
      ! set up parameters for the radiator
      call set_process_and_parameters(resm_opts, process, inv_mass,
     &   muR, muF, coupling, muResum, pow_sup, jet_radius, observable) 

      ! initialise parameters needed in the radiator
      call init_proc(resm_opts)

      ! rescale factorisation scale for PDFs
      if (part .eq. 'LL') then
        facscaleLtilde = facscale
        return
      else if (part .eq. 'NLL') then
        facscaleLtilde = facscale * exp(-Ltilde(ptveto/resm_opts%Q,
     &                                 resm_opts%p))
      else if (part .eq. 'NNLL') then
        facscaleLtilde = facscale * exp(-Ltilde(ptveto/resm_opts%Q,
     &                                 resm_opts%p))
      else
        write(*,*) 'something wrong...'
        stop
      endif

      end
      
