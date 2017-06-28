      subroutine resumset(p)
!      (rscalestart,fscalestart,resumscalestart,p)
!    wrapper routine to set up the resummation
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
      
      character(len=5) collider
      character(len=2) process
      character(len=6) observable
      real(dp) muR, muF, muResum
      real(dp) inv_mass
      real(dp) coupling, pow_sup, jet_radius
      real(dp) p(mxpart,4), mu0
      

      
      mu0=(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
     &     -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     &     -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
     &     -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2       
        mu0=dsqrt(dabs(mu0))

!     check in scale is dynamic
!     if it is make sure to set 
      
c      if     (dynstring .eq. 'm(34)') then
c        call scaleset_m34(p,mu0)
c      elseif (dynstring .eq. 'm(345)') then
c        call scaleset_m345(p,mu0)
c      elseif (dynstring .eq. 'm(3456)') then
c        call scaleset_m3456(p,mu0)
c      elseif (dynstring .eq. 'sqrt(M^2+pt34^2)') then
c        call scaleset_Msqpt34sq(p,mu0)
c      elseif (dynstring .eq. 'sqrt(M^2+pt345^2)') then
c        call scaleset_Msqpt345sq(p,mu0)
c      elseif (dynstring .eq. 'sqrt(M^2+pt5^2)') then
c        call scaleset_Msqpt5sq(p,mu0)
c      elseif (dynstring .eq. 'sqrt(M^2+ptj1^2)') then
c        call scaleset_Msqptj1sq(p,mu0)
c      elseif (dynstring .eq. 'sqrt(M^2+sumptj^2)') then
c        call scaleset_Msqsumptjsq(p,mu0)
c      elseif (dynstring .eq. 'sqrt(m(34)^2+sumptj^2)') then
c        call scaleset_m34sqsumptjsq(p,mu0)
c      elseif (dynstring .eq. 'pt(photon)') then
c        call scaleset_ptphoton(p,mu0)
c      elseif (dynstring .eq. 'HT') then
c        call scaleset_HT(p,mu0)
c      elseif (dynstring .eq. 'DDIS') then
c        call scaleset_ddis(p,mu0)
c      elseif (dynstring .eq. 's-hat') then
c        call scaleset_shat(p,mu0)
c!      else
c!        write(6,*) 'Dynamic scale choice not recognized'
c!        write(6,*) '   dynamicscale = ',dynstring
c!        stop
c      endif
      
      collider   = 'pp'
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
      call set_process_and_parameters(cs, collider, process, inv_mass,
     &   muR, muF, coupling, muResum, pow_sup, jet_radius, observable) 

      ! initialise parameters needed in the radiator
      call init_proc(cs)

      ! rescale factorisation for PDFs
      if (part .eq. 'LL') then
        return
      else if (part .eq. 'NLL') then
        facscale = facscale * exp(-Ltilde(ptveto/cs%Q, cs%p))
      else if (part .eq. 'NNLL') then
        facscale = facscale * exp(-Ltilde(ptveto/cs%Q, cs%p))
      else
        write(*,*) 'something wrong...'
        stop
      endif

      end
      
