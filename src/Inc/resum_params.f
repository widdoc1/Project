!     assorted parameters needed for the resummation

!     parts for computing the radiator
      type(process_and_parameters) :: resm_opts
!     the modified factorisation scale needed for PDFs
!     muF -> muF*exp(-Ltilde) ~ pt_jet (max)
      real(dp) :: facscaleLtilde
      
      common/resum_params/resm_opts,facscaleLtilde
