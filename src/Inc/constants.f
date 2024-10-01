!      - comments start with "!" rather than "C"
!      - continuation lines are avoided
   
      real(dp),parameter::zip=0._dp
      real(dp),parameter::zero=0._dp
      real(dp),parameter::one=1._dp
      real(dp),parameter::two=2._dp
      real(dp),parameter::three=3._dp
      real(dp),parameter::four=4._dp
      real(dp),parameter::five=5._dp
      real(dp),parameter:: six=6._dp
      real(dp),parameter::seven=7._dp
      real(dp),parameter::eight=8._dp
      real(dp),parameter::nine=9._dp
      real(dp),parameter::ten=10._dp
      real(dp),parameter::eleven=11._dp
      real(dp),parameter::twelve=12._dp
      real(dp),parameter::sixteen=16._dp
      real(dp),parameter::half=0.5_dp
      real(dp),parameter::quarter=0.25_dp

      real(dp),parameter::pi=3.14159265358979311599796346854418516_dp
      real(dp),parameter::pisq=pi*pi
      real(dp),parameter::pisqo6=pisq/six
      real(dp),parameter::twopi=two*pi
      real(dp),parameter::fourpi=four*pi
      real(dp),parameter::pion4=pi/four
      real(dp),parameter::pion10=pi/ten
      real(dp),parameter::pisqm8=pisq-eight

      real(dp),parameter::rt2=1.41421356237309504880168872420969798_dp
      real(dp),parameter::twort2=two*rt2
      real(dp),parameter::fourrt2=four*rt2
      real(dp),parameter::rt2onpi=0.797884560802865371431347722487877591
! sqrt(two/pi)
!-----------------------------------------------------
!-----------------------------------------------------
      complex(dp),parameter::im=(zip,one)
      complex(dp),parameter::impi=(zip,pi)
      complex(dp),parameter::czip=(zip,zip)
      complex(dp),parameter::cone=(one,zip)
      complex(dp),parameter::ctwo=(two,zip)
      complex(dp),parameter::chalf=(half,zip)
!-----------------------------------------------------


      real(dp),parameter::cf=four/three
      real(dp),parameter::ca=three
      real(dp),parameter::xn=three
      real(dp),parameter::Nc=three
      real(dp),parameter::Ncinv=one/three
      real(dp),parameter::xnsq=nine
      real(dp),parameter::v=eight
      real(dp),parameter::tr=half
      real(dp),parameter::Von4=two
      real(dp),parameter::ninth=one/nine
      real(dp),parameter::xn4=xnsq-four
      real(dp),parameter::qu=two/three
      real(dp),parameter::qd=-one/three
      real(dp),parameter::qe=-one
      real(dp),parameter::spinave=one/four
      real(dp),parameter::aveqq=spinave/xnsq
      real(dp),parameter::aveqg=spinave/xn/v
      real(dp),parameter::avegg=spinave/v**2
      real(dp),parameter::aem=one/137.035989_dp

!-----------------------------------------------------

      real(dp),parameter::nbGeV2=0.389379e6_dp
      real(dp),parameter::pbGeV2=0.389379e9_dp
      real(dp),parameter::fbGeV2=0.389379e12_dp
!----decifemtobarns
      real(dp),parameter::dfbGeV2=0.389379e13_dp
      real(dp),parameter::overa=pbGeV2/xn/256._dp/pi
      integer,parameter:: nloop=2,fn=-5

