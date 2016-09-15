      subroutine calcLtilde(muresum,Ltilde)
      implicit none

      include 'resumscale.f'
      include 'ptveto.f'
c      include 'ptvetomax.f'
      double precision muresum, jetpow, Ltilde
      
      jetpow = 5d0
      Ltilde=log( (muresum/ptveto)**jetpow + 1 ) /jetpow 
      
      return
      end
