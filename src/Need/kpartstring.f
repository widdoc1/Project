c--- converts integer kpart to corresponding 4-character string
      function kpartstring(k)
      implicit none
      include 'kpart.f'
      character*15 kpartstring
      integer k

      if     (k == klord) then
        kpartstring='lo'
      elseif (k == kvirt) then
        kpartstring='virt'
      elseif (k == kreal) then
        kpartstring='real'
      elseif (k == ktota) then
        kpartstring='nlo'
      elseif (k == kfrag) then
        kpartstring='frag'
      elseif (k == ktodk) then
        kpartstring='todk'
      elseif (k == ksnlo) then
        kpartstring='snlo'
      elseif (k == knnlo) then
        kpartstring='nnlo'
      elseif (k == kll) then
        kpartstring='ll'
      elseif (k == knll) then
        kpartstring='nll'
      elseif (k == knnll) then
        kpartstring='nnll'
      elseif (k == knllexpd) then
         kpartstring='nllexpd'
      elseif (k == knnllexpd) then
         kpartstring='nnllexpd'
      elseif (k == klumi) then
         kpartstring='lumi'
      elseif (k == klumi0) then
         kpartstring='lumi0'
      elseif (k == klumi1) then
         kpartstring='lumi1'
      else
        write(6,*) 'Unexpected kpart in kpartstring: ',k
        stop
      endif

      if (coeffonly) then
        kpartstring=trim(kpartstring)//'coeff'
      endif

      return
      end

      
