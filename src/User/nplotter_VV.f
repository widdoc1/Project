      subroutine nplotter_VV(p,wt,wt2,switch,nd)
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c--- switch:  an integer equal to 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      implicit none
      include 'vegas_common.f'
      include 'constants.f'
      include 'histo.f'
      include 'jetlabel.f'
      include 'masses.f'
      include 'outputflags.f'
      include 'nqcdjets.f'      
      include 'VVcut.f'
      include 'ptveto.f'
      
      double precision sqrts 
      common/energy/sqrts
      
      double precision p(mxpart,4),wt,wt2, wt_tmp, wt2_tmp
      double precision etarap,pt
      double precision y4,y5,pt3,pt4,pt5,pt6
      double precision pt34,pttwo
      double precision pt45,pt56,m45,mt45,mtrans
      double precision et_vec(4),etmiss,r2,delphi,m34,m56,m3456
      double precision mthl,mthu
      integer switch,n,nplotmax,nd 
      character*4 tag
      integer i

      logical, save::first=.true.
      common/nplotmax/nplotmax

      logical passcuts, passveto
      double precision etaj,ptj,ptmiss,rjl1,rjl2,r,eta4,eta5,ptll
      double precision dphi,ptrel,pt36(2) 
      integer ptbin,ptbinmin
      
ccccc!$omp threadprivate(first,/nplotmax/,y4,pt3,pt4,y5)

!=====info for "peak" transverse mass plot
      mthu=hmass*1d0
      mthl=0.75d0*hmass
      
************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************


      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag='book'
        y4=1d3
        pt3=1d7
        pt4=1d7
        pt34=1d7
c---Intiailise photon 
        y5=1d3
        pt5=1d7
        mt45=1d7
        m45=1d7
        pt45=1d7
        pt34=1d7
        pt56=1d7
c----Initialise jet values will not pass cuts in there is an NLO jet
        pt6=1d7
        jets=nqcdjets
        goto 99
      else
c--- Add event in histograms
        tag='plot'
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
************************************************************************

!     121 '  f(p1)+f(p2) --> H(--> W^+(nu(p3)+e^+(p4)) + W^-(e^-(p5)+nu~(p6))) [top, bottom loops, exact]' 'L'
!     122 '  f(p1)+f(p2) --> H(--> W^+(nu(p3)+e^+(p4)) + W^-(e^-(p5)+nu~(p6))) [above + interf. with gg->WW]' 'L'
      
         pt3=pt(3,p)
         pt4=pt(4,p) 
         pt5=pt(5,p)
         pt6=pt(6,p) 
         y4=etarap(4,p)
         y5=etarap(5,p) 
         pt45=pttwo(4,5,p)
         pt34=pttwo(3,4,p) 
         pt56=pttwo(5,6,p)
         r2= (p(4,1)*p(5,1)+p(4,2)*p(5,2))
     .        /dsqrt((p(4,1)**2+p(4,2)**2)*(p(5,1)**2+p(5,2)**2))
         if (r2 .gt. +0.9999999D0) r2=+1D0
         if (r2 .lt. -0.9999999D0) r2=-1D0
         delphi=dacos(r2)
!--- mll cut 
         m45=(p(4,4)+p(5,4))**2 
         do i=1,3
           m45=m45-(p(4,i)+p(5,i))**2
         enddo
         m45=dsqrt(m45)
         m34=(p(3,4)+p(4,4))**2 
         m3456=(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
         do i=1,3
           m34=m34-(p(3,i)+p(4,i))**2
           m3456=m3456-(p(3,i)+p(4,i)+p(5,i)+p(6,i))**2
        enddo
        m34=dsqrt(m34)
        m3456=dsqrt(m3456)
        m56=(p(5,4)+p(6,4))**2 
        do i=1,3
           m56=m56-(p(5,i)+p(6,i))**2
        enddo
        m56=dsqrt(m56)

      mt45=0d0 
      mt45=(dsqrt(dsqrt(pttwo(4,5,p)**2+m45**2)+etmiss(p,et_vec))**2)
 !    endif

C     PFM & GZ place cuts for WW
      passcuts = .true. 
      passveto = .true.

C define ptmiss 
      ptmiss = pttwo(3,6,p) 

C define ptrel 
      pt36(1) = p(3,1)+p(6,1)
      pt36(2) = p(3,2)+p(6,2)
      dphi = acos((p(4,1)*pt36(1)+p(4,2)*pt36(2))/pt4/ptmiss)
      dphi = min(dphi,
     C     acos((p(5,1)*pt36(1)+p(5,2)*pt36(2))/pt5/ptmiss))

      if (jets > 0.and.ptj > ptveto) then 
         dphi = min(dphi,
     C        acos((p(7,1)*pt36(1)+p(7,2)*pt36(2))/pt(7,p)/ptmiss))
      endif
      if (dphi > pi/2d0) then 
         ptrel = ptmiss 
      else
         ptrel = ptmiss * sin(dphi) 
      endif

C define ptj, etaj, rjl1, rjl2
      if (jets > 0) then 
         ptj = sqrt(p(7,1)**2+p(7,2)**2)
         etaj = etarap(7,p)
         rjl1 = r(p,4,7)
         rjl2 = r(p,5,7)
      endif

C rapidities of leptons      
      eta4 = etarap(4,p)
      eta5 = etarap(5,p)

C pt of leptons     
      ptll = pttwo(4,5,p) 

C     8 TEV CUTS       
C     mumu  cuts of ATLAS CONF 014-33
      if (abs(sqrts-8000.d0) < 1d0) then 
         
         if (VVcut .eq. 1) then 
            if (abs(eta4).gt.2.4) passcuts=.false.
            if (abs(eta5).gt.2.4) passcuts=.false.
            if (m45 .lt. 15) passcuts = .false. 
            if (abs(m45-91.188d0) < 15) passcuts = .false. 
            if (ptmiss < 45 ) passcuts = .false. 

            if(max(pt(4,p),pt(5,p)).lt.25) passcuts=.false.
            if(min(pt(4,p),pt(5,p)).lt.20) passcuts=.false.
            
            if (jets > 0) then 
               if (ptj > ptveto .and. abs(etaj) < 4.5) 
     C 	       	 passveto = .false. 
            endif
            if (ptrel < 45d0) passcuts = .false.


         elseif (VVcut .eq. 2) then 
C     ee  cuts of ATLAS CONF 014-33
            if (abs(eta4).gt.2.47) passcuts=.false.
            if (abs(eta5).gt.2.47) passcuts=.false.
            if (abs(eta4).gt.1.37 .and. abs(eta4).lt.1.52) 
     C           passcuts=.false. 
            if (abs(eta5).gt.1.37 .and. abs(eta5).lt.1.52) 
     C           passcuts=.false. 
            if (m45 .lt. 15) passcuts = .false. 
            if (abs(m45-91.188d0) < 15) passcuts = .false. 
            if (ptmiss < 45 ) passcuts = .false. 

            if(max(pt(4,p),pt(5,p)).lt.25) passcuts=.false.
            if(min(pt(4,p),pt(5,p)).lt.20) passcuts=.false.

            if (jets > 0) then 
               if (ptj > ptveto .and. abs(etaj) < 4.5 .and. 
     C              rjl1 > 0.3d0 .and. rjl2 > 0.3d0 ) passveto = .false. 
            endif
            if (ptrel < 45d0) passcuts = .false.



         elseif (VVcut .eq. 3) then 
C     emu  cuts of ATLAS CONF 014-33
            if (abs(eta5).gt.2.4) passcuts=.false.
            if (abs(eta4).gt.2.47) passcuts=.false.
            if (abs(eta4).gt.1.37 .and. abs(eta4).lt.1.52) 
     C           passcuts=.false. 
            if (m45 .lt. 10) passcuts = .false. 
            if (ptmiss < 20 ) passcuts = .false. 

            if(max(pt(4,p),pt(5,p)).lt.25) passcuts=.false.
            if(min(pt(4,p),pt(5,p)).lt.20) passcuts=.false.

            if (jets > 0) then 
               if (ptj > ptveto .and. abs(etaj) < 4.5 .and. 
     C              rjl1 > 0.3d0) passveto = .false. 
            endif
            if (ptrel < 15d0) passcuts = .false.

            

         elseif (VVcut .eq. 4) then 
C     mue  cuts of ATLAS CONF 014-33
            if (abs(eta4).gt.2.4) passcuts=.false.
            if (abs(eta5).gt.2.47) passcuts=.false.
            if (abs(eta5).gt.1.37 .and. abs(eta5).lt.1.52) 
     C           passcuts=.false. 
            if (m45 .lt. 10) passcuts = .false. 
            if (ptmiss < 20 ) passcuts = .false. 

            if(max(pt(4,p),pt(5,p)).lt.25) passcuts=.false.
            if(min(pt(4,p),pt(5,p)).lt.20) passcuts=.false.

            if (jets > 0) then 
               if (ptj > ptveto .and. abs(etaj) < 4.5 .and. 
     C              rjl2 > 0.3d0) passveto = .false. 
            endif
            if (ptrel < 15d0) passcuts = .false.

C   only jet veto cuts
         elseif (VVcut .eq. 5) then 
            if (jets > 0) then 
               if (ptj > ptveto) passcuts = .false. 
            endif

         else
            stop 'VVcuts not set' 
         endif
      elseif (abs(sqrts-7000.d0) < 1d0) then 
         if (VVcut .eq. 1) then 
C==========================================
C     7 Tev cuts 1210.2979 
C     mumu 
            if (abs(eta4).gt.2.4) passcuts=.false.
            if (abs(eta5).gt.2.4) passcuts=.false.
            if (m45 .lt. 15) passcuts = .false. 
            if (abs(m45-91.188d0) < 15) passcuts = .false.
            if (ptll < 30 ) passcuts = .false. 
            
            if(max(pt(4,p),pt(5,p)).lt.25) passcuts=.false.
            if(min(pt(4,p),pt(5,p)).lt.20) passcuts=.false.

            if (jets > 0) then 
               if (ptj > ptveto .and. abs(etaj) < 4.5) 
     C             passveto = .false. 
            endif
            if (ptrel < 45d0) passcuts = .false.

         elseif (VVcut .eq.2) then 
C     ee 
            if (abs(eta4).gt.2.47) passcuts=.false.
            if (abs(eta5).gt.2.47) passcuts=.false.
            if (abs(eta4).gt.1.37 .and. abs(eta4).lt.1.52) 
     C           passcuts=.false. 
            if (abs(eta5).gt.1.37 .and. abs(eta5).lt.1.52) 
     C           passcuts=.false. 
            if (m45 .lt. 15) passcuts = .false. 
            if (abs(m45-91.188d0) < 15) passcuts = .false. 
            if (ptll < 30 ) passcuts = .false. 
            
            if(max(pt(4,p),pt(5,p)).lt.25) passcuts=.false.
            if(min(pt(4,p),pt(5,p)).lt.20) passcuts=.false.

            if (jets > 0) then 
               if (ptj > ptveto .and. abs(etaj) < 4.5 .and. 
     C              rjl1 .gt. 0.3 .and. rjl2 .gt. 0.3 ) 
     C              passveto = .false. 
            endif
            if (ptrel < 45d0) passcuts = .false.


         elseif (VVcut .eq. 3) then 
C     emu 
            if (abs(eta5).gt.2.4) passcuts=.false.
            if (abs(eta4).gt.2.47) passcuts=.false.
            if (abs(eta4).gt.1.37 .and. abs(eta4).lt.1.52) 
     C           passcuts=.false. 
            if (m45 .lt. 10) passcuts = .false.
            if (ptll < 30) passcuts = .false.  

            if(max(pt(4,p),pt(5,p)).lt.25) passcuts=.false.
            if(min(pt(4,p),pt(5,p)).lt.20) passcuts=.false.

            if (jets > 0) then 
               if (ptj > ptveto .and. abs(etaj) < 4.5 .and. 
     C              rjl1 > 0.3d0) passveto = .false. 
            endif
            if (ptrel < 25d0) passcuts = .false.


         elseif (VVcut .eq. 4) then 
C     mue 
            if (abs(eta4).gt.2.4) passcuts=.false.
            if (abs(eta5).gt.2.47) passcuts=.false.
            if (abs(eta5).gt.1.37 .and. abs(eta5).lt.1.52) 
     C           passcuts=.false. 
            if (m45 .lt. 10) passcuts = .false.
            if (ptll < 30) passcuts = .false.  

            if(max(pt(4,p),pt(5,p)).lt.25) passcuts=.false.
            if(min(pt(4,p),pt(5,p)).lt.20) passcuts=.false.

            if (jets > 0) then 
               if (ptj > ptveto .and. abs(etaj) < 4.5 .and. 
     C              rjl2 > 0.3d0) passveto = .false. 
            endif
            if (ptrel < 25d0) passcuts = .false.

C   only jet veto cuts
         elseif (VVcut .eq. 5) then 
            if (jets > 0) then 
               if (ptj > ptveto) passcuts = .false. 
            endif

         else
            stop 'VVcut not set' 
         endif
      else
         stop 'Energy not OK for this study' 
      endif


C      if (.not. passcuts) return 

!     set removebr=.true. to print the branching ratio
!      write(*,*) 'branching', BrnRat
!      stop



************************************************************************
*                                                                      *
*     FILL HISTOGRAMS                                                  *
*                                                                      *
************************************************************************

c--- Call histogram routines
   99 continue

c--- Book and fill ntuple if that option is set, remembering to divide
c--- by # of iterations now that is handled at end for regular histograms
      if (creatent .eqv. .true.) then
        call bookfill(tag,p,wt/dfloat(itmx))  
      endif

c--- "n" will count the number of histograms
      n=nextnplot

c--- Syntax of "bookplot" routine is:
c
c---   call bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot)
c
c---        n:  internal number of histogram
c---      tag:  "book" to initialize histogram, "plot" to fill
c---   titlex:  title of histogram
c---      var:  value of quantity being plotted
c---       wt:  weight of this event (passed in)
c---      wt2:  weight of this event (passed in)
c---     xmin:  lowest value to bin
c---     xmax:  highest value to bin
c---       dx:  bin width
c---   llplot:  equal to "lin"/"log" for linear/log scale


C     divide out branching ratios
C      wt = wt/(0.108**2*1000)
C      wt2 = wt2/(0.108**2*1000)**2

C     inclusive cross section
      call bookplot(n,tag,'sigincl',0.5d0,wt,wt2,0d0,1d0,1d0,'lin')
      n=n+1

C     cross section with cuts on leptons
      if (passcuts) then
         wt_tmp=wt
         wt2_tmp=wt2
      else
         wt_tmp=0d0 
         wt2_tmp=0d0
      endif 
      call bookplot(n,tag,'sigcutlep',0.5d0,wt_tmp,wt2_tmp,
     C     0d0,1d0,1d0,'lin')
      n=n+1


C     cross section with cuts on jets only
      if (passveto) then
         wt_tmp=wt
         wt2_tmp=wt2
      else
         wt_tmp=0d0 
         wt2_tmp=0d0
      endif 
      call bookplot(n,tag,'sigcutjet',0.5d0,wt_tmp,wt2_tmp,
     C     0d0,1d0,1d0,'lin')
      n=n+1
      

C     cross section with cuts on both jets and leptons
      if (passveto.and.passcuts) then
         wt_tmp=wt
         wt2_tmp=wt2
      else
         wt_tmp=0d0 
         wt2_tmp=0d0
      endif 
      call bookplot(n,tag,'sigcuts',0.5d0,wt_tmp,wt2_tmp,
     C     0d0,1d0,1d0,'lin')
      n=n+1

C     rapidity lepton 1
      call bookplot(n,tag,'eta4',eta4,
     C     wt,wt2,-4d0,4d0,0.1d0,'lin')
      n=n+1

C     rapidity lepton 2
      call bookplot(n,tag,'eta5',eta5,
     C     wt,wt2,-4d0,4d0,0.1d0,'lin')
      n=n+1


C     histogram inclusive ptjet
      call bookplot(n,tag,'lnptj',log(ptj/(2d0*wmass)),
     C     wt,wt2,-6d0,4d0,0.1d0,'lin')
      n=n+1


C     histogram ptjet with lepton cuts
      if (passcuts) then
         call bookplot(n,tag,'lnptjlept',log(ptj/(2d0*wmass)),
     C        wt,wt2,-6d0,4d0,0.1d0,'lin')
      end if
      n=n+1

c--- Plot dsigma/dmWW for all values of mWW where events are likely

c    Need to remove weights that are cut or vetoed
c    so do this as above

      if (passveto.and.passcuts) then
         wt_tmp=wt
         wt2_tmp=wt2
      else
         wt_tmp=0d0 
         wt2_tmp=0d0
      endif 

c--- Coarse grid
      call bookplot(n,tag,'0 < m(3456) < 8000',
     & m3456,wt_tmp,wt2_tmp,0d0,8000d0,80d0,'log')
      n=n+1

c--- Finer grid
      call bookplot(n,tag,'200 < m(3456) < 1200',
     & m3456,wt_tmp,wt2_tmp,200d0,1200d0,20d0,'log')
      n=n+1

c      call bookplot(n,tag,'1100 < m(3456) < 2000',
c     & m3456,wt_tmp,wt2_tmp,1100d0,2000d0,20d0,'log')
c      n=n+1

c      call bookplot(n,tag,'pt_nu_1',pt3,wt,wt2,0d0,100d0,2.5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'pt_e_1',pt4,wt,wt2,0d0,100d0,2.5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'y_e_1',y4,wt,wt2,-4d0,4d0,0.2d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'pt_W_1',pt34,wt,wt2,0d0,100d0,2.5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'pt_e_2',pt5,wt,wt2,0d0,100d0,2.5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'y_e_2',y5,wt,wt2,-4d0,4d0,0.2d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'pt_nu_2',pt6,wt,wt2,0d0,100d0,2.5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'pt_W_2',pt56,wt,wt2,0d0,100d0,2.5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'pt_ll',pt45,wt,wt2,0d0,100d0,2.5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'m_ll',m45,wt,wt2,0d0,100d0,2d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'m_ll',m34,wt,wt2,0d0,250d0,5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'m_ll_2',m56,wt,wt2,0d0,250d0,5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'m_4l',m3456,wt,wt2,hmass-0.5d0,hmass+0.5d0,
c     & 0.05d0,'lin')
c      n=n+1
c       call bookplot(n,tag,'m_4l',m3456,wt,wt2,0d0,2000d0,20d0
c     &,'lin')
c      n=n+1
c      call bookplot(n,tag,'mt',mt45,wt,wt2,0d0,200d0,5d0,'lin')
c      n=n+1
c      call bookplot(n,tag,'delphi',delphi,wt,wt2,0d0,3.14d0,0.1d0,'lin')
c      n=n+1
  
************************************************************************
*                                                                      *
*     FINAL BOOKKEEPING                                                *
*                                                                      *
************************************************************************

c--- We have over-counted the number of histograms by 1 at this point
      n=n-1

c--- Ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(n)

c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
c
      return 
      end
