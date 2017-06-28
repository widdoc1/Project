      subroutine check_resm_order(resm_order)
c--- checks the value of nproc and part against a list of processes that
c--- are calculable only at LO. If the calculation is not possible,
c--- writes an error message and aborts
      implicit none
      include 'frag.f'
      include 'part.f'
      include 'nproc.f'
      character*4 resm_order
      
c--- catch if the resummation cannot be performed
      if (resm_order .eq. '') then 
        write(6,*)
        write(6,*)'This process cannot be resummed,'
        write(6,*)'it is not a colour singlet final state'
      stop
      endif  

c--- if the resummation exists and we're calculating it 
c--- at leading order in alpha_s
      if (part .eq. 'LL') return
      if (part .eq. 'NLL') return

c--- otherwise, we're performing a resummation at next to 
c--- leading order alpha_s, or NNLL resummation 
      if (resm_order .ne. 'NNLL') then 
        write(6,*)
        write(6,*)'This process cannot be resummed beyond NLL - please'
        write(6,*)'check the values of nproc and part then try again'
        stop
      endif
       
      return
      end
  
