module random
   implicit none
   integer::q1=1,q2=104
   integer ir(670)
   integer::sem(670)

   real(8)::nmax=2*(2**30-1)+1
   real(8) r
end module random




Program GenerandoGauss
   use random
   implicit none

   integer rea,k,i,nrea   

   real(8) med,sig
   real(8) x,y
   
   
   character(5) Var01
   character(5) Var02 
   
   
   med   =5d0   !Media
   sig   =5d0   !Dispersion 
   nrea  =100000   !Nro de realizaciones
   
   write(Var01,'(F5.2)')med
   call zeros(Var01)
   
   write(Var02,'(F5.2)') sig
   call zeros(Var02)
   
   
   call initialize_random   


   open(1,file="Numeros_med_"//Var01//"_disp_"//Var02//".dat")
   realiz: do rea=1,nrea
      call gauss(x,y,sig,med)
      write(1,*)x      
   enddo realiz
   close(1)      
 
end Program GenerandoGauss

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gauss(x,y,sig,med)
  use random
  implicit none
  real(8) med,sig,x,y
  real(8) phi,r0
  real(8), parameter::pi=4*atan(1d0)

  call rand
  do while(r == 0d0)
     call rand
  enddo
  x=r
  
  call rand
  do while(r == 0d0)
     call rand
  enddo
  y=r
  
  phi=2.D0*pi*x
  r0=dsqrt(-dlog(y)*2.D0*sig)
  x=r0*dcos(phi)+med
  y=r0*dsin(phi)+med
end

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_random
  use random
  implicit none

  integer i,see(33)
  integer hour

  CALL SYSTEM_CLOCK(hour)		
  !hour=87647444
  see=hour
  
  CALL RANDOM_SEED(put=see)

  do i=1,670
     call random_number(r)
     sem(i)   =r*nmax
     ir(i)    =i+1
  enddo
  ir(670)=1
  return
end subroutine initialize_random

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rand
  use random
  implicit none

 1 q1=ir(q1)
  q2=ir(q2)
  sem(q1)= IEOR(sem(q1),sem(q2))
  r=dfloat(sem(q1))/nmax
  if(r==1d0) go to 1
  return
end subroutine rand


subroutine zeros(b)
 implicit none
 
 character(*):: b
 integer i
 
 do i=1,len(b)
   if (b(i:i)==' ') then  
      b(i:i)='0'
   endif
 end do
end subroutine zeros






