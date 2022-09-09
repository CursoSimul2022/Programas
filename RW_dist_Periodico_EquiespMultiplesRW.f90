module random
   implicit none
   integer::q1=1,q2=104
   integer ir(670)
   integer::sem(670)

   real(8)::nmax=2*(2**30-1)+1
   real(8) sig
   real(8) r
end module random

module globals
   implicit none
   integer, allocatable, dimension(:)::der,izq   
end module globals

module temporal
   implicit none
   integer tptos
   !integer tmax
   real(8),allocatable,dimension(:)::time
end module temporal


Program Numeros
   use random
   use temporal
   use globals
   implicit none

   integer rea,k,i,nrea   
   integer L,N,x0,nw,sel,pos
   integer tmedir
   
   real(8) tt
   real(8) rho  !densidad de particulas o caminantes. Consideraremos que 0<rho<=1.
   

   integer, allocatable, dimension(:)::position !posicion de caminantes.
   
   real(8), allocatable, dimension(:)::rho_t

   
   L    =1000   !Longitud de la red o cuadricula 1D en la que se mueve
   N    =10000   !Nro de pasos que realizar'a cada caminante
   rho  =0.1d0   !Densidad de caminantes   
   nrea =100   !Nro de realizaciones
   
   allocate(position(L))  !Como rho<=1, el numero maximo de caminantes ser'a L
   
   call initialize_random   
   call tiempo(N)  
   call bc(L) 
   
   allocate(rho_t(0:tptos))
   
   rho_t =0d0
   
   realiz: do rea=1,nrea
   
      if(mod(rea,10)==0)print*,rea
      !condicion inicial
      tmedir   =0
      nw       =0
      position =0
      do i=1,L
         call rand
         if(r<rho)then
            nw           =nw+1
            position(nw) =i
         endif
      enddo      
      
      tt       =0
   
      temp: do while(tt<=N)
         call rand
         sel  =nw*r+1
         pos  =position(sel)
         
         call rand
         if(r<0.5d0)then
            position(sel)  =der(pos)
         else
            position(sel)  =izq(pos)
         endif
         tt   =tt+1d0/nw
        
         if(tt>=time(tmedir))then    
            rho_t(tmedir)  =rho_t(tmedir)+nw/dble(L)
            tmedir         =tmedir+1
            if(tmedir>tptos) exit
         endif         
      enddo temp         
   enddo realiz
   
   open(1,file="Rho_Evol.dat")
      do k=0,tptos
         write(1,*) time(k),rho_t(k)/nrea
      enddo
   close(1)      
 
end Program Numeros

!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tiempo(N)
  use temporal
  implicit none
  integer N,i
  
  !equiespaciado en potencias de base b
  tptos=log(1.d0*N)/log(2d0)
  allocate(time(0:tptos))
  do i=0,tptos
     time(i)=2d0**i
  enddo
end subroutine tiempo


!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bc(L)
  use globals
  implicit none
  integer i,L

  allocate(der(L),izq(L))
  
  do i=1,L
    der(i)   =i+1
    izq(i)   =i-1
  enddo
   
  der(L)=1
  izq(1)=L

end subroutine bc

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





