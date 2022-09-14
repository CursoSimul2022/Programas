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


Program rw
   use random
   use temporal
   use globals
   implicit none

   integer rea,k,i,nrea   
   integer L,x0,nw,sel,pos,posN
   integer tmedir
   
   real(8) tt
   real(8) rho  !densidad de particulas o caminantes. Consideraremos que 0<rho<=1.
   
   
   character(4) Var01
   character(5) Var02 
   

   integer, allocatable, dimension(:)::position,Lat !posicion de caminantes.
   
   real(8), allocatable, dimension(:)::rho_t
   
   L    =10   !Longitud de la red o cuadricula 1D en la que se mueve
   rho  =0.40d0   !Densidad de caminantes   
   nrea =100   !Nro de realizaciones
   
   write(Var01,'(F4.2)')rho
   call zeros(Var01)
   
   write(Var02,'(I5)') L
   call zeros(Var02)
   
   allocate(position(L),Lat(L))  !Como rho<=1, el numero maximo de caminantes ser'a L
   
   call initialize_random   
   call tiempo 
   call bc(L) 
   
   allocate(rho_t(0:tptos))
   
   rho_t =0d0
   
   realiz: do rea=1,nrea
   
      if(mod(rea,10)==0)print*,rea
      !condicion inicial
      tmedir   =0
      nw       =0
      position =0
      Lat      =0
      do i=1,L
         call rand
         if(r<rho)then
            nw           =nw+1
            position(nw) =i
            Lat(i)       =nw   
         endif
      enddo      
      
      tt       =0

      temp: do
         call rand
         sel  =nw*r+1
         pos  =position(sel)
         
         call rand
         if(r<0.5d0)then
            posN=der(pos)
            if(Lat(posN)==0)then
               Lat(pos)       =0
               position(sel)  =posN
               Lat(posN)      =sel
            endif
         else
            posN=izq(pos)
            if(Lat(posN)==0)then
               Lat(pos)       =0
               position(sel)  =posN
               Lat(posN)      =sel
            endif               
         endif
         tt   =tt+1d0/nw
         !print*,Lat
         !read(*,*)
         if(tt>=time(tmedir))then    
            rho_t(tmedir)  =rho_t(tmedir)+nw/dble(L)
            tmedir         =tmedir+1
            if(tmedir>tptos) exit
         endif         
      enddo temp         
   enddo realiz
   
   open(1,file="Rho_Evol_rho_"//Var01//"_L_"//Var02//".dat")
      do k=0,tptos
         write(1,*) time(k),rho_t(k)/nrea
      enddo
   close(1)      
 
end Program rw

!!!!!!!!!!!!!!!!!!!!!!!!!


!****************************
subroutine tiempo
  use temporal
  implicit none
  integer i
  real(8) xmin,xmax,base,dx
  
  base=10d0
  xmin=-2 !exponente minimo
  xmax=3!exponente maximo
  dx  =0.1d0
  tptos=(xmax-xmin)/dx !numero de puntos con equiespaciado log con equ
  
  allocate(time(0:tptos))
  do i=0,tptos
    time(i)=base**xmin
    xmin=xmin+0.1d0
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






