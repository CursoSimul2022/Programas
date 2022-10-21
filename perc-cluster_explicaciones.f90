module global
implicit none
save
integer,allocatable,dimension(:,:)::latt
integer,allocatable,dimension(:)::mas1,menos1
integer,allocatable,dimension(:)::posx,posy
integer,allocatable,dimension(:)::cluster
integer,allocatable,dimension(:)::mass
integer max_mass
integer npoints,L
integer s2
end module global

module random
   implicit none
   integer::q1=1,q2=104
   integer ir(670)
   integer::sem(670)

   real(8)::nmax=2*(2**30-1)+1
   real(8) r
end module random



program perc
use global
use random
implicit none

INTEGER see(20),hour
integer rea,nrea
integer i,j,k
integer x,y
integer label
real(8) p
integer cluster_number
character(80) output
integer nptos
real(8) pmin,pmax,dp
integer jk
integer,allocatable,dimension(:)::histo_mass
real(8),allocatable,dimension(:)::Pinf,Pi,pp,Smed

write(*,*) 'nrea?'
read(*,*) nrea
write(*,*) 'L?'!longitud de la red
read(*,*) L
write(*,*) 'pmin?'
read(*,*) pmin
write(*,*) 'pmax?'
read(*,*) pmax
write(*,*) 'dp?'
read(*,*) dp
write(*,*) 'output file? '
read(*,*) output

call initialize_random  


allocate(mas1(L),menos1(L))
allocate(latt(0:L+1,0:L+1),posx(L*L),posy(L*L))
allocate(cluster(L*L),mass(L*L))

call bc

nptos=(pmax-pmin)/dp+1
allocate(pp(nptos),Pinf(nptos),Smed(nptos))
do i=1,nptos
   pp(i)=pmin
   pmin=pmin+dp
enddo
!if(pmax==1) pp(nptos)=0.9999999999

Pinf=0
Smed=0
ite:do rea=1,nrea
   do jk=1,nptos
      p=pp(jk)
      label=0
      latt=0 !aca voy a guardar  LxL ptos ocupados
      npoints=0 !que empiezan en 0
      do i=1,L
         do j=1,L
            call rand
            if(r <=p) then
               npoints=npoints+1 
               latt(i,j)=npoints
               posx(npoints)=i !aca sestan las posiciones x e y de los sitios ocupados
               posy(npoints)=j
            endif
         enddo
      enddo
      call cluster_id(cluster_number) !identifico los clusters
      do i=1,npoints
         x=posx(i) !digo su posicion en x 
         y=posy(i) !digo su posicion en y
         latt(x,y)=cluster(i) ! a cada elemento de la matriz le signo al cluster al cual pertenece
      enddo
      call spanning(label) !se fija si expande, label sale o 0 o con el numero de cluster al cual esta asociado ese sitio
      if(label.ne.0)  Pinf(jk)=Pinf(jk)+mass(label) !calculo el tamanio de la componente que expande
      do i=1,cluster_number
         if(i==label)cycle
         Smed(jk)  =Smed(jk)+mass(i)**2d0         
      enddo
   enddo
   if(mod(rea,10)==0) then !grabo cada 10
      open(2,file=output)
      print*,rea
      do i=1,nptos
         write(2,*) pp(i), Pinf(i)/rea/L/L,Smed(i)/rea/L/L
      enddo
      close(2)
   endif
enddo ite

stop
end program perc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spanning(label)
use global
implicit none
integer i,j,label
integer,allocatable,dimension(:)::aux
allocate(aux(L))

do i=1,L
   if(latt(i,1).ne.0) then 
      aux(i)=latt(i,1)
  endif
enddo
do i=1,L
   do j=1,L
      if(latt(j,L).ne.0.and.aux(i)==latt(j,L)) then
         !print*,' percolo en el eje vertical' 
         label=aux(i)
         exit
      endif
   enddo
enddo
deallocate(aux)
return
end subroutine spanning

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cluster_id(cluster_number)
use global
implicit none
integer,allocatable,dimension(:)::ocupp,w
integer cluster_number,nburn,cs
integer i,j,k,m,x,y,xx,yy
allocate(ocupp(npoints),w(0:npoints+1))

!!! Variables: 
!!!! Cluster_number= # de clusters
!!! mass = numero de sitios en el cluster = cluster_number
!!!! cluster= cluster al que pertenece el punto i=1,npoints

cluster_number=0
cluster=0
ocupp=1 !No esta nada asignado
mass=0
max_mass=0
do i=1,npoints
   if(cluster(i) .eq. 0.and.ocupp(i)==1) then !No le fue asignado un cluster y no fue visitado
      cluster_number=cluster_number+1 ! se le asigna el numero de cluster
      w(0)=i                          !primer sitio quemado 
      nburn=1
      cs=0                            !pongo cs=masa en 0
      do while(nburn.ne.0)
         nburn=nburn-1                !contador para la lista a quemar notar que empieza del maximo valor de nburn
         j=w(nburn)                   !j = sitio quemado
         cluster(j)=cluster_number    !le asigno el cluster #
         cs=cs+1                      !una vez que lo quemo cs - cs+1
         x=posx(j)                    !coordenada x
         y=posy(j)                    !coordenada y
         xx=mas1(x)                    !miro vecinos
         yy=y
         m=latt(xx,yy)               !m es el punto numero latt(xx,yy)
         if(m.ne.0) call vecinos(m,ocupp,w,nburn)     !si m/=0 (o sea que esta ocupado llamo a vecinos
         xx=menos1(x)
         m=latt(xx,yy)
         if(m.ne.0) call vecinos(m,ocupp,w,nburn) 
         xx=x
         yy=mas1(y)
         m=latt(xx,yy)
         if(m.ne.0) call vecinos(m,ocupp,w,nburn)   
         yy=menos1(y)
         m=latt(xx,yy)
         if(m.ne.0) call vecinos(m,ocupp,w,nburn) 
      enddo !este loop termina cuando no incorporo más sitios a la lista de quemados
      mass(cluster_number)=cs      
      if(cs > max_mass) max_mass=cs !cluster + grande
   endif
enddo
!deallocate(ocupp,w)
return
end subroutine cluster_id

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vecinos(m,ocupp,w,nburn)
use global
implicit none
integer ocupp(npoints),w(0:npoints+1)
integer m,nburn
if(cluster(m) == 0.and.ocupp(m)==1) then !si al punto m no le asigne 1 cluster y no lo visite
   ocupp(m)=0                            !digo que lo visite
   w(nburn)=m                            !lo pongo en la lista a quemar
   nburn=nburn+1                         !aumento el contador en uno
endif
return
end subroutine vecinos


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bc
use global
implicit none

integer i
do i=1,L
   mas1(i)  =i+1
   menos1(i)=i-1
enddo
return
end subroutine bc
!!!!!!!!!!!!!!!!!!!!!!!!!!!

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










