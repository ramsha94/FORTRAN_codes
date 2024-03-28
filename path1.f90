PROGRAM mdcode
IMPLICIT NONE

real*8 :: ke0,rand,gauss,sigma
integer :: i, ii
REAL*8::dsdx,dsdy,dzdx,dzdy,s,z

integer :: steps_max,step, pfrq,npoints,ndim
real*8  :: dt,u,u1,u2,te,ke,pe
real*8  :: dummy,t0,t_inst, avg_t
real*8  :: umbr_mean,umbr_k
logical :: us,ext, lang
real*8  :: t0_s, t_s_inst,ke_s,pe_s,k_s,avg_t_s, ddt

integer :: j,ndof,ndof_s, ios, ios2, imts, nmts
real*8  :: glang, glang_s, kbt0, kbt0_s,lambda

integer :: meta_max,freq_meta,nmtd,meta_cv,umb_cv
real*8  :: width,h0,cvar_s,cvar_z,cv_ist,cv_ist1,alpha,wtdt,gausspot
logical :: meta,hilladd,restart,mts

real*8 :: lang_params(6), lang_params_s(6)

real*8, allocatable :: x(:),v(:),mass(:),s_z(:)
real*8,allocatable  :: f(:),x_i(:,:)
real*8,allocatable  :: height(:),cv_history_s(:),cv_history_z(:)
real*8, allocatable ::x_s(:),v_s(:),mass_s(:),f_harm(:), f_s(:), rnd(:), rnd_s(:), &
                      v_old(:), v_s_old(:)
real*8, parameter :: kb=3.16e-6  !Boltman constant in a.u. K^-1
real*8, parameter :: amu_to_au=1822.d0
real*8, parameter :: pi=4.d0*atan(1.d0)
integer, parameter :: nd=2, ns=2, ncv=2            !number of dimensions
character(len=20) :: myfilestat, myfilepos
integer:: ir


npoints=5
ndim=2

mts=.false.
nmts=0

print *, "pi is =", pi
print *, "kb in HK^-1 =", kb
print *, "amu to au=", amu_to_au

OPEN(50,file="reco.dat")
open( unit =20, file = 'input',status='unknown')
OPEN(unit=1,file='output')
avg_t=0.0d0
avg_t_s=0.0d0

allocate(x(nd),v(nd),s_z(nd))
allocate(x_s(ns))
allocate(v_s(ns))
allocate(mass(nd))
allocate(mass_s(ns))
allocate(f(nd))
allocate(f_s(ns))
allocate(f_harm(ns))
allocate(rnd(nd))
allocate(rnd_s(ns))
allocate(v_s_old(ns))
allocate(v_old(nd))
ALLOCATE(x_i(npoints,ndim))

read(20,*)x(1:2)
read(20,*)mass(1:2)
read(20,*)dt
read(20,*)t0
read(20,*)t0_s
read(20,*)lang
read(20,*)glang
read(20,*)glang_s
read(20,*)ext
read(20,*)us
read(20,*)meta
read(20,*)k_s
read(20,*)mass_s(1:2)
read(20,*)umbr_k
read(20,*)umbr_mean
read(20,*)steps_max
read(20,*)freq_meta
read(20,*)h0
read(20,*)width
read(20,*)wtdt !deltaT parameter in energy unit (a.u.)
read(20,*)pfrq
read(20,*)restart
read(20,*)nmts

open(31,file="prob.dat")
write(31,*) steps_max

if(nmts.gt.0)mts=.true.

write(1,*)'x=',x(1:2)
write(1,*)'mass=',mass(1:2)
write(1,*)'dt=',dt
write(1,*)'t0=',t0
write(1,*)'t0_s=',t0_s
if(lang)then 
  write(1,*)'lang is true'
else
  write(1,*)'lang is false'
end if
write(1,*)'glang=',glang
write(1,*)'glang_s=',glang_s
if(ext)then
  write(1,*)'ext is true'
else 
  write(1,*)'ext is false'
end if
if(us)then
  write(1,*)'us is true'
else 
  write(1,*)'us is false'
end if
if(meta)then
  write(1,*)'meta is true'
else 
  write(1,*)'meta is false'
end if
write(1,*)'k_s=',k_s
write(1,*)'mass_s=',mass_s(1:2)
write(1,*)'umbr_k=',umbr_k
write(1,*)'umbr_mean=',umbr_mean
write(1,*)'steps_max=',steps_max
write(1,*)'freq_meta',freq_meta
write(1,*)'h0=',h0
write(1,*)'width=',width
write(1,*)'wtdt (a.u.)=',wtdt !deltaT parameter in energy unit (a.u.)
write(1,*)'wtdt (K)=',wtdt/kb !deltaT parameter in K is printed
write(1,*)'pfrq=',pfrq !printing frequency
if(restart)then 
  write(1,*)'restart is true'
else
  write(1,*)'restart is false'
end if
if(restart)then 
  write(1,*)'mts is true; nmts=',nmts
else
  write(1,*)'mts is false'
end if

if(restart)then
   myfilestat='old'
   myfilepos='append'
else
   myfilestat='unknown'
   myfilepos='asis'
end if

open( unit =21, file = 'TRAJ',status=trim(myfilestat),position=trim(myfilepos))
open( unit =22, file = 'ENERGIES',status=trim(myfilestat),position=trim(myfilepos))
open( unit =23, file = 'CV_VAL',status=trim(myfilestat),position=trim(myfilepos))
open( unit =24, file = 'AVG_TEMP',status=trim(myfilestat),position=trim(myfilepos))
open( unit =25, file = 'MTD',status=trim(myfilestat),position=trim(myfilepos))

meta_max=steps_max/freq_meta+1   !max mtd steps
print *, 'meta_max =', meta_max
WRITE(50,*) meta_max-1
WRITE(50,*) width
if(meta_max.le.0)stop 'error meta_max <=0'

umb_cv=1          !define cv on which US works 
meta_cv=2         !define cv on which MTD works

allocate(height(meta_max))
allocate(cv_history_s(meta_max),cv_history_z(meta_max))

print *, 'allocatiions done'

mass(1:nd)=mass(1:nd)*amu_to_au !mass of particle     
mass_s(1:ns)=mass_s(1:ns)*amu_to_au
!==================================
x_i(1,1)=1.d0
x_i(1,ndim)=-1.d0

x_i(npoints,1)=1.d0
x_i(npoints,ndim)=1.d0
!----------------------------------

step=0          !strting step md

ndof=nd*1
ndof_s=ns*1

kbt0=kb*t0
kbt0_s=kb*t0_s

ke0=0.5d0*kb*t0*dfloat(nd)
ke_s=0.5d0*kb*t0_s*dfloat(ns)

nmtd=0
!--------------------------------------
if(restart)then
  open(55,file='restart',status='old')
  read(55,*)v(1:2),v_s(1:2)
  read(55,*)x(1:2),x_s(1:2)
  read(55,*)step
  close(55)
  print *, 'finished reading restart file'
!  open(55,file='MTD',status='old',IOSTAT=ios)
  rewind(25)
!  if(ios.eq.0)then
    nmtd=0
    do
      nmtd=nmtd+1
      if(nmtd.le.meta_max)&
       read(25,*,iostat=ios2) i,j,cv_history_s(nmtd),height(nmtd),gausspot !*alpha
      if(ios2.ne.0)exit
      if(nmtd.gt.meta_max)stop 'something is wrong...nmtd>metamax'
    end do
    nmtd=nmtd-1
    print *, 'finished reading MTD file for restart: #hills=',nmtd
    close(25)
    open( unit =25, file = 'MTD',status='old',position='APPEND')
!  else
!    print *, 'cant find MTD file..skipping restart'
!    nmtd=0          !initial nmtd
!  end if
!----------------------------------------------------------------------------------------
else
!initial velocities for the physical degrees of freedom
   print *, 'initial velocity to assign 1'
   do i=1,nd
     call random_number(u1)
     call random_number(u2)
     v(i)= dsqrt(-2.d0*dlog(u1))*dcos(2.d0*pi*u2)*dsqrt(2.d0*ke0/mass(i)/dfloat(nd))
   end do
   call v_scal(nd,v,ke0,mass)
   !initial velocities for the auxiliary degrees of freedom
   print *, 'initial velocity to assign 2'
   do i=1,ns
     call random_number(u1)
     call random_number(u2) 
     v_s(i)=dsqrt(-2.d0*dlog(u1))*dcos(2.d0*pi*u2)*dsqrt(2.d0*ke_s/mass_s(i)/dfloat(ns))
   end do
   call v_scal(ncv,v_s,ke_s,mass_s)
   write(*,*) 'initial velocity done'

CALL init_traj(x_i,npoints,ndim)

CALL lambda_cal(x_i,npoints,ndim,lambda)

!do i=1,5
!write(6,*)x_i(i,1),x_i(i,2),lambda
!enddo


!open(16,file="test_com1.dat")
!do i=1,100
! do ir=1,100
!  x(1)=-1.5d0+(i-1)*0.03 
!  x(2)=-1.5d0+(ir-1)*0.03 

CALL path(x,x_i,lambda,npoints,ndim,s,z,dsdx,dsdy,dzdx,dzdy)
s_z(1)=s
s_z(2)=z

!  write(16,*)i,ir,x(:),s,z,dsdx,dsdy,dzdx,dzdy
! enddo
!enddo

!do i=1,5
!write(6,*)x_i(i,1),x_i(i,2),lambda
!enddo

!stop
!PRINT*,s,z
!setting initial position of auxiliary variables to physical variables
!   do i=1,ns
 !    do j=1,nd
  !     if (i.eq.j) x_s(i)=x(j)
!x(1)=s
!x(2)=z
!PRINT*,"sz",s,z,s_z

!checked upto here
 
x_s(1)=s_z(1) !setting initial values of auxilliary variables to physical variables
x_s(2)=s_z(2) 
  !   end do
  ! end do
end if
print *, 'initialize done '
!--------------------------------------------------------------------------
!PRINT*,"xs",x_s,s_z
!get the cv value
cvar_s=x_s(1)   !initializiing cv values
cvar_z=x_s(2)

!PRINT*,cvar_s,cvar_z
!!call cv_value(ns,meta_cv,x_s,cvar)
!!call cv_value(nd,meta_cv,x,cv_ist)
!!call cv_value(ns,umb_cv,x_s,cvar1)
!!call cv_value(nd,meta_cv,x,cv_ist1)
!!PRINT*,"CVAR",cvar,cvar1

print *, 'initial cvvalues done '
!--------------------------------------------------------------------------
!get initial gaussian position
if(nmtd.eq.0)then 
  cv_history_s(1)=cvar_s
  cv_history_z(1)=cvar_z
  height(1)=h0
  nmtd=1
end if

!PRINT*,"main",cvar_s,cvar_z
!--------------------------------------------------------------------------
print *, 'initial force'
f_s=0.d0
pe_s=0.d0

!s_z(1)=s
!s_z(2)=z

CALL force1(nd,x,x_s,pe,f,ext,f_s,pe_s,ns,x_i,npoints,ndim,lambda,s_z,k_s,meta,height,&
      width,nmtd,cv_history_s,cv_history_z,h0,wtdt,gausspot,us,umbr_mean,umbr_k,umb_cv,meta_cv)
!!call force1(nd,x,x_s,pe,f,ext,f_s,pe_s,ns,x_i,npoints,ndim,lambda,s_z,k_s,meta,height,&
 !!                         width,nmtd,cv_history_s,cv_history_z,h0,wtdt,gausspot) !,dsdx,dsdy,dzdx,dzdy)

!--------------------------------------------------------------------------
!!PRINT*,x,x_s,pe,f,ext,f_harm,pe_s,x_i,lambda,s_z,k_s
!if (ext) call force_ext(nd,ns,x,x_s,pe_s,f_s,k_s,nmtd,height,width,cvar,cvar1,cv_history,cv_history1, &
!                        meta,f_harm,h0,wtdt,gausspot,us,umbr_k,umbr_mean,x_i,npoints,ndim,s_z)
!call force(nd,x,pe,f,ncv,ext,f_harm,pe_s,ns,x_i,npoints,ndim)

print *, 'starting MD'
md_loop : do
step= step + 1
ke=0.d0
te=0.d0
pe=0.d0

if(lang)then

call lang_v_update1_new3(nd,dt,mass,kbt0,x,v,f,rnd,glang)
call lang_x_update_new3(nd,dt,mass,kbt0,x,v,f,rnd,glang)

   if(ext)then
      call lang_v_update1_new3(ns,dt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
      call lang_x_update_new3 (ns,dt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
  end if
else
   call verlet_x_update(nd,dt,mass,x,v,f)
   if(ext)call verlet_x_update(ns,dt,mass_s,x_s,v_s,f_s)
   call verlet_v_update(nd,dt,mass,x,v,f)
   if(ext)call verlet_v_update(ns,dt,mass_s,x_s,v_s,f_s)
end if
!--------------------------------------------------------------------------------
!CALL path(x_s,x_i,lambda,npoints,ndim,s,z,dsdx,dsdy,dzdx,dzdy)

cvar_s=x_s(1)
cvar_z=x_s(2)
!PRINT*,"sa",x_s,s_z
!!call cv_value(ns,meta_cv,x_s,cvar)
!!call cv_value(nd,meta_cv,x,cv_ist)
!!call cv_value(ns,umb_cv,x_s,cvar1)
!!call cv_value(nd,umb_cv,x,cv_ist1)

if(meta)then
if(mod(step,freq_meta).eq.0) then
  nmtd=nmtd+1
  
  cv_history_s(nmtd)=cvar_s  !for s
  cv_history_z(nmtd)=cvar_z  !for z
  height(nmtd)= h0*dexp(-gausspot/wtdt)

!PRINT*,nmtd,cv_history_s(nmtd),cvar_s ,cv_history_s(nmtd),cvar_s
!stop
!PRINT*,height(nmtd),gausspot,wtdt
  !print *, "updating hills at ", step, " Hill #=", nmtd, " Height =", height(nmtd)
!!  write(25,'(2I16,4F16.6)') nmtd,step,cv_history_s(nmtd),cv_history_z(nmtd),height(nmtd),gausspot !*alpha
write(25,'(2I16,4F16.6)') nmtd,step,cv_history_z(nmtd),height(nmtd),gausspot
end if
end if
!---------------------------------------------------------------------------------------


!if (ext)call force_ext(nd,ns,x,x_s,pe_s,f_s,k_s,nmtd,height,width,cvar,cvar1,&
!            cv_history,cv_history1,meta,f_harm,h0,wtdt,gausspot,us,umbr_k,umbr_mean,x_i,npoints,ndim,s_z)
!call force(nd,x,pe,f,ncv,ext,f_harm,pe_s,ns,x_i,npoints,ndim)

CALL force1(nd,x,x_s,pe,f,ext,f_s,pe_s,ns,x_i,npoints,ndim,lambda,s_z,k_s,meta,height,&
    width,nmtd,cv_history_s,cv_history_z,h0,wtdt,gausspot,us,umbr_mean,umbr_k,umb_cv,meta_cv)
!!call force1(nd,x,x_s,pe,f,ext,f_s,pe_s,ns,x_i,npoints,ndim,lambda,s_z,k_s,meta,height,&
! !                         width,nmtd,cv_history_s,cv_history_z,h0,wtdt,gausspot) !,dsdx,dsdy,dzdx,dzdy)


if(lang)then
  call lang_v_update2_new3(nd,dt,mass,kbt0,x,v,f,rnd,glang)
  if(ext)call lang_v_update2_new3(ns,dt,mass_s,kbt0_s,x_s,v_s,f_s,rnd_s,glang_s)
else 
   call verlet_v_update(nd,dt,mass,x,v,f)
   if(ext)call verlet_v_update(ns,dt,mass_s,x_s,v_s,f_s)
end if

!calculate total energy and temperature
ke=0.0d0
do i=1,nd 
   ke=ke+0.5d0*mass(i)*v(i)*v(i)
end do
t_inst=2.d0*ke/(dfloat(nd)*kb)


ke_s=0.0d0
if(ext) then
  do i=1,ns 
    ke_s=ke_s+0.5d0*mass_s(i)*v_s(i)*v_s(i)
  end do
  t_s_inst=2.d0*ke_s/(dfloat(ns)*kb)
  ke=ke+ke_s
end if

te=ke+pe !pe has pe_s included (inside force routine)

avg_t=avg_t+t_inst
avg_t_s=avg_t_s+t_s_inst

!================================================

if(mod(step,pfrq).eq.0) then
  write(1,'(I16,6F16.6)')step,t_inst,t_s_inst,pe,ke,te,gausspot
  write(24,*)step, avg_t/real(step), avg_t_s/real(step)
  write(22,'(I16,5F16.6)')step,t_inst,t_s_inst,pe,ke,te
  write(21,'(I16,5F16.8)')step,x(1:nd),v(1:nd),pe!,v(1:nd)
  write(23,'(I16,4F16.4)')step,x_s(1:ns),s_z(1:ns)!,v_s(1:ns)
  open(55,file='restart')
   write(55,'(4e16.8)')v(1:2),v_s(1:2)
   write(55,'(4e16.8)')x(1:2),x_s(1:2)
   write(55,*)step
  close(55)
  !!print *, 'wrote retart file'
end if

!if(step.eq.2)stop 'stopping'

if (step.ge.steps_max) exit md_loop
end do md_loop
!---------------------------------------------------
  open(55,file='restart',status='old')
   write(55,'(4e16.8)')v(1:2),v_s(1:2)
   write(55,'(4e16.8)')x(1:2),x_s(1:2)
   write(55,*)step
  close(55)
  print *, 'wrote retart file'

end program mdcode

!===============================================================================
SUBROUTINE v_scal(nd,v,ke_init,mass)
IMPLICIT NONE

integer :: nd,i
real*8 :: v(nd),ke_init,mass(nd),ke,scal

ke=0.0d0
do i=1,nd
 ke=ke+0.5d0*mass(i)*v(i)**2
end do
scal=dsqrt(ke_init/ke)
DO i=1,nd
 v(i)=v(i)*scal
END DO

END SUBROUTINE 
!=====================================================================
SUBROUTINE force_ext(nd,ns,x,x_s,pe_s,f_s,k_s,nmtd,height,width,cvar,cvar1,&
    cv_history,cv_history1,meta,f_harm,h0,wtdt,gausspot,us,umbr_k,umbr_mean,x_i,npoints,ndim,s_z)
IMPLICIT NONE

integer :: nd,ns,i,j,it
real*8  :: x(nd),f(nd),pe_s,k_s,f_s(ns),x_s(ns), gausspot, umbr_k, umbr_mean,x_i(npoints,ndim)

real*8  :: expn,cv_history(*),cv_history1(*)
integer :: meta_cv, nmtd,npoints,ndim

real*8::cvar,cvar1,diff,diff1,width,height(*)
logical :: meta,us

real*8  :: f_harm(ns),h0,wtdt,step_gauss,s,z,s_z(*)
integer :: umb_cv

!REAL*8::lambda,dsdx,dsdy,dzdx,dzdy


pe_s=0.d0
f_s=0.d0

do i=1,ns
  !do j=1,nd
    ! if (i.eq.j) then
       f_s(i)=f_s(i)-k_s*(s_z(i)-x_s(i)) !v=-kx as f=0.5*k*x**2  !x_s=>s,x=>S
       pe_s=pe_s+0.5d0*k_s*(s_z(i)-x_s(i))**2
    !end if
   !end do
end do

do i=1,ns
  f_harm(i)=f_s(i)
end do

!force contribution from umbrella.
!if(us) then
  umb_cv=1
  meta_cv=2
  !do i=1,nd
!     if (i.eq.umb_cv) then
  !!     f_s(i)=f_s(i)-umbr_k*(x_s(i)-umbr_mean) !v=-kx as f=0.5*k*x**2
  !!   pe_s=pe_s+0.5d0*umbr_k*(x_s(i)-umbr_mean)**2
cvar1=x_s(umb_cv)  
cvar=x_s(meta_cv)
 !!  end if
  !!end do
!!end if

!force contribution from bias.
!meta=.false.
!!meta_cv=2
gausspot=0.d0
!!if(meta) then
  do it= 1,nmtd-1
     diff=(cvar-cv_history(it)) !for z ==>y
     diff1=cvar1-cv_history1(it)  !for s ==>x
     step_gauss=height(it)*dexp(-0.5d0*((diff/width)**2+ (diff1/width)**2))
      gausspot= gausspot+step_gauss
!     do i=1,nd
!!       if (i.eq.meta_cv) then
!!         f_s(i)=f_s(i)+(diff/(width)**2.0d0)*step_gauss
  f_s(1)=f_s(1)+(diff1/(width)**2.0d0)*step_gauss
  f_s(2)=f_s(2)+(diff/(width)**2.0d0)*step_gauss
!!     end if 
 !    end do
  end do
  pe_s=pe_s+gausspot
!TODO
!  height(nmtd)= h0*dexp(-gausspot/wtdt)
!!end if

end subroutine    
!==========================================================================================
SUBROUTINE force(nd,x,pe,f,ncv,ext,f_harm,pe_s,ns,x_i,npoints,ndim)
IMPLICIT NONE

INTEGER :: nd,ncv,ns,i,j
REAL*8  :: x(nd),f(nd),pe,f_s(ns),pe_s,f_harm(*)
LOGICAL :: us, ext
REAL*8,PARAMETER:: v0=0.002d0, a=1.d0, v1=0.1d0 !in atomic units
REAL*8::dsdx,dsdy,dzdx,dzdy,s,z,lambda
INTEGER::npoints,ndim
REAl*8::x_i(npoints,ndim)
f(1) = -v0*4.0d0*x(1)*(x(1)**2-a**2)/a**4 
f(2) = -v0*4.0d0*x(2)*(x(2)**2-a**2)/a**4 
pe=v0*(x(1)**2-a**2)**2 + v0*(x(2)**2-a**2)**2

!force contribution due to extended part
IF(ext)THEN
!  do i=1,nd
 !    do j=1,ns
  !     if (i.eq.j) then
   !      f(i)=f(i)-f_harm(j) !v=-kx as f=0.5*k*x**2
!CALL lambda_cal(x_i,npoints,ndim,lambda)
!Call path(x,x_i,lambda,npoints,ndim,s,z,dsdx,dsdy,dzdx,dzdy) 
f(1)=f(1)-f_harm(1)*dsdx-f_harm(2)*dzdx
f(2)=f(2)-f_harm(1)*dsdy-f_harm(2)*dzdy  
!PRINT*,"f",f
 !    end if
   !  end do
  !end do
  pe=pe+pe_s
end if

END SUBROUTINE force
!===========================================================================
SUBROUTINE cv_value(nd,meta_cv,x,cvar)
IMPLICIT NONE

INTEGER:: i,icv,nd,meta_cv
REAL*8::cvar,x(*)

DO i=1,nd
  IF (i.EQ.meta_cv) cvar = x(i)
END DO
RETURN

END SUBROUTINE cv_value
!============================================================================================
SUBROUTINE verlet_v_update(nd,dt,mass,x,v,f)
IMPLICIT NONE

INTEGER :: nd,i
REAL*8 :: dt, mass(nd), x(nd), v(nd), f(nd)

DO i=1,nd
  v(i)=v(i)+dt*0.5d0*f(i)/mass(i)
END DO

END SUBROUTINE verlet_v_update
!====================================================================================
SUBROUTINE verlet_x_update(nd,dt,mass,x,v,f)
IMPLICIT NONE

INTEGER :: nd,i
REAL*8 :: dt, mass(nd), x(nd), v(nd), f(nd)

DO i=1,nd
  x(i)=x(i)+dt*v(i)+0.5d0*dt*dt*f(i)/mass(i)
END DO

END SUBROUTINE verlet_x_update
!====================================================================================
SUBROUTINE lang_x_update_new3(nd,dt,mass,kbt,x,v,f,rnd,gamma)
IMPLICIT NONE

INTEGER:: nd,i
REAL*8:: x(nd),v(nd),f(nd),rnd(nd),mass(nd), dt, gamma, kbt
REAL*8,PARAMETER :: pi=4.d0*atan(1.d0)

DO i=1,nd
   x(i)=x(i)+dt*v(i)
END DO

END SUBROUTINE lang_x_update_new3
!=============================================================================
SUBROUTINE lang_v_update1_new3(nd,dt,mass,kbt,x,v,f,rnd,gamma)
IMPLICIT NONE

INTEGER :: nd
REAL*8    :: x(nd),v(nd),f(nd),rnd(nd),mass(nd), dt, gamma, kbt
REAL*8    :: FACT,z1,z2,u1,u2,fact2
INTEGER   :: i
REAL*8,PARAMETER:: pi=4.d0*atan(1.d0)

DO i=1,nd
   call random_number(u1)
   call random_number(u2) 
   z1=dsqrt(-2.d0*dlog(u1))*dcos(2.d0*pi*u2) 
   z2=dsqrt(-2.d0*dlog(u1))*dsin(2.d0*pi*u2) 
   rnd(i)=z1
   fact=dt/(2.d0*mass(i))
   v(i)=v(i)+fact*f(i)
   v(i)=v(i)-mass(i)*v(i)*gamma &
       *fact+dsqrt(fact*2.d0*mass(i))*0.5d0*dsqrt(2.d0 &
       *kbt/mass(i)*gamma)*z1
END DO

END SUBROUTINE lang_v_update1_new3
!=================================================================================
SUBROUTINE lang_v_update2_new3(nd,dt,mass,kbt,x,v,f,rnd,gamma)
IMPLICIT NONE

INTEGER::nd,i
REAL*8::x(nd),v(nd),f(nd),rnd(nd),mass(nd), dt, gamma, kbt
REAL*8:: FACT,z1,z2,u1,u2,fact2
REAL*8,PARAMETER::pi=4.d0*atan(1.d0)

DO i=1,nd
   fact=dt/(2.d0*mass(i))
   v(i)=v(i)+fact*f(i)
   v(i)=v(i)-mass(i)*v(i)*gamma &
       *fact+dsqrt(fact*2.d0*mass(i))*0.5d0*dsqrt(2.d0 &
       *kbt/mass(i)*gamma)*rnd(i)
END DO

END SUBROUTINE lang_v_update2_new3
!==================================================================================
SUBROUTINE  path(x,x_i,lambda,npoints,ndim,s,z,dsdx,dsdy,dzdx,dzdy)
IMPLICIT NONE

INTEGER:: npoints,ndim,i

REAL*8::x_i(npoints,ndim),x(ndim),z,s,num_s,den,num_dsdx,num_dsdy,num_dzdx,num_dzdy,dsdx, &
       dsdy,dzdx,dzdy,lambda,l

l=lambda
num_s=0.d0
den=0.d0
num_dsdx=0.d0
num_dsdy=0.d0
num_dzdx=0.d0
num_dzdy=0.d0


!CALL init_traj(x_i,npoints,ndim)

!CALL lambda_cal(x_i,npoints,ndim,l)

DO i=1,npoints
 num_s=num_s+(i-1)*dexp(-l*((x(1)-x_i(i,1))**2+(x(2)-x_i(i,2))**2))
 den=den+dexp(-l*((x(1)-x_i(i,1))**2+(x(2)-x_i(i,2))**2))
END DO

s=num_s/den/(npoints-1)
z=-dlog(abs(den))/l

DO i=1,npoints
 num_dsdx=num_dsdx+2.d0*(x(1)-x_i(i,1))*dexp(-l*((x(1)-x_i(i,1))**2+(x(2)-x_i(i,2))**2))* &
                   (-s*(npoints-1)+(i-1))
 num_dsdy=num_dsdy+2.d0*(x(2)-x_i(i,2))*dexp(-l*((x(1)-x_i(i,1))**2+(x(2)-x_i(i,2))**2))* &
                   (-s*(npoints-1)+(i-1))

 num_dzdx=num_dzdx-l*dexp(-l*((x(1)-x_i(i,1))**2+(x(2)-x_i(i,2))**2))*2.d0*(x(1)-x_i(i,1))
 num_dzdy=num_dzdy-l*dexp(-l*((x(1)-x_i(i,1))**2+(x(2)-x_i(i,2))**2))*2.d0*(x(2)-x_i(i,2))
END DO

dsdx=l*num_dsdx/den/(1-npoints)
dsdy=l*num_dsdy/den/(1-npoints)
dzdx=-num_dzdx/den/l
dzdy=-num_dzdy/den/l

END SUBROUTINE path

!============================================================================
SUBROUTINE init_traj(x_i,npoints,ndim)
IMPLICIT NONE

INTEGER::j,npoints,ndim

REAL*8::x_i(npoints,ndim)

!x_i(2,1)=0.05d0
!x_i(2,2)=-1.05d0

!x_i(3,1)=-1.05d0
!x_i(3,2)=-1.05d0

!x_i(4,1)=-1.05d0
!x_i(4,2)=0.05d0
DO j=2,npoints-1
 x_i(j,1:ndim)=x_i(1,1:ndim)+(j-1)*(x_i(npoints,1:ndim)-x_i(1,1:ndim))/(npoints-1)
!PRINT*,x_i(j,1:ndim)
END DO
!stop
!do j=1,5
!write(6,*)x_i(j,1),x_i(j,2)
!enddo
!stop
END SUBROUTINE init_traj
!=============================================================================
SUBROUTINE lambda_cal(x_i,npoints,ndim,l)
IMPLICIT NONE

REAL*8::l,d,x_i(npoints,ndim),dx,dy
INTEGER::i,npoints,ndim

d=0.d0
DO i=1,npoints-1
dx=(x_i(i,1)-x_i(i+1,1))
dy=(x_i(i,2)-x_i(i+1,2))

d=d+dsqrt(dx**2+dy**2)
!write(6,*)i,dsqrt(dx**2+dy**2)
END DO

l=2.3d0*(npoints-1)/d
!PRINT*,l,d
!stop
END SUBROUTINE 
!============================================================================
SUBROUTINE force1(nd,x,x_s,pe,f,ext,f_s,pe_s,ns,x_i,npoints,ndim,lambda,s_z,k_s,meta,height,&
       width,nmtd,cv_history_s,cv_history_z,h0,wtdt,gausspot,us,umbr_mean,umbr_k,umb_cv,meta_cv) !,dsdx,dsdy,dzdx,dzdy)
IMPLICIT NONE

INTEGER :: nd,ns,i,j,npoints,ndim

INTEGER::nmtd,it !Metadynamics Variables
INTEGER::umb_cv,meta_cv

REAL*8  :: x(nd),f(nd),pe,f_s(ns),pe_s,f_harm(nd),k_s,x_s(ns),s_z(ndim)

LOGICAL ::  ext,meta,us

REAL*8,PARAMETER:: v0=0.002d0, a=1.d0, v1=0.1d0 !in atomic units

REAL*8::dsdx,dsdy,dzdx,dzdy,s,z,lambda,x_i(npoints,ndim)

REAL*8::height(*),width,cv_history_s(*),cv_history_z(*),gausspot,diff_s,diff_z,cvar_s, &
        wtdt, h0,cvar_z,step_gauss  !variables for metadynamics


REAl*8::umbr_mean,umbr_k


f(1) = -v0*4.0d0*x(1)*(x(1)**2-a**2)/a**4
f(2) = -v0*4.0d0*x(2)*(x(2)**2-a**2)/a**4
pe=v0*(x(1)**2-a**2)**2 + v0*(x(2)**2-a**2)**2
pe_s=0.d0
f_s=0.d0

!force contribution due to extended part
IF(ext)THEN
!PRINT*,x_i,lambda
!stop
Call path(x,x_i,lambda,npoints,ndim,s_z(1),s_z(2),dsdx,dsdy,dzdx,dzdy) 
!PRINT*,x,s_z(1),s_z(2)
f_harm(1)=k_s*(s_z(1)-x_s(1))
!PRINT*,s_z,x_s
f_harm(2)=k_s*(s_z(2)-x_s(2))


f(1)=f(1)-f_harm(1)*dsdx-f_harm(2)*dzdx
f(2)=f(2)-f_harm(1)*dsdy-f_harm(2)*dzdy

!force due to extended part
do i=1,ns
       f_s(i)=f_s(i)+f_harm(i) !v=-kx as f=0.5*k*x**2  !x_s(1)=>s,x_s(2)=>z
       pe_s=pe_s+(0.5d0*f_harm(i)**2)/k_s
end do

cvar_z=x_s(2)
gausspot=0.d0

!PRINT*,"umb_cv",umb_cv
!PRINT*,"meta_cv",meta_cv
!umb_cv=1
!meta_cv=2

!PRINT*,us,meta

!force contribution from umbrella
!US along s
if(us)then
f_s(umb_cv)=f_s(umb_cv)-umbr_k*(x_s(umb_cv)-umbr_mean)
pe_s=pe_s+0.5d0*umbr_k*(x_s(umb_cv)-umbr_mean)**2
end if

!force contribution from MTD
!MTD along z

if(meta)then
do it=1,nmtd-1
  diff_z=cvar_z-cv_history_z(it)
  step_gauss=height(it)*dexp(-0.5d0*(diff_z/width)**2)
  gausspot=gausspot+step_gauss
  f_s(meta_cv)=f_s(meta_cv)+(diff_z/width**2)*step_gauss
end do
  !f_s(meta_cv)=f_s(meta_cv)
end if


!PRINT*,h0
!!cvar_s=x_s(1)
!!cvar_z=x_s(2)
!!gausspot=0.d0

!!if(meta) then
!! do it= 1,nmtd-1
   
  !!  diff_s=(cvar_s-cv_history_s(it)) !for s ==>x
   !! diff_z=(cvar_z-cv_history_z(it))  !for z ==>y 
   !! step_gauss=height(it)*dexp(-0.5d0*(diff_s*diff_s+diff_z*diff_z)/width**2)
   !! gausspot= gausspot+step_gauss
    !!f_s(1)=f_s(1)+(diff_s/(width)**2)*step_gauss
    !!f_s(2)=f_s(2)+(diff_z/(width)**2)*step_gauss
!! end do

!!end if

pe=pe+pe_s+gausspot

END IF

END SUBROUTINE force1
!===============================================================================


