!need to change integrate_ggmt
PROGRAM mdcode
IMPLICIT NONE

REAL*8 :: ke0,glang,kb
REAL*8  :: dt,te,ke,pe,t0,t_inst, avg_t
REAL*8 :: x,v,mass,f

INTEGER :: i,ndof,steps_max,step, pfrq

LOGICAL :: ggmt

REAL*8, PARAMETER :: kbt0=1.d0,k=1.d0!3.16e-6  !Boltman constant in a.u. K^-1

REAL*8::cons
CHARACTER(LEN=20) :: myfilestat, myfilepos
!=======================================
!GGMT variable declaration
REAL*8::Q1,Q2,eta1,eta2,v1,v2
REAL*8,DIMENSION(4)::dt_ggmt
REAL*8,PARAMETER::tau=1.d0
INTEGER,PARAMETER::n_respa_ggmt=1
!=========================================
OPEN( unit =20, file = 'input',status='unknown')
!=====================================================
avg_t=0.0d0
!==================================
!GGMT variables allocation
!================================
READ(20,*)x
READ(20,*)mass
READ(20,*)dt
READ(20,*)t0
READ(20,*)steps_max
READ(20,*)pfrq

kb=1.d0
!================================================= =
!GGMT variables initialization
Q1= 1.d0*kbt0*tau * tau
Q2= (8.d0/3.d0) *(kbt0**3.d0)*tau *tau 

dt_ggmt(1)  = dt/n_respa_ggmt
dt_ggmt(2) = (1.d0 / (2.d0 - 2.d0**(1.d0/3.d0))) * dt_ggmt(1)
dt_ggmt(3)  = dt_ggmt(1) -  2.d0*dt_ggmt(2)
dt_ggmt(4)  = dt_ggmt(2)

eta1= 0.d0
eta2= 0.d0

v1=1.d0/Q1!-dsqrt( kb/Q1)
v2=-1.d0/Q2!-dsqrt( kb/Q2)
!=============================================
   myfilestat='unknown'
   myfilepos='asis'
open( unit =21, file = 'TRAJ',status=trim(myfilestat),position=trim(myfilepos))
open( unit =22, file = 'ENERGIES',status=trim(myfilestat),position=trim(myfilepos))
open( unit =24, file = 'AVG_TEMP',status=trim(myfilestat),position=trim(myfilepos))
OPEN(11,file="ggmt.dat")
!==================================
step=0          !strting step md
ndof=1
ke0=0.5d0*kbt0

v=1.0/mass

f=-k*x   !force calculation
print *, 'starting MD'
md_loop : do
step= step + 1
ke=0.d0
te=0.d0
pe=0.d0
!=================================================================================================
CALL integrate_ggmt(mass,kbt0,n_respa_ggmt,dt_ggmt,v,v1,v2,eta1,eta2,Q1,Q2) !call properly here

v=v+0.5d0*dt*f/mass
x=x+dt*v

f=-k*x   !force

v=v+0.5*dt*f/mass
CALL integrate_ggmt(mass,kbt0,n_respa_ggmt,dt_ggmt,v,v1,v2,eta1,eta2,Q1,Q2) !call properly here
!==================================================================================================
ke=ke+0.5d0*mass*v**2.d0
pe=0.5d0*k*x**2.d0
t_inst=2.d0*ke/kb
te=ke+pe !pe has pe_s included (inside force routine)
cons=te+ (v1*Q1)**2/(2.d0*Q1) +eta1 + (v2*Q2)**2/(2.d0*Q2) +  eta2
avg_t=avg_t+t_inst
!================================================
if(mod(step,pfrq).eq.0) then
  write(24,*)step, avg_t/real(step)
  write(22,'(I16,5F16.6)')step,t_inst,pe,ke,te,cons
  write(21,'(I16,5F16.8)')step,x,v,pe !,v(1:nd)
  write(11,'(I16,8F16.4)')step,v1,v2,eta1,eta2
end if

if (step.ge.steps_max) exit md_loop
end do md_loop

end program mdcode
!===============================================================================
!===========================================================================
SUBROUTINE integrate_ggmt(mass,kt,n_respa_ggmt,dt,vfict,v1,v2,eta1,eta2,Q1,Q2) !dt is array one here
IMPLICIT NONE

INTEGER iii, jjj,n_respa_ggmt
REAL*8:: aa, bb,d,kt
REAL*8:: G1,G2,Q1,Q2
REAL*8:: dt2, dt4, dt8

REAL*8::vfict,v1,v2,eta2,eta1
REAL*8::mass
REAL*8,dimension(4)::dt

d=mass*(vfict**2.d0)   !p**2/2m for ith DOF

DO iii=1,n_respa_ggmt
  DO jjj=2,4

       dt2 = 0.5*dt(jjj)  ! half time step
       dt4 = 0.25*dt(jjj)  ! quarter time step
       dt8 = 0.125*dt(jjj) ! 1/8 time step
!!PRINT*,dt2,dt4,dt8
!stop
      G1   = (d - kt)/Q1
      G2   = ((d**2)/3.d0 - (kt**2))/Q2
      v1   =v1 + dt4 * G1 
      v2  =v2 + dt4 * G2 

      aa = dexp(-dt8 * (v1 + kt*v2))
      vfict   = vfict*aa   
      d=(vfict**2)*mass

      bb         = d*(v2/3.d0)
      vfict   = vfict*dsqrt(1.d0/(1.d0 + dt2*bb))  

      vfict  = vfict*aa !change to vfict
      d=(vfict**2)*mass

      eta1= eta1 + dt2 * v1
      eta2 =eta2 + dt2 * v2*(d + kt)

      aa         = dexp(-dt8 * (v1 + kt*v2))
      vfict   = vfict*aa 
      d=(vfict**2)*mass
      
      bb         = d*(v2/3.d0)
      vfict   = vfict*dsqrt(1.d0/(1.d0 + dt2 * bb))  

      vfict  = vfict*aa !change to vfict
      d=(vfict**2)*mass

      G1   = (d - kt)/Q1
      G2   = ((d**2)/3.d0 - (kt**2))/Q2
      v1 =v1 + dt4 * G1 
      v2 =v2 + dt4 * G2 
!PRINT*,v1,v2,eta1,eta2
!stop 
   END DO
END DO
END SUBROUTINE
