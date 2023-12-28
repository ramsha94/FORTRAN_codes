PROGRAM histogram
IMPLICIT NONE

REAL*8:: grid_min,grid_max,grid_width,x,y,intx
REAL*8,ALLOCATABLE:: Prob(:),marks(:)

INTEGER::isteps,nsteps,igrid,ngrids,i
grid_min=0.d0
grid_max=100.d0
grid_width=10.d0
nsteps=10

OPEN(1,file="trajectory.xy")
!READ(1,*) x,y

!CALL read_trajectory( R )
!PRINT*,"Enter the starting,ending and width of grid respetively" 
!READ *, grid_min, grid_max, grid_width
!PRINT*,"Staring value of grid:",grid_min
!PRINT*,"Ending value of grid:",grid_max
!PRINT*,"Width of grid:",grid_width

ngrids=nint((grid_max-grid_min)/grid_width + 0.5) 
!Print*,(grid_max-grid_min)/grid_width
PRINT*,"ngride=",ngrids
ALLOCATE(Prob(ngrids))
ALLOCATE(marks(nsteps))

DO i=1,nsteps
OPEN(1,file="trajectory.xy")
READ(1,*) x,y
marks(i)=y
END DO
!Print*,"a",x,y,i
DO i=1,nsteps
 DO isteps=1,nsteps
! READ(1,*) x,y
 !PRINT*,"x ",x,y
 !PRINT*,"i",isteps,grid_width*(isteps-1) ,grid_width*isteps
y=marks(isteps) 
PRINT*,"y",y,(grid_width*(i-1)),grid_width*(i)  
IF((y.ge.(grid_width*(i-1))).AND.(y.lt.(grid_width*(i))))THEN
     PRINT*,"IF",y
     igrid=int((y-grid_min)/grid_width) + 1
     Prob(igrid)=Prob(igrid)+1.d0
!     PRINT*,"s",Prob(igrid)
   END IF
 END DO
END DO
CLOSE(1)
!END DO

intx=0.d0

DO igrid=1,ngrids
 intx=intx+Prob(igrid)*grid_width
END DO

Prob(1:ngrids)=Prob(1:ngrids)/ngrids

DO igrid=1,ngrids
PRINT*,igrid,Prob(igrid)
END DO

DEALLOCATE(Prob)
DEALLOCATE(marks)
END PROGRAM
