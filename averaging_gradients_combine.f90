!For combining gradients from different files and Averaging the Gradients
PROGRAM test
IMPLICIT NONE

REAL*8,ALLOCATABLE :: prob1(:,:,:),prob2(:,:,:),prob_comb1(:,:),prob_comb2(:,:),sum1(:,:),x(:,:,:),y(:,:,:)
REAL*8,ALLOCATABLE :: gridmin(:),gridmax(:),griddiff(:),sum2(:,:),fes(:,:,:),fes_comb(:,:) 
REAL*8,ALLOCATABLE :: lower_wall(:),upper_wall(:)

REAL*8:: dum

CHARACTER(LEN=50), ALLOCATABLE :: filename(:)

INTEGER:: i,j,k,ir,nr,ndim
INTEGER::index1,index2

INTEGER,ALLOCATABLE::nbin(:)

OPEN(44,FILE='input',STATUS='old')
READ(44,*) nr,ndim
ALLOCATE(gridmin(ndim),gridmax(ndim),griddiff(ndim))
ALLOCATE(nbin(ndim))
ALLOCATE(lower_wall(nr),upper_wall(nr))
ALLOCATE(filename(nr))
READ(44,*) gridmin(1:ndim)
READ(44,*) gridmax(1:ndim)
READ(44,*) nbin(1:ndim)

DO ir=1,nr
READ(44,*) lower_wall(ir),upper_wall(ir)
   READ(44,'(a)')filename(ir)
END DO
!nbin(1:ndim) = NINT((gridmax(1:ndim)-gridmin(1:ndim))/griddiff(1:ndim))+1
griddiff(1:ndim) = (gridmax(1:ndim)-gridmin(1:ndim))/(nbin(1:ndim)-1)
PRINT*,"GRIDMIN ",gridmin(1:ndim)
PRINT*,"GRIDMAX ",gridmax(1:ndim)
PRINT*,"GRIDDIF ",griddiff(1:ndim)
PRINT*,"NBIN ",nbin(1:ndim)  !checked
!stop
!--------------Reading the probability from different files-----------------------------------------
ALLOCATE(x(nbin(1),nbin(2),nr))
ALLOCATE(y(nbin(1),nbin(2),nr))
ALLOCATE(prob1(nbin(1),nbin(2),nr),prob2(nbin(1),nbin(2),nr))
ALLOCATE(fes(nbin(1),nbin(2),nr))
!

! OPEN(111,FILE="test-read.dat")
DO ir=1,nr
 OPEN(11,FILE=filename(ir),STATUS='old')
 PRINT *, '...reading PROBABILITY (',ir,')'
 DO j=1,nbin(1)
    DO k=1,nbin(2)
       READ(11,*) x(j,k,ir),y(j,k,ir),fes(j,k,ir),prob1(j,k,ir),prob2(j,k,ir)  !Checked
!WRITE(111,*) x(j,k,ir),y(j,k,ir),prob(j,k,ir)   !This is okay
 END DO
END DO
END DO
!CLOSE(111)
!stop
!----Storing probabilty in a single array and writing in file---------------------
!
ALLOCATE(prob_comb1(nbin(1),nbin(2)))
ALLOCATE(prob_comb2(nbin(1),nbin(2)))
ALLOCATE(fes_comb(nbin(1),nbin(2)))
ALLOCATE(sum1(nbin(1),nbin(2)))  
ALLOCATE(sum2(nbin(1),nbin(2)))  

sum1=1.d0   !modified
sum2=1.d0
!PRINT*,"After here important i,j,x,prob"
!stop

prob_comb1=0.d0
prob_comb2=0.d0
DO ir=1,nr
  DO j=1,nbin(1)
   DO k=1,nbin(2)
   index1=NINT((x(j,k,ir)-gridmin(1))/griddiff(1)) + 1
   index2=NINT((y(j,k,ir)-gridmin(2))/griddiff(2)) + 1
!PRINT*,index1,index2,x(j,k,ir),y(j,k,ir)
!PRINT*,x(j,k,ir),y(j,k,ir),index1,index2
!DO dim=1,ndim
 ! IF((x(j,k,dim).LE.upper_wall(dim)).AND.(x(j,k,dim).GE.upper_wall(dim)))THEN
 !       prob_comb1(index1,index2)=prob_comb1(index1,index2) + prob1(j,k,1)
!END DO
    IF((x(j,k,ir).GE.lower_wall(ir)).AND.(x(j,k,ir).LE.upper_wall(ir)))THEN   !TODO for more Walls properly
       ! prob_comb1(index1,index2)=prob_comb1(index1,index2) + prob1(j,k,ir)  
       ! prob_comb2(index1,index2)=prob_comb2(index1,index2) + prob2(j,k,ir)
       
      IF(ir.GE.2.AND.(x(j,k,ir).LE.upper_wall(ir-1).AND.(x(j,k,ir).GE.lower_wall(ir))))THEN  !condition for averaging
       PRINT*,x(j,k,ir),prob1(j,k,ir-1),prob1(j,k,ir)  !Checked 
         prob_comb1(index1,index2) =prob_comb1(index1,index2) + ( prob1(j,k,ir-1) + prob1(j,k,ir))
         prob_comb2(index1,index2) =prob_comb2(index1,index2) + ( prob2(j,k,ir-1) + prob1(j,k,ir))
         prob_comb1(index1,index2)=prob_comb1(index1,index2)/2.d0
         prob_comb2(index1,index2)=prob_comb2(index1,index2)/2.d0
       END IF
        fes_comb(index1,index2)=fes(j,k,ir)
     !   IF (prob1(j,k,ir) .NE. 0.d0) sum1(index1,index2)=sum1(index1,index2) + 1
    END IF
!    IF(x(j,k,2).GT.0.15d0)THEN    !TODO
 !       prob_comb1(index1,index2)=prob_comb1(index1,index2) + prob1(j,k,2)
      ! IF (prob2(j,k,ir) .NE. 0.d0) sum2(index1,index2)=sum2(index1,index2) + 1
  !  END IF
   !prob_comb2(index1,index2)=prob_comb2(index1,index2) + prob2(j,k,ir)
!   IF (prob1(j,k,ir) .NE. 0.d0) sum1(index1,index2)=sum1(index1,index2) + 1
   END DO
  END DO
END DO
!
OPEN(77,file="mean_force.in")
DO i=1,nbin(1)

  DO j=1,nbin(2)
    WRITE(77,*) gridmin(1)+(i-1)*griddiff(1),gridmin(2)+(j-1)*griddiff(2), &
            prob_comb1(i,j),prob_comb2(i,j),fes_comb(i,j)
  END DO
!WRITE(77,*)
END DO
!
CLOSE(11)
CLOSE(44)
CLOSE(77)

DEALLOCATE(gridmin,gridmax,griddiff)
DEALLOCATE(filename,nbin,x,prob1,prob2)

END PROGRAM
