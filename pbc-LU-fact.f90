!Solving Ax=B matrix using LU factorisation technique for periodic system
PROGRAM LU
IMPLICIT NONE

REAL*8,ALLOCATABLE:: x(:,:),fes(:),A(:,:)
REAL*8,ALLOCATABLE::A_temp(:,:),fes_temp(:),Y(:)
INTEGER:: i,j,k
INTEGER::M,N
REAL*8:: gridmin1,gridmax1,griddif1
REAL*8:: gridmin2,gridmax2,griddif2
REAL*8::dum,width_x,width_y ! Read this
INTEGER,ALLOCATABLE:: IPIV(:)
INTEGER:: INFO,i_temp

REAL*8,ALLOCATABLE::sum(:),sum_g(:)
REAL*8:: griddif,gridmin,gridmax
REAL*8:: diff_1, diff_2

REAL*8:: x_i,x_j,x_k
INTEGER:: nbin

gridmin1=-3.14d0
gridmax1=3.14d0
griddif1= 0.33
gridmin2=-3.14d0
gridmax2=3.14d0
griddif2=0.33d0

M = NINT((gridmax1-gridmin1)/griddif1)+1  !Read later 
N = NINT((gridmax2-gridmin2)/griddif2)+1  

PRINT*,M,N,M*N
width_x=0.33d0     
width_y=0.33d0    

OPEN(11,file="free_energy",status="unknown")
OPEN(13,file="output-coeff",status="unknown")
OPEN(15,file="output-A-matrix",status="unknown")
OPEN(17,file="output-fes",status="unknown")
OPEN(19,file="fes100",status="unknown")
!OPEN(11,file="free_energy-Cg.dat",status="unknown")

ALLOCATE(IPIV(M*N))
ALLOCATE(sum(M*N))
ALLOCATE(sum_g(100*100))
ALLOCATE(x(M*N,2),fes(M*N),fes_temp(M*N),Y(M*N))

!DO i=1,M
!  DO j=1,N
!i_temp=(i-1)*N+j
! READ(11,*) x(i_temp,1),x(i_temp,2),fes(i_temp)
!  END DO
!READ(11,*)
!END DO

Do i=1, M*N
 READ(11,*) x(i,1),x(i,2),fes(i)
  IF( x(i,1) .gt.  3.14d0)  x(i,1) = x(i,1) - 6.28d0
     IF( x(i,1) .lt. -3.14d0 ) x(i,1) = x(i,1) + 6.28d0
     IF( x(i,2) .gt.  3.14d0)  x(i,2) = x(i,2) - 6.28d0
     IF( x(i,2) .lt. -3.14d0 ) x(i,2) = x(i,2) + 6.28d0
END DO

ALLOCATE(A(M*N,M*N),A_temp(M*N,M*N))
fes_temp=fes

DO i=1,M*N
  Do j=1,M*N
        diff_1=x(i,1)-x(j,1)
          if (diff_1 .gt. 3.14d0 ) diff_1 =diff_1 - 6.28d0
          if (diff_1 .lt.-3.14d0 ) diff_1 =diff_1 + 6.28d0
        diff_2=x(i,2)-x(j,2)
          if (diff_2 .gt. 3.14d0 ) diff_2 =diff_2 - 6.28d0
          if (diff_2 .lt.-3.14d0 ) diff_2 =diff_2 + 6.28d0
       dum=(diff_1)**2/width_x + (diff_2)**2/width_y
    A(i,j)=dexp(-0.5d0*dum)
  END DO
END DO

A_temp=A

DO i=1,M*N
 WRITE(15,*),(A(i,j), J =1, M*N)
END DO

IPIV=0
CALL dgesv(M*N,1,A,M*N,IPIV,fes,M*N,INFO)
!CALL dposv('L',M*N,1,A,M*N,fes,M*N,INFO)
!CALL dgev('L',M*N,1,A,M*N,fes,M*N,INFO)

CALL dgemv('N',M*N,M*N,1.d0,A_temp,M*N,fes,1,0.d0,Y,1)
PRINT*,"INFO= ",INFO

!DO i=1,M*N
! PRINT*,(A(i,j), J =1, M*N)
!END DO


DO i=1,M*N
 WRITE(13,*), fes(i),fes_temp(i),Y(i),IPIV(i)
END DO

sum=0.d0

DO i=1,M*N
!sum=0.d0
  DO j=1,N*M
 sum(i)=sum(i)+fes(j)*A_temp(i,j)
  END DO
!WRITE(17,*) sum
END DO

DO i=1,M
  DO j=1,N
    WRITE(17,*),(i-1)*N+j, gridmin1+griddif1*dfloat(i-1),gridmin1+griddif1*dfloat(j-1),sum((i-1)*N+j)
  END DO
WRITE(17,*)
END DO

nbin=20
gridmin=-3.14d0
gridmax=3.14d0
griddif=(gridmax - gridmin)/dfloat(nbin-1)
 griddif=0.33d0!stop
sum_g=0.d0
!open(21,file="out.dat")
DO i=1,nbin
x_i=gridmin+griddif*dfloat(i-1) 
  DO k=1,nbin
    x_k=gridmin+griddif*dfloat(k-1) 
!write(21,*),x_i,x_k
       DO j=1,20*20

        diff_1=x(j,1)-x_i
          if (diff_1 .gt. 3.14d0 ) diff_1 =diff_1 - 6.28d0
          if (diff_1 .lt.-3.14d0 ) diff_1 =diff_1 + 6.28d0
        diff_2=x(j,2)-x_k
          if (diff_2 .gt. 3.14d0 ) diff_2 =diff_2 - 6.28d0
          if (diff_2 .lt.-3.14d0 ) diff_2 =diff_2 + 6.28d0
!print*,x(j,1),x(j,2)
          sum_g((i-1)*nbin+k)=sum_g((i-1)*nbin+k) + fes(j)*dexp(-0.5d0*((diff_1)**2/width_x + (diff_2)**2/width_y))
       END DO
   END DO
END DO

DO i=1,nbin
  DO j=1,nbin
    WRITE(19,*),(i-1)*nbin+j, gridmin+griddif*dfloat(i-1),gridmin+griddif*dfloat(j-1),sum_g((i-1)*nbin+j)
  END DO
WRITE(19,*)
END DO

DEALLOCATE(x,fes,A,IPIV)

CLOSE(11)
CLOSE(13)
CLOSE(15)
CLOSE(17)
END PROGRAM
