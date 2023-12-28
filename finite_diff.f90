PROGRAM check
IMPLICIT NONE

REAL*8,ALLOCATABLE:: s(:)!,z(:,:)
INTEGER::i,j,ir,nbin1,nbin2
!INTEGER,PARAMETER::dim=100
real*8:: dum

REAL*8::gridmin1, gridmax1, griddiff1
!REAL*8::gridmin2, gridmax2, griddiff2

!REAL*8, PARAMETER :: kb=3.16e-6  !in a.u K^-1 !1.9872041E-3
!REAL*8 :: T

OPEN(11,file="PROB_1D_wall_5.0_4.4")
OPEN(1,file="input-fd")

!READ(1,*) T
READ(1,*) gridmin1, gridmax1, griddiff1
!READ(1,*) gridmin2, gridmax2, griddiff2


nbin1 = NINT((gridmax1-gridmin1)/griddiff1)+1
!nbin2 = NINT((gridmax2-gridmin2)/griddiff2)+1

PRINT*,nbin1!,nbin2
!stop
ALLOCATE(s(nbin1))!,z(dim,dim))

DO i=1,nbin1
! do j=1,nbin2
  READ(11,*)dum,s(i)!,z(i,j)
! enddo
!READ(11,*)
 END DO


OPEN(3,file="mean_force-5.0_4.4.in")
DO i=2,nbin1-1
!DO j=2,nbin2-1
!write(3,*)gridmin1+(i-1)*griddiff1,gridmin2+(j-1)*griddiff2,-kb*T*dlog(max(1.d-32,(s(i+1,j)-s(i-1,j))/2.d0/griddiff1)), &
!                                                      -kb*T*dlog(max(1.d-32,(s(i,j+1)-s(i,j-1))/2.d0/griddiff2))
write(3,*)gridmin1+(i-1)*griddiff1,(s(i+1)-s(i-1))/2.d0/griddiff1
!END DO
!write(3,*)
END DO

END PROGRAM
