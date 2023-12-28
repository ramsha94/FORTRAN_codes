PROGRAM fitting_mean_force

  USE mympi
  IMPLICIT NONE
  REAL*8 :: width,cutoff
  INTEGER :: np, ncv, ncpu, icpu, i
  REAL*8, POINTER :: a(:,:), b(:), grid(:,:),mf(:,:)
  LOGICAL, POINTER :: nn_list(:,:)
  LOGICAL :: ionode
  CHARACTER (LEN=15) :: myfmt

  INTEGER,POINTER:: IPIV(:)
  ionode=.FALSE.
  CALL MPI_Start
  CALL MPI_get_ncpu(ncpu)
  CALL MPI_get_cpuid(icpu)
  IF(icpu==0)ionode=.TRUE.
  IF(ionode)THEN
    PRINT *, "Number of CPUs      =", ncpu

    PRINT *, 'ENTER FUNCTION WIDTH'
    READ(*,*) width

    PRINT *, 'ENTER NUMBER OF CVs'
    READ(*,*) ncv

    PRINT *, 'ENTER CUTOFF FACTOR (cutoff*width will be the actual cutoff)'
    READ(*,*) cutoff
  END IF
  CALL MPI_IBcast1(ncv) ; CALL MPI_RBcast1(width) ; CALL MPI_RBcast1(cutoff)
 
  cutoff=cutoff*width
  
  IF(ncv.eq.0)STOP 'ncv cannot be zero'

  !print *, 'reading ', icpu
  CALL read_mean_force
  print *, 'read...', SIZE(grid,1), SIZE(grid,2), SIZE(mf,1), SIZE(mf,2)
  !print *, 'done reading ', icpu

  IF(ionode)THEN
    PRINT *, 'number of given points =', np
    PRINT *, 'number of CVs          =', ncv
    PRINT *, 'width                  =', width
    PRINT *, 'cutoff                 =', cutoff
  END IF
  
 
   print *, 'call prepare-nn-list', icpu
  CALL prepare_nn_list
   print *, 'done prepare-nn-list', SIZE(nn_list,1), SIZE(nn_list,2)
  CALL prepare_lin_eq(width) ! compute a, and b matrices
  print *, 'going to solve'
  PRINT *, 'a size (outside): ', SIZE(a,1), size(a,2)
!  CALL solve_linear_equation ! solve ax=b 
!  print *, 'done  solve'

  CALL solve_linear_equation_dgel ! solve ax=b 
  print *, 'goting to print '
  IF(ionode)THEN
    OPEN(100,file='weights.log')
  print *, 'goting to get_fmt '
    CALL get_fmt(myfmt)
  print *, 'starting to write ', myfmt
    DO i=1,np
      WRITE(100,trim(myfmt))grid(1:ncv,i),width,cutoff,-b(i)
    END DO
  print *, 'done  write '
    CLOSE(100)
    WRITE(*,*)'Finished writing weights!'
  END IF

  CALL MPI_Stop

CONTAINS

  SUBROUTINE get_fmt(myfmt)
     IMPLICIT NONE
     CHARACTER(LEN=*) :: myfmt
     SELECT CASE(ncv)
     CASE(1)
       myfmt='(4E16.6)'
     CASE(2)
       myfmt='(5E16.6)'
     CASE(3)
       myfmt='(6E16.6)'
     CASE(4)
       myfmt='(7E16.6)'
     CASE DEFAULT
       myfmt='*'
     END SELECT
  END SUBROUTINE get_fmt

  SUBROUTINE prepare_nn_list
     IMPLICIT NONE
     REAL*8 :: dx
     REAL*8, POINTER :: xi(:), xj(:), dxij(:)
     INTEGER :: i, j
     
     ALLOCATE(nn_list(np,np))
     ALLOCATE(xi(ncv),xj(ncv),dxij(ncv))
     nn_list(:,:)=.FALSE.
     DO i=1,np
       xi(:)=grid(:,i)
       DO j=i,np
         xj(:)=grid(:,j)
         dxij(:)=difference(xi,xj)
         dx=DOT_PRODUCT(dxij(:),dxij(:))
         IF(DSQRT(dx)<=cutoff)THEN
           nn_list(i,j)=.TRUE.
           nn_list(j,i)=.TRUE.
         END IF
       END DO
     END DO
     DEALLOCATE(xi,xj,dxij)
     print *, 'done prepare_nn_list', icpu
  END SUBROUTINE prepare_nn_list

  SUBROUTINE solve_linear_equation
     IMPLICIT NONE
!     INTEGER,POINTER:: IPIV(:)
     INTEGER :: info
     REAL*8 :: time1

     time1=start_clock1()

     ALLOCATE(IPIV(np))
     IPIV=0
!TODO USE Scalapack 
     print *, 'calling dgesv'
     print *, 'size a', SIZE(a,1), SIZE(a,2)
     print *, 'size ipv', SIZE(ipiv,1)
     print *, 'size b', SIZE(b,1)
     CALL dgesv(np,1,a,np,IPIV,b,np,INFO)  !to get solutions of x matrix in fes
     print *, 'done dgesv'
     DEALLOCATE(ipiv)
     time1=stop_clock(time1)
     CALL print_cputime(time1,'solve_linear_equation')
  END SUBROUTINE solve_linear_equation

  SUBROUTINE solve_linear_equation_dgel
     IMPLICIT NONE     
!     INTEGER ::          M, N, NRHS
     INTEGER,PARAMETER :: M = 800, N = 400, NRHS = 1 
!     INTEGER   ::        LDA, LDB
     INTEGER,PARAMETER :: LDA = M, LDB = M 
!     INTEGER  ::        LWMAX
     INTEGER,PARAMETER  ::      LWMAX = 100 
     REAL*8 :: time1
     INTEGER :: INFO, LWORK
     REAL*8,ALLOCATABLE :: WORK(:) 
     time1=start_clock1()

     ALLOCATE(IPIV(np))
     IPIV=0
!TODO USE Scalapack 
     print *, 'calling dgelsv'
     print *, 'size a', SIZE(a,1), SIZE(a,2)
     print *, 'size ipv', SIZE(ipiv,1)
     print *, 'size b', SIZE(b,1)
     LWORK = 2*SIZE(a,2)
     ALLOCATE(WORK(LWORK))
!     CALL dgesv(np,1,a,np,IPIV,b,np,INFO)  !to get solutions of x matrix in fes
!     CALL DGELS( 'No transpose', M, N, NRHS, a, LDA, b, LDB, WORK,LWORK, INFO )
      CALL DGELS('N', SIZE(a,1),SIZE(a,2),1,a,SIZE(a,1),b,SIZE(b,1),WORK,LWORK,INFO )
     print *, 'done dgelsv'
     DEALLOCATE(ipiv)
     time1=stop_clock(time1)
     CALL print_cputime(time1,'solve_linear_equation')
  END SUBROUTINE solve_linear_equation_dgel
!
  SUBROUTINE prepare_lin_eq(ds)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: ds
  INTEGER :: i, j, k, ij, ii, kk
  REAL*8, POINTER :: xi(:), xj(:), xk(:),df(:), deriv(:,:),d(:)
  REAL*8 :: time1,time2,time3

print*,"check"
  ALLOCATE(xi(ncv),xj(ncv),xk(ncv),df(ncv))
  ALLOCATE(d(ncv))
!   Compute A and b for Ax=b 

  time1=start_clock1()

  ALLOCATE(a(ncv*np,np))
  a(:,:)=0.d0
  DO i=1,np
   xi(:)=grid(:,i)
   ii=0
   DO k=1,np
     IF(nn_list(i,k))THEN
       xk(:)=grid(:,k)       
       !d(:)=d_vonmises(ds,xk,xi)
       d(:)=d_gaussian(ds,xk,xi)   !For using Gaussians to get weights and fit F(s) !RJ
!       print*,"Working Check"
       DO kk=1,ncv
         a(ii+kk,i)=d(kk)
       END DO
       ii=ii+ncv
     ELSE
       ii=ii+ncv
     END IF
   END DO
  END DO 
print*,"A generated"
  print *, 'prepare lin eq', icpu

  ALLOCATE(b(ncv*np))

  time3=start_clock1()

  ii=0
  DO i=1,np
    d(:)=-mf(:,i)
!    ii=ii+ncv
    do kk=1,ncv
      b(ii+kk)=d(kk)
!      print*,"ii+kk",ii+kk
    END DO
    ii=ii+ncv
!    print*,"ii",ii
  END DO
  time3=stop_clock(time3)
  time1=stop_clock(time1)
print*,"B generated"
  DEALLOCATE(xi,xj,xk,df)
  CALL print_cputime(time1,'prepare_lin_eq')
  CALL print_cputime(time2,'setup A matrix')
  CALL print_cputime(time3,'setup b matrix')
  print *, 'size a (end)=', size(a,1),size(a,2)
  print *, 'done prepare lin eq', icpu

  END SUBROUTINE prepare_lin_eq
!
  SUBROUTINE parallel_grid(k,kmin,kmax)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: k
    INTEGER, INTENT(OUT) ::  kmin, kmax
    INTEGER :: k_r, k_p
    k_r=MOD(k,ncpu)  !reminder of grid after equal partition
    k_p=k/ncpu       !chunk of the grid per processor
    kmin=icpu*k_p+1  !minimum value of the grid
    kmax=kmin+k_p-1  !maximum value of the grid
    IF(icpu==ncpu-1)kmax=kmax+k_r ! Any remaining grids added to the last processor
  END SUBROUTINE parallel_grid

!
  FUNCTION d_gaussian(ds,xi,xj) RESULT(dg)
    IMPLICIT NONE
    REAL*8, INTENT(IN), POINTER :: xi(:), xj(:)
    REAL*8, INTENT(IN) :: ds
    REAL*8 :: dg(SIZE(xi))
    REAL*8 :: e
!  
    dg(:)=-difference(xi(:),xj(:))/ds/ds*gaussian(ds,xi,xj)
  

  END FUNCTION d_gaussian
!
  FUNCTION gaussian(ds,xi,xj) RESULT(f)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: ds
    REAL*8, INTENT(IN),POINTER :: xi(:), xj(:)
    REAL*8 :: f
    REAL*8 :: e
    REAL*8 :: dx(ncv)

    dx(:)=difference(xi(:),xj(:))
    e=DOT_PRODUCT(dx,dx)
    f=DEXP(-0.5*e/ds/ds)
  END FUNCTION gaussian

!
  FUNCTION vonmises(ds,xi,xj) RESULT(f)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: ds
    REAL*8, INTENT(IN),POINTER :: xi(:), xj(:)
    REAL*8 :: f

    REAL*8 :: variance
    INTEGER :: i 
    REAL*8 :: e
    REAL*8 :: dx(ncv),expn(ncv)
    REAL*8,PARAMETER :: factor=1.0d0
  
    variance = ds**2    
    dx(:)=difference(xi(:),xj(:))
     
     e=0.d0
     DO i=1,size(dx)
      e=e+(cos(factor*dx(i))-1.0d0)/(factor*factor*variance)
     END DO
     f =exp(e)
  END FUNCTION vonmises

!
  FUNCTION d_vonmises(ds,xi,xj) RESULT(dg)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: ds
    REAL*8, INTENT(IN),POINTER :: xi(:), xj(:)
    REAL*8 :: variance
    REAL*8 :: dg(SIZE(xi))

    REAL*8 :: deriv(SIZE(xi))
    REAL*8 :: dx(SIZE(xi)),expn(SIZE(xi))
    REAL*8,PARAMETER :: factor=1.0d0

    variance = ds**2 

    dx(:)=difference(xi(:),xj(:))
    DO i=1,size(dx)
      deriv(i)= -sin(factor*dx(i))/(factor*variance)
    END DO
    dg(:) = vonmises(ds,xi,xj)*deriv(:)
!    print*,"xi,xj,difference",xi(:),xj(:),dx(:)
!    print*,"dg",dg(:) 
  END FUNCTION d_vonmises

!
  FUNCTION difference(xi,xj) RESULT(diff)
    IMPLICIT NONE
    REAL*8, INTENT(IN), POINTER :: xi(:), xj(:)
    REAL*8 :: diff(SIZE(xi))

    REAL*8, PARAMETER :: twopi=8.d0*ATAN(1.d0)
    INTEGER :: i

    diff(:)=xj(:)-xi(:)
!    DO i=1,SIZE(diff)
!      IF(diff(i)>twopi)diff(i)=diff(i)-twopi
!      IF(diff(i)>twopi)diff(i)=diff(i)-twopi
!    END DO 
  END FUNCTION difference

!
  SUBROUTINE read_mean_force
    IMPLICIT NONE
    INTEGER :: i

    IF(ionode)THEN
      OPEN(100,FILE='mean_force.in')
      CALL get_lines(100,np)    ; REWIND(100)
      if(np.eq.0)STOP 'np cannot be zero'
    END IF
    CALL MPI_IBcast1(np)
!    CALL get_columns(100,ncv) ; REWIND(100) ; ncv=ncv/2 ! grid_x, grid_y, mf_x,
!    mf_y
    ALLOCATE(grid(ncv,np)) ; ALLOCATE(mf(ncv,np))
    IF(ionode)THEN
      DO i=1,np
       READ(100,*)grid(:,i), mf(:,i)
      END DO
   !   READ(*,*)
      CLOSE(100)
    END IF
    CALL MPI_RBcast(grid,ncv*np)
    CALL MPI_RBcast(mf,ncv*np)
    print *, 'done read -', icpu

  END SUBROUTINE read_mean_force

!
  SUBROUTINE get_lines(ir,nl)
    IMPLICIT NONE

    INTEGER :: ir, nl
    INTEGER :: ios

    nl=0
    DO
      READ(ir,*,IOSTAT=ios)
      IF(ios.ne.0)EXIT
        nl=nl+1
    END DO
  END SUBROUTINE get_lines

!
  SUBROUTINE get_columns(ir,nc)
    IMPLICIT NONE
    INTEGER :: ir, nc
    INTEGER :: ios,i
    INTEGER, PARAMETER :: max_column_width=120
    CHARACTER (LEN=max_column_width) :: line
    REAL*8 :: dummy(max_column_width)
    READ(ir,'(A120)')line
    print *, 'line =',line
    do i=1,max_column_width
      READ(line, *, iostat=ios)dummy(1:i)
      if(ios==0) exit
      print *, i, dummy(1:i)
    enddo
    nc=i-1
  END SUBROUTINE get_columns


  FUNCTION start_clock1() RESULT(time)
    IMPLICIT NONE
    REAL*8 :: time
    CALL cpu_time(time)
  END FUNCTION start_clock1

  FUNCTION stop_clock(time0) RESULT(time)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: time0
    REAL*8 :: time
    CALL cpu_time(time)
    time=time-time0
  END FUNCTION stop_clock

  SUBROUTINE print_cputime(time,routine)
   IMPLICIT NONE
   REAL*8 :: time
   CHARACTER(LEN=*) :: routine
   PRINT '(3A,F6.3)', 'CPU Time [',TRIM(routine),'] (s)=',time
  END SUBROUTINE print_cputime 
END PROGRAM fitting_mean_force

 
