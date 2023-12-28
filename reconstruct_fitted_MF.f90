PROGRAM fitting_mean_force

  USE mympi
  IMPLICIT NONE

  REAL*8 :: width,cutoff
  INTEGER :: np, ncv, ncpu, icpu, icv, ifunc
  REAL*8, POINTER ::  b(:), grid(:,:), gridmin(:), gridmax(:), gridw(:)
  INTEGER, POINTER :: ngrid(:)
  LOGICAL :: ionode
  CHARACTER (LEN=15) :: myfmt

  ionode=.FALSE.
  CALL MPI_Start
  CALL MPI_get_ncpu(ncpu)
  CALL MPI_get_cpuid(icpu)
  IF(icpu==0)ionode=.TRUE.
  IF(ionode)THEN
    PRINT *, "Number of CPUs      =", ncpu

    PRINT *, 'ENTER THE FUNCTION CODE: 1 - Von Mises ; 2 - Gaussian'
    READ(*,*) ifunc

    PRINT *, 'ENTER NUMBER OF CVs'
    READ(*,*) ncv
  END IF
  CALL MPI_IBcast1(ncv) ; CALL MPI_IBcast1(ifunc)
 
  IF(ncv.eq.0)STOP 'ncv cannot be zero'

  !print *, 'reading ', icpu
  CALL read_weights
  !print *, 'read...', SIZE(grid,1), SIZE(grid,2), SIZE(b,1)
  !print *, 'done reading ', icpu

  ALLOCATE(ngrid(ncv))  ;  ALLOCATE(gridmin(ncv)) 
  ALLOCATE(gridmax(ncv));  ALLOCATE(gridw(ncv)) 

  DO icv=1,ncv
    PRINT *, 'ENTER grid min, grid max,and grid width for cv=',icv
    READ(*,*) gridmin(icv),gridmax(icv),gridw(icv)
    ngrid(icv)=INT((gridmax(icv)-gridmin(icv))/gridw(icv))+1
!    PRINT*,ngrid
  END DO
!ngrid(1)=60
!ngrid(2)=19
  print *, 'ngrids =', ngrid(:)
  DO icv=1,ncv
    IF(ngrid(icv)<=0)STOP 'error..ngrid<=0'
  END DO

  IF(ionode)THEN
    PRINT *, 'number of given points =', np
    PRINT *, 'number of CVs          =', ncv
    PRINT *, 'width                  =', width
    PRINT *, 'cutoff                 =', cutoff
  END IF
  
  CALL print_potential
  CALL MPI_Stop

CONTAINS
  SUBROUTINE print_potential
    IMPLICIT NONE
    INTEGER :: i, j, k, kmin, kmax
    REAL*8  :: g,g1
    REAL*8, POINTER :: xi(:), xk(:)
    
    CALL parallel_grid(np,kmin,kmax)
    print *, 'kmin =', kmin, ' kmax=',kmax
!    print *, 'size b =', size(b)
    ALLOCATE(xi(ncv),xk(ncv))
    IF(ionode)OPEN(100,file='free_energy.dat',status="replace")

    IF(ncv==2)THEN
      DO i=1,ngrid(1)
        xi(1)=gridmin(1)+DFLOAT(i-1)*gridw(1)
        !print *, 'xi created i=', i
        DO j=1,ngrid(2)
          xi(2)=gridmin(2)+DFLOAT(j-1)*gridw(2)
        !print *, 'xi(2) created j=', j
          g=0.d0
          DO k=kmin,kmax
            !print *, 'i =',i, 'j=',j,'k=',k
            xk(:)=grid(:,k)
            !print *, 'xk copied'
            g=g+b(k)*gaussian(width,xi,xk)   !For fitting with gaussian RJ
        !    g=g+b(k)*vonmises(width,xk,xi)
            !print *, 'g-calculated'
          END DO
          g1=0.d0
!          CALL MPI_GlobSumR1(g,g1)
            !print *, 'glosum calculated'
         ! IF(ionode)
 !        print*,"grid-FES-V",xi(1),xi(2),g
        ! WRITE(100,'(3F16.6)')xi(1),xi(2),g
         WRITE(100,*)xi(1),xi(2),g
            !print *, 'printed '
        END DO
        !IF(ionode)
        WRITE(100,*)
      END DO
    ELSE
       STOP 'ncv .ne.2,not implemented'
    END IF

    IF(ionode)CLOSE(100)

  END SUBROUTINE print_potential


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

!  time2=start_clock1()
!  time2=stop_clock(time2)
!  CALL print_cputime(time2,'setup A matrix')


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
    REAL*8 :: e,dd
    REAL*8, POINTER :: dx(:)
    REAL*8,PARAMETER :: factor=1.0d0

    ALLOCATE(dx(ncv))
  
    variance = ds*ds    
    dx(:)=difference(xi(:),xj(:))

    dd=DSQRT(DOT_PRODUCT(dx(:),dx(:)))
    
    f=0.d0
    IF(dd<=cutoff)THEN 
      e=0.d0
      DO i=1,size(dx)
        e=e+(cos(factor*dx(i))-1.0d0)/(factor*factor*variance)
      END DO
      f =exp(e)
!      print*,"x(i),x(j),VM_val",xi(:),xj(:),f
    END IF
    DEALLOCATE(dx)
  END FUNCTION vonmises

!
!
  FUNCTION difference(xi,xj) RESULT(diff)
    IMPLICIT NONE
    REAL*8, INTENT(IN), POINTER :: xi(:), xj(:)
    REAL*8 :: diff(SIZE(xi))

    REAL*8, PARAMETER :: twopi=8.d0*ATAN(1.d0)
    INTEGER :: i

    diff(:)=xj(:)-xi(:)
  !  DO i=1,SIZE(diff)
  !    IF(diff(i)>twopi)diff(i)=diff(i)-twopi
  !    IF(diff(i)>twopi)diff(i)=diff(i)-twopi
   ! END DO 
  END FUNCTION difference

!
  SUBROUTINE read_weights
    IMPLICIT NONE
    INTEGER :: i

    IF(ionode)THEN
      OPEN(100,FILE='weights.log')
      CALL get_lines(100,np)    ; REWIND(100)
      if(np.le.0)STOP 'np cannot be zero'
    END IF

    CALL MPI_IBcast1(np)

    ALLOCATE(grid(ncv,np)) ; ALLOCATE(b(np))

    IF(ionode)THEN
      DO i=1,np
        READ(100,*)grid(1:ncv,i),width,cutoff,b(i)
      END DO
      CLOSE(100)
    END IF
    CALL MPI_RBcast(grid,ncv*np)
    CALL MPI_RBcast(b,np)

    !print *, 'done read -', icpu

  END SUBROUTINE read_weights

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
    !print *, 'line =',line
    do i=1,max_column_width
      READ(line, *, iostat=ios)dummy(1:i)
      if(ios==0) exit
      !print *, i, dummy(1:i)
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

 
