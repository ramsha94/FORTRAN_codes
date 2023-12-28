MODULE mympi
  CONTAINS

  SUBROUTINE MPI_Start()
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: i_err
!#if defined (_MPI)
  CALL mpi_init(i_err)
!#endif
  END SUBROUTINE MPI_Start

  SUBROUTINE MPI_Stop()
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: i_err
!#if defined (_MPI)
  call MPI_FINALIZE(i_err)
!#endif
  END SUBROUTINE MPI_Stop

  SUBROUTINE MPI_get_ncpu(ncpu)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: ncpu, i_err
  ncpu=1
!#if defined (_MPI)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,i_err)
!#endif
  END SUBROUTINE MPI_get_ncpu

  SUBROUTINE MPI_get_cpuid(icpu)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: icpu, i_err
  icpu=0
!#if defined (_MPI)
  CALL  MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
  END SUBROUTINE MPI_get_cpuid

  SUBROUTINE MPI_IBcast(myint,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: leng, myint(*), i_err
!#if defined (_MPI)
  CALL MPI_BCAST(myint,leng,MPI_INTEGER,0,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_IBcast

  SUBROUTINE MPI_IBcast1(myint)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: leng, myint, i_err
!#if defined (_MPI)
  leng=1
  CALL MPI_BCAST(myint,leng,MPI_INTEGER,0,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_IBcast1


  SUBROUTINE MPI_RBcast(myreal,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  REAL*8 :: myreal(*)
  INTEGER :: leng, i_err
!#if defined (_MPI)
  CALL MPI_BCAST(myreal,leng,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_RBcast

  SUBROUTINE MPI_RBcast1(myreal)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  REAL*8 :: myreal
  INTEGER :: leng, i_err
!#if defined (_MPI)
  leng=1
  CALL MPI_BCAST(myreal,leng,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_RBcast1


  SUBROUTINE MPI_Sync_Procs
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER i_err
!#if defined (_MPI)
  CALL MPI_Barrier(MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_Sync_Procs



  SUBROUTINE Set_IOnode(ionode)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  LOGICAL :: ionode
  INTEGER :: icpu, i_err
  ionode=.false.
  icpu=0
!#if defined (_MPI)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
  IF(icpu.EQ.0)ionode=.TRUE.
  END SUBROUTINE Set_IOnode

  SUBROUTINE MPI_GlobSumR1(myreal_in,myreal_out)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  REAL*8 :: myreal_in, myreal_out
  INTEGER :: leng,i_err
  leng=1
!#if defined (_MPI)
  CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_GlobSumR1


  SUBROUTINE MPI_GlobSumR(myreal_in,myreal_out,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  REAL*8 :: myreal_in(*), myreal_out(*)
  INTEGER :: leng,i_err
!#if defined (_MPI)
  CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_GlobSumR

  SUBROUTINE MPI_GlobSumI1(myint_in,myint_out)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: myint_in, myint_out
  INTEGER :: leng,i_err
!#if defined (_MPI)
  leng=1
  CALL MPI_Allreduce(myint_in,myint_out,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_GlobSumI1



  SUBROUTINE MPI_GlobSumI(myint_in,myint_out,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: myint_in(*), myint_out(*)
  INTEGER :: leng,i_err
!#if defined (_MPI)
  CALL MPI_Allreduce(myint_in,myint_out,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_GlobSumI

  SUBROUTINE MPI_GlobSumC(mycmplx_in,mycmplx_out,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  COMPLEX*16 :: mycmplx_in(*), mycmplx_out(*)
  INTEGER :: leng,i_err
!#if defined (_MPI)
  CALL MPI_Allreduce(mycmplx_in,mycmplx_out,leng,MPI_DOUBLE_COMPLEX, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_GlobSumC

END MODULE mympi
