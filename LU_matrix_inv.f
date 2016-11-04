C     *******************************************************
      SUBROUTINE INMATLU
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 

      complex*16,allocatable,dimension(:,:) :: A_MAT
      complex*16,allocatable,dimension(:)   :: WORK
      integer,allocatable,dimension(:)      :: IPIV
      integer i,j,info,error,k, KJ
      
      real*16 CRRR,CIII,DDD, QA_MATI, QA_MATR
      complex*16 MYY
      
      
      INCLUDE 'PRIVCOM.FOR'
      
      M=NCLL

      allocate(A_MAT(M,M),WORK(M),IPIV(M),stat=error)
      if (error.ne.0)then
        print *,"!!!!!!!!!!!!!error:not enough memory"
        stop
      end if
  
      DO 1 K=1,NCLL
       K1=(K-1)*NCLL
       DO 1 L=1,NCLL
        KL=K1+L
        A_MAT(K,L)=dcmplx(ABR(KL),ABI(KL))
    1 CONTINUE  
#ifdef LAPACK
      call ZGETRF(M,M,A_MAT,M,IPIV,info)
#else
      print *, 'non MKL'
#endif
      
      if(info.ne.0) then
       write(*,*)"!!!!!!!!!!!!!!!ZGETRF failed"
      end if
  
#ifdef LAPACK
      call ZGETRI(M,A_MAT,M,IPIV,WORK,M,info)
#else
      print *, 'non MKL'
#endif
      if(info .ne. 0) then
       write(*,*)"!!!!!!!!!!!!!!!!!ZGETRI failed"
      end if
  
c     CRRR=0.Q0
c     CIII=0.Q0
c     DDD=0.Q0
c     MYY=dcmplx(0,0)
c     DO j=1,NCLL
c         DO i=1,NCLL
c             DO k=1,NCLL
c                 KJ=(k-1)*NCLL+j
c                   C=C+A_MAT(i,k)*A(k,j)
c                 QA_MATR=DBLE(A_MAT(i,k))
c                 QA_MATI=DIMAG(A_MAT(i,k))
c                 CRRR=CRRR+QA_MATR*ABR(KJ)
c    *                 -QA_MATI*ABI(KJ)
c                 CIII=CIII+QA_MATR*ABI(KJ)
c    *                 +QA_MATI*ABR(KJ)
c                 MYY = MYY + A_MAT(i,k)*dcmplx(ABR(KJ),ABI(KJ))
c             END DO
c             IF(i.eq.j) CRRR=CRRR-1.Q0
c             IF(i.eq.j) MYY=MYY-dcmplx(1,0)
c             DDD = DDD + CDABS(MYY)!QABS(CRRR)+QABS(CIII)!CRRR*CRRR + CIII*CIII
c             CRRR=0.Q0
c             CIII=0.Q0
c         END DO
c     END DO
c     
c     PRINT *, '!!!INVERSION ACCURACY = ', DDD

      
      
      DO 2 K=1,NCLL
       K1=(K-1)*NCLL
       DO 2 L=1,NCLL
        KL=K1+L
        ABR(KL)=REAL(A_MAT(K,L))
        ABI(KL)=AIMAG(A_MAT(K,L))
    2 CONTINUE       

      deallocate(A_MAT,IPIV,WORK,stat=error)
      if (error.ne.0)then
        print *,"!!!!!!!!!!!!!!!!!!!error:fail to release"
        stop
      end if   
      
      
!      IF(MEPRI.LT.98) PRINT 24,SUMN,SUMD,INFOR,ITER
!      WRITE(21,24)SUMN,SUMD,INFOR,ITER
!   24 FORMAT(10X,'WARNING! MATRIX IS POORLY INVERTED'/
!     *5X,'SUM OF NON-DIAG. ELEM-S=',D11.5,
!     *', DIAGONAL=',D11.5,',INFOR=',I2,',ITER=',I2)
     
      RETURN
      END      