C     *******************************************************
      SUBROUTINE INMATLU
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 

      complex*16,allocatable,dimension(:,:) :: A_MAT
      complex*16,allocatable,dimension(:)   :: WORK
      integer,allocatable,dimension(:)      :: IPIV
      integer i,j,info,error,k

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
  
      call ZGETRF(M,M,A_MAT,M,IPIV,info)
      if(info.ne.0) then
       write(*,*)"!!!!!!!!!!!!!!!ZGETRF failed"
      end if
  
      call ZGETRI(M,A_MAT,M,IPIV,WORK,M,info)
      if(info .ne. 0) then
       write(*,*)"!!!!!!!!!!!!!!!!!ZGETRI failed"
      end if
  
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