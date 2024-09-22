C     *******************************************************
      SUBROUTINE INMATLU
C     *******************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 

      complex*16,allocatable,dimension(:,:) :: A_MAT
      complex*16,allocatable,dimension(:,:) :: B_Z
      complex*16,allocatable,dimension(:,:) :: C_Z
      complex*16,allocatable,dimension(:,:) :: C_MM
      complex*16,allocatable,dimension(:)   :: WORK
      integer,allocatable,dimension(:)      :: IPIV
      integer i,j,info,error,k, KJ

      complex*16 Cinv,Cdir
      
      real*16 CRRR,CIII,DDR,DDI,NDDR,NDDI, QA_MATI, QA_MATR
      real*16 CRRR2,CIII2,DDR2,DDI2,NDDR2,NDDI2, QA_MATI2, QA_MATR2
    
      logical flag
      
      INCLUDE 'PRIVCOM.FOR'
      
      M=NCLL

      allocate(A_MAT(M,M),B_Z(M,M),C_Z(M,M),C_MM(M,M),
     * WORK(M),IPIV(M),stat=error)
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
      !print *, 'non MKL'
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
c
c     DDR=0.Q0
c     DDI=0.Q0
c     NDDR=0.Q0
c     NDDI=0.Q0
c     DDR2=0.Q0
c     DDI2=0.Q0
c     NDDR2=0.Q0
c     NDDI2=0.Q0
c     EPS=1.Q-1
c     DO j=1,NCLL
c         DO i=1,NCLL
c             CRRR=0.Q0
c             CIII=0.Q0
c             CRRR2=0.Q0
c             CIII2=0.Q0
c             DO k=1,NCLL
c                 KJ=(k-1)*NCLL+j
c                 IK=(i-1)*NCLL+k
c                 QA_MATR=DBLE(A_MAT(i,k))
c                 QA_MATI=DIMAG(A_MAT(i,k))
c                 QA_MATI2=DIMAG(A_MAT(k,j))
c                 QA_MATR2=DBLE(A_MAT(k,j))
c                 CRRR=CRRR+QA_MATR*ABR(KJ)
c    *                 -QA_MATI*ABI(KJ)
c                 CIII=CIII+QA_MATR*ABI(KJ)
c    *                 +QA_MATI*ABR(KJ)
c                 CRRR2=CRRR2+QA_MATR2*ABR(IK)
c    *                 -QA_MATI2*ABI(IK)
c                 CIII2=CIII2+QA_MATR2*ABI(IK)
c    *                 +QA_MATI2*ABR(IK)
c
c             END DO
c             IF(i.eq.j) then 
c                 DDR=DDR+ABS(CRRR-1.Q0)
c                 DDI=DDI+ABS(CIII)
c                 DDR2=DDR2+ABS(CRRR2-1Q0)
c                 DDI2=DDI2+ABS(CIII2)
c             else
c                 NDDR=NDDR+ABS(CRRR)
c                 NDDI=NDDI+ABS(CIII)
c                 NDDR2=NDDR2+ABS(CRRR2)
c                 NDDI2=NDDI2+ABS(CIII2)
c             end if
c             !if (ABS(CRRR)+ABS(CRRI).GT.EPS.and.i.ne.j) then
c             !    print *,'XA failed***********',i,j,ABS(CRRR),ABS(CIII)
c             !endif
c             !if (ABS(CRRR2)+ABS(CRRI2).GT.EPS.and.i.ne.j) then
c             !  print *,'AX failed***********',i,j,ABS(CRRR2),ABS(CIII2)
c             !endif
c         END DO
c     END DO
c     
c     flag = DDR.GT.EPS.OR.DDI.GT.EPS.OR.NDDR.GT.EPS.OR.
c    * NDDI.GT.EPS
c     
c     IF(flag) then
c         print *, '######## Direct Re and Im ##########'
c         print *, 'DDR=',DDR,' DDI=',DDI,
c    *     ' NDDR=',NDDR,    ' NDDI=',NDDI
c         print *, '####################################'
c     end if
c
c     IF(flag) then
c         print *, '####### Inverse Re and Im ##########'
c         print *, 'DDR=',DDR2,' DDI=',DDI2,
c    *     ' NDDR=',NDDR2,    ' NDDI=',NDDI2
c         print *, '####################################'
c     end if
c         
ccccccccccccccccccccccccccccccccccccccccccc      
c     
c     DO K=1,NCLL
c      K1=(K-1)*NCLL
c      DO L=1,NCLL
c       KL=K1+L
c       B_Z(K,L)=dcmplx(ABR(KL),ABI(KL))
c      END DO
c     END DO
c           
c     call ZGEMM('N', 'N', NCLL, NCLL, NCLL, 1.d0, A_MAT, NCLL, B_Z,
c    * NCLL, 0.d0, C_Z, NCLL)
c
c     DDR=0.Q0
c     DDI=0.Q0
c     NDDR=0.Q0
c     NDDI=0.Q0
c     DO iK=1,NCLL
c      DO j=1,NCLL
c             IF(i.eq.j) then 
c                 DDR=DDR+ABS(DBLE(C_Z(i,j))-1.Q0)
c                 DDI=DDI+ABS(DIMAG(C_Z(i,j)))
c             else
c                 NDDR=NDDR+ABS(DBLE(C_Z(i,j)))
c                 NDDI=NDDI+ABS(DIMAG(C_Z(i,j)))
c             end if
c      END DO
c     END DO
c     
c     IF(flag) then
c         print *, '########## Direct ZGEMM ############'
c         print *, 'DDR=',DDR,' DDI=',DDI,
c    *     ' NDDR=',NDDR,    ' NDDI=',NDDI
c         print *, '####################################'
c     end if
c     
c     C_Z=(0.d0,0.d0)
c     call ZGEMM('N', 'N', NCLL, NCLL, NCLL, 1.d0, B_Z, NCLL, A_MAT,
c    * NCLL, 0.d0, C_Z, NCLL)
c
c     DDR=0.Q0
c     DDI=0.Q0
c     NDDR=0.Q0
c     NDDI=0.Q0
c     DO iK=1,NCLL
c      DO j=1,NCLL
c             IF(i.eq.j) then 
c                 DDR=DDR+ABS(DBLE(C_Z(i,j))-1.Q0)
c                 DDI=DDI+ABS(DIMAG(C_Z(i,j)))
c             else
c                 NDDR=NDDR+ABS(DBLE(C_Z(i,j)))
c                 NDDI=NDDI+ABS(DIMAG(C_Z(i,j)))
c             end if
c      END DO
c     END DO
c     
c     IF(flag) then
c         print *, '########## Inverse ZGEMM ###########'
c         print *, 'DDR=',DDR,' DDI=',DDI,
c    *     ' NDDR=',NDDR,    ' NDDI=',NDDI
c         print *, '####################################'
c     end if
c     
cccccccccccccccccccccccccccccccccccccccccccc      
c     
c     C_MM = MATMUL(A_MAT, B_Z)
c     
c     DDR=0.Q0
c     DDI=0.Q0
c     NDDR=0.Q0
c     NDDI=0.Q0
c     DO iK=1,NCLL
c      DO j=1,NCLL
c             IF(i.eq.j) then 
c                 DDR=DDR+ABS(DBLE(C_MM(i,j))-1.Q0)
c                 DDI=DDI+ABS(DIMAG(C_MM(i,j)))
c             else
c                 NDDR=NDDR+ABS(DBLE(C_MM(i,j)))
c                 NDDI=NDDI+ABS(DIMAG(C_MM(i,j)))
c             end if
c      END DO
c     END DO      
c     
c     IF(flag) then
c         print *, '########## Direct MATMUL ###########'
c         print *, 'DDR=',DDR,' DDI=',DDI,
c    *     ' NDDR=',NDDR,    ' NDDI=',NDDI
c         print *, '####################################'
c     end if
c
c     C_MM = MATMUL(B_Z, A_MAT)
c     
c     DDR=0.Q0
c     DDI=0.Q0
c     NDDR=0.Q0
c     NDDI=0.Q0
c     DO iK=1,NCLL
c      DO j=1,NCLL
c             IF(i.eq.j) then 
c                 DDR=DDR+ABS(DBLE(C_MM(i,j))-1.Q0)
c                 DDI=DDI+ABS(DIMAG(C_MM(i,j)))
c             else
c                 NDDR=NDDR+ABS(DBLE(C_MM(i,j)))
c                 NDDI=NDDI+ABS(DIMAG(C_MM(i,j)))
c             end if
c      END DO
c     END DO      
c     
c     IF(flag) then
c         print *, '########## Inverse MATMUL ##########'
c         print *, 'DDR=',DDR,' DDI=',DDI,
c    *     ' NDDR=',NDDR,    ' NDDI=',NDDI
c         print *, '####################################'
c     end if
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     
c     DDR=0.Q0
c     DDI=0.Q0
c     NDDR=0.Q0
c     NDDI=0.Q0
c     DDR2=0.Q0
c     DDI2=0.Q0
c     NDDR2=0.Q0
c     NDDI2=0.Q0
c     DO j=1,NCLL
c         DO i=1,NCLL
c             Cdir =(0.d0,0.d0)
c             Cinv=(0.d0,0.d0)
c             DO k=1,NCLL
c                 Cdir=Cdir+A_MAT(i,k)*B_Z(k,j)
c                 Cinv=Cinv+B_Z(i,k)*A_MAT(k,j)
c             END DO
c             IF(i.eq.j) then 
c                 DDR=DDR+ABS(DBLE(Cdir-(1.d0,0.d0)))
c                 DDI=DDI+ABS(DIMAG(Cdir))
c                 DDR2=DDR2+ABS(DBLE(Cinv-(1.d0,0.d0)))
c                 DDI2=DDI2+ABS(DIMAG(Cinv))
c             else
c                 NDDR=NDDR+ABS(DBLE(Cdir))
c                 NDDI=NDDI+ABS(DIMAG(Cdir))
c                 NDDR2=NDDR2+ABS(DBLE(Cinv))
c                 NDDI2=NDDI2+ABS(DIMAG(Cinv))                  
c             end if
c         END DO
c     END DO
c     
c     IF(flag) then
c         print *, '########### Direct Cplx ############'
c         print *, 'DDR=',DDR,' DDI=',DDI,
c    *     ' NDDR=',NDDR,    ' NDDI=',NDDI
c         print *, '####################################'
c     end if
c
c     IF(flag) then
c         print *, '########### Inverse Cplx ###########'
c         print *, 'DDR=',DDR2,' DDI=',DDI2,
c    *     ' NDDR=',NDDR2,    ' NDDI=',NDDI2
c         print *, '####################################'
c     end if      
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      !PRINT *, '!!!INVERSION ACCURACY = ', DDD

      
      
      DO 2 K=1,NCLL
       K1=(K-1)*NCLL
       DO 2 L=1,NCLL
        KL=K1+L
        ABR(KL)=REAL(A_MAT(K,L))
        ABI(KL)=AIMAG(A_MAT(K,L))
    2 CONTINUE       

      deallocate(A_MAT,B_Z,C_Z,C_MM,IPIV,WORK,stat=error)
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