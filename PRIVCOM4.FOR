C     *******************************************************
C     START of privcom4  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK

C     FROM SUBROUTINE RACAH  
      COMMON/LOFAC/A(800)
       
!$OMP THREADPRIVATE(/LOFAC/)


C     *******************************************************
C     END of pricom4
C     *******************************************************
