C     *******************************************************
C     START of privcom14  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK


C     FROM SUBROUTINE COPHA
      COMMON/COUL/CONZ,ETA,COPH(40,90)
      
       
!$OMP THREADPRIVATE(/COUL/)
 
C     *******************************************************
C     END of privcom14
C     *******************************************************
