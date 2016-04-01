C     *******************************************************
C     START of privcom20  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK



C     C     FROM SUBROUTINE KNDITD        
      COMMON/WEIGHT/NKIWE(40),KWE(16,40),WEIGH(16,40)
      COMMON/EFFD/EFFDEF(40,40,6),DEFNUL
      
!$OMP THREADPRIVATE(/WEIGHT/,/EFFD/) 
      

C     *******************************************************
C     END of pricom20
C     *******************************************************
