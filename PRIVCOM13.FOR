C     *******************************************************
C     START of privcom13  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK


C     FROM SUBROUTINE OVLAG
C     COMMON/ENA/EN,EL(40),BET(10),NUR,NMAX,NPD,LAS 
      COMMON/ENAa/NPD,LAS 
      COMMON/ENA/EN,EL(40),NMAX 
      COMMON/ENAi/BET(10),NUR
 
!$OMP THREADPRIVATE(/ENA/,/ENAi/,/ENAa/)
 
C     *******************************************************
C     END of privcom13
C     *******************************************************
