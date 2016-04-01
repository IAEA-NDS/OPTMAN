C     *******************************************************
C     START of privcom3  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK


C     FROM SUBROUTINE POTET        
      COMMON/POTEB/R,DE,VP,WP
      
!$OMP THREADPRIVATE(/POTEB/)
      

C     *******************************************************
C     END of pricom3
C     *******************************************************
