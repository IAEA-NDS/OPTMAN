C     *******************************************************
C     START of privcom6  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK


C     FROM SUBROURINE RIPAT
      COMMON/POT1/VR3,VR2,VR1,VR0,WC1,WC0,WD1,WD0,VS,AC0,AR0,AW0,AD0,AS0
      COMMON/RIP/ST,NH,NAN1,NAN2    
      
C     /POT1/ is input
!$OMP THREADPRIVATE(/POT1/)

!$OMP THREADPRIVATE(/RIP/)

C     *******************************************************
C     END of privcom6
C     *******************************************************
