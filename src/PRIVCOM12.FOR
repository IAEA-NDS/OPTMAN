C     *******************************************************
C     START of privcom12  
C     *******************************************************
C==========================================================================
C     PRIVATE COMMONS FOR EACH PARALLEL THREAD TO BE USED AS INCLUDE BLOCK


C     FROM SUBROUTINE OVLAG
      COMMON/SHEM1/HW,AMG0,AMB0,GAM0,BET0,BET4,GAMDE,GSHAPE 

!$OMP THREADPRIVATE(/SHEM1/)

C     *******************************************************
C     END of privcom12
C     *******************************************************
